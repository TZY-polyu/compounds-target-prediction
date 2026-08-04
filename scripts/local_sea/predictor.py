"""
SEA+TC 预测引擎

基于本地 ECFP4 指纹库和校准后的统计模型，预测化合物的蛋白靶标。

算法: SEA (Similarity Ensemble Approach) + TC (Tanimoto Coefficient) 增强
     SEA+TC = (P-value < cutoff) OR (MaxTc >= cutoff)
     参考: Irwin et al. (2018) J. Chem. Inf. Model.
"""

import math
from typing import Dict, List, Optional

from rdkit import DataStructs
from rdkit import RDLogger

from .fingerprints import compute_ecfp4
from .calibration import DEFAULT_FIT_PARAMS

RDLogger.logger().setLevel(RDLogger.ERROR)


def evd_pvalue(z: float) -> float:
    """
    Extreme Value Distribution (Gumbel) P-value 计算。

    公式来源: Keiser et al. (2007) Nature Biotechnology, Eq. 13-14

    Args:
        z: Z-score

    Returns:
        P-value (0~1)，越小越显著
    """
    if z <= 0:
        return 1.0

    gamma = 0.577215665  # 欧拉常数
    pi = math.pi

    x = -math.exp(-z * pi / (math.sqrt(6) - gamma))

    if z <= 28:
        # 直接计算
        p = 1.0 - math.exp(x)
    else:
        # 泰勒展开，避免 e^z 溢出（z > 28 时 e^z > 1.4e12）
        p = -x - (x ** 2) / 2.0 - (x ** 3) / 6.0

    # Cap at 1e-300 to avoid exact 0 (display/JSON issues)
    p = max(p, 1e-300)
    return min(p, 1.0)


def pvalue_to_probability(pvalue: float, maxtc: float) -> float:
    """
    将 P-value + MaxTC 转换为 0-1 概率分数。

    与现有 sea_target_predictor.py 保持一致的计算方式：
    - log P-value 映射到 0-1 (cap at -log10(1e-20) = 20)
    - 50% 权重 P-value + 50% 权重 MaxTC

    Args:
        pvalue: SEA P-value
        maxtc: 最大 Tanimoto 系数

    Returns:
        组合概率分数 (0~1)
    """
    if pvalue <= 0:
        return 1.0

    # P-value >= 1.0 → scaled P-value contribution = 0.0, but MaxTC is preserved
    if pvalue >= 1.0:
        scaled_log_p = 0.0
    else:
        log_p = -math.log10(pvalue)
        scaled_log_p = min(log_p / 20.0, 1.0)

    combined = 0.5 * scaled_log_p + 0.5 * maxtc
    return min(combined, 1.0)


def predict_targets(
    query_smiles: str,
    target_fps: Dict[str, List],
    fit_params: Optional[Dict] = None,
    pvalue_cutoff: float = 0.05,
    maxtc_cutoff: float = 0.4,
    top_n: Optional[int] = None,
) -> List[Dict]:
    """
    使用 SEA+TC 算法预测化合物的蛋白靶标。

    Args:
        query_smiles: 查询化合物的 SMILES
        target_fps: {target_id: [ECFP4_bitvector, ...]}
        fit_params: 背景模型参数 {"TS": 0.57, "mu": ..., "phi": ..., "eta": ...}
                    若为 None，使用 DEFAULT_FIT_PARAMS
        pvalue_cutoff: SEA P-value 阈值（默认 0.05）
        maxtc_cutoff: MaxTc 阈值（默认 0.4，用于 SEA+TC 增强）
        top_n: 返回前 N 个结果（None = 全部）

    Returns:
        预测靶标列表，按 P-value 升序排列
    """
    if fit_params is None:
        fit_params = dict(DEFAULT_FIT_PARAMS)

    # 1. 计算查询指纹
    query_fp = compute_ecfp4(query_smiles, nBits=2048)
    if query_fp is None:
        raise ValueError(f"Invalid SMILES: {query_smiles}")

    TS = fit_params.get("TS", 0.57)
    mu = fit_params.get("mu", DEFAULT_FIT_PARAMS["mu"])
    phi = fit_params.get("phi", DEFAULT_FIT_PARAMS["phi"])
    eta = fit_params.get("eta", DEFAULT_FIT_PARAMS["eta"])

    results = []

    for target_id, ligand_fps in target_fps.items():
        if not ligand_fps:
            continue

        # 2. 批量 Tanimoto 计算
        tcs = DataStructs.BulkTanimotoSimilarity(query_fp, ligand_fps)
        tcs_filtered = [tc for tc in tcs if tc < 0.9999]
        maxtc = max(tcs_filtered) if tcs_filtered else 0.0

        # 3. Raw Score
        rawscore = sum(tc for tc in tcs_filtered if tc >= TS)

        # 4. Z-score
        n_ligands = len(tcs_filtered)
        expected_raw = mu * n_ligands
        std_raw = phi * (n_ligands ** eta)

        if std_raw > 1e-10:
            z = (rawscore - expected_raw) / std_raw
        else:
            z = 0.0

        # 5. P-value (EVD)
        pvalue = evd_pvalue(z)

        # 6. SEA+TC 双条件判断
        if pvalue <= pvalue_cutoff or maxtc >= maxtc_cutoff:
            prob = pvalue_to_probability(pvalue, maxtc)
            results.append({
                "target_id": target_id,
                "pvalue": pvalue,
                "maxtc": round(maxtc, 4),
                "probability": round(prob, 4),
                "n_ligands": n_ligands,
                "rawscore": round(rawscore, 4),
                "z_score": round(z, 4) if math.isfinite(z) else 0.0,
            })

    # 按 Probability 降序，P-value 升序
    results.sort(key=lambda x: (-x["probability"], x["pvalue"]))

    if top_n is not None:
        results = results[:top_n]

    return results
