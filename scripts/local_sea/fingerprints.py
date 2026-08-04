"""
分子指纹计算模块

使用 RDKit 计算 ECFP4 (Morgan, radius=2) 2048-bit 指纹，
支持批量计算、存储和相似度搜索。
"""

import pickle
import os
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)


def compute_ecfp4(
    smiles: str,
    radius: int = 2,
    nBits: int = 2048,
) -> Optional[DataStructs.ExplicitBitVect]:
    """
    计算 ECFP4 (Morgan) 指纹。

    ECFP4 = Extended Connectivity Fingerprint, diameter 4 (radius 2).
    这是 SEA 论文推荐的最优指纹类型。

    Args:
        smiles: 分子 SMILES
        radius: Morgan 指纹半径（默认 2，对应直径 4）
        nBits: 指纹长度（默认 2048）

    Returns:
        RDKit ExplicitBitVect，失败返回 None
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)


def compute_fp_matrix(
    target_ligands: Dict[str, List[Dict]],
    nBits: int = 2048,
) -> Tuple[Dict[str, List[DataStructs.ExplicitBitVect]], Dict]:
    """
    为所有靶标的所有配体计算指纹。

    Args:
        target_ligands: {target_id: [{"smiles": ..., ...}, ...]}
        nBits: 指纹长度

    Returns:
        (target_fps, stats)
        target_fps: {target_id: [ECFP4_bitvector, ...]}
        stats: 统计信息
    """
    target_fps = {}
    total_success = 0
    total_failed = 0

    for tid, ligands in target_ligands.items():
        fps = []
        for lig in ligands:
            fp = compute_ecfp4(lig["smiles"], nBits=nBits)
            if fp is not None:
                fps.append(fp)
                total_success += 1
            else:
                total_failed += 1
        if fps:
            target_fps[tid] = fps

    stats = {
        "total_success": total_success,
        "total_failed": total_failed,
        "n_targets": len(target_fps),
        "nBits": nBits,
    }

    return target_fps, stats


def save_fingerprints(target_fps: Dict, output_path: str):
    """保存指纹数据库为 pickle 文件"""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "wb") as f:
        pickle.dump(target_fps, f, protocol=pickle.HIGHEST_PROTOCOL)
    size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"Fingerprints saved to: {output_path} ({size_mb:.1f} MB)")


def load_fingerprints(input_path: str) -> Dict:
    """从 pickle 文件加载指纹数据库"""
    with open(input_path, "rb") as f:
        return pickle.load(f)
