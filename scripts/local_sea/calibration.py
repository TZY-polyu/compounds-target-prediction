"""
背景统计模型校准模块

基于随机配体集合配对，校准 SEA 的统计参数：
- TS: 最优 Tanimoto 阈值
- μ: 期望 Raw Score 斜率 (F_mean = μ·s)
- φ, η: 标准差参数 (F_sd = φ·s^η)

算法来源: Keiser et al. (2007); Wang et al. (2016) Algorithm 1
"""

import json
import os
import random
import math
from typing import Dict, List, Optional

import numpy as np
from rdkit import DataStructs
from rdkit import RDLogger


RDLogger.logger().setLevel(RDLogger.ERROR)


# 论文中的经验默认值（跳过完整校准时的回退方案）
# 来源: Keiser et al. 2007; Hert et al. 2008
DEFAULT_FIT_PARAMS = {
    "TS": 0.57,
    "mu": 0.0123,
    "phi": 0.0456,
    "eta": 0.89,
}


def _gumbel_pdf(z: np.ndarray) -> np.ndarray:
    """Gumbel (极值 I 型) 分布 PDF"""
    return np.exp(-z - np.exp(-z))


def calibrate_background(
    all_fps: List[DataStructs.ExplicitBitVect],
    ts_range: Optional[List[float]] = None,
    set_size_range: Optional[List[int]] = None,
    n_repeats: int = 30,
    seed: int = 42,
    verbose: bool = True,
) -> Dict:
    """
    背景模型校准。

    从所有 ChEMBL 配体指纹中随机抽样，模拟随机配体集合之间的
    相似度分布，拟合极值分布参数。

    Args:
        all_fps: 所有 ChEMBL 配体的指纹列表（作为随机抽样池）
        ts_range: TS 搜索范围（默认 0.00 ~ 0.99, 步长 0.01）
        set_size_range: 集合大小范围（默认 10, 20, ..., 60, 100, ..., 600, 1000）
        n_repeats: 每个 TS × 每个集合大小对的重复次数
        seed: 随机种子
        verbose: 输出详细进度

    Returns:
        fit_params: {"TS": float, "mu": float, "phi": float, "eta": float}
    """
    random.seed(seed)
    np.random.seed(seed)

    if ts_range is None:
        ts_range = [round(i * 0.01, 2) for i in range(0, 100)]
    if set_size_range is None:
        # 小集合密集采样，大集合稀疏采样
        set_size_range = list(range(10, 70, 10)) + list(range(100, 700, 100)) + [1000]

    pool_size = len(all_fps)
    if verbose:
        print("Calibrating background model...")
        print(f"  Pool size: {pool_size} fingerprints")
        print(f"  TS range: {ts_range[0]:.2f} ~ {ts_range[-1]:.2f} ({len(ts_range)} steps)")
        print(f"  Set sizes: {len(set_size_range)} ({set_size_range[0]} ~ {set_size_range[-1]})")
        print(f"  Repeats per config: {n_repeats}")
        total = len(ts_range) * len(set_size_range) ** 2 * n_repeats
        print(f"  Total comparisons: ~{total:.0e}")
        print()

    # Phase 1: 收集所有 (size_product, rawscore) 数据
    # 使用固定的随机种子确保可重现
    all_data = []  # [(size_product, rawscore), ...]

    # 为每个 TS 值预先收集
    for ts_idx, ts in enumerate(ts_range):
        if verbose and ts_idx % 10 == 0:
            print(f"  TS={ts:.2f} ({ts_idx + 1}/{len(ts_range)})...")

        for size_a in set_size_range:
            for size_b in set_size_range:
                s_product = size_a * size_b
                actual_a = min(size_a, pool_size)
                actual_b = min(size_b, pool_size)

                for _ in range(n_repeats):
                    # 随机抽样
                    indices_a = random.sample(range(pool_size), actual_a)
                    indices_b = random.sample(range(pool_size), actual_b)

                    A = [all_fps[i] for i in indices_a]
                    B = [all_fps[i] for i in indices_b]

                    # 计算 Raw Score
                    rawscore = 0.0
                    for fp_a in A:
                        for fp_b in B:
                            tc = DataStructs.TanimotoSimilarity(fp_a, fp_b)
                            if tc >= ts:
                                rawscore += tc

                    all_data.append((s_product, rawscore, ts))

    if verbose:
        print()

    # Phase 2: 对每个 TS，拟合参数并评估拟合优度
    best_chi2 = float("inf")
    best_fit = None
    best_ts = None

    unique_ts = sorted(set(d[2] for d in all_data))

    for ts in unique_ts:
        ts_data = [(s, r) for s, r, t in all_data if t == ts]
        sizes = np.array([d[0] for d in ts_data])
        raws = np.array([d[1] for d in ts_data])

        # 按 size_product 分组统计
        unique_sizes = sorted(set(sizes))
        mean_by_size = np.array([np.mean(raws[sizes == s]) for s in unique_sizes])
        std_by_size = np.array([np.std(raws[sizes == s], ddof=1) for s in unique_sizes])
        sizes_u = np.array(unique_sizes)

        # 过滤掉 std=0 的组
        valid = std_by_size > 0
        if valid.sum() < 5:
            continue

        sizes_valid = sizes_u[valid]
        means_valid = mean_by_size[valid]
        stds_valid = std_by_size[valid]

        # 线性拟合 F_mean(s) = μ · s（过原点）
        mu = np.sum(sizes_valid * means_valid) / np.sum(sizes_valid ** 2)

        # 幂律拟合 F_sd(s) = φ · s^η（对数空间线性回归）
        log_sizes = np.log(sizes_valid)
        log_stds = np.log(stds_valid)
        slope, intercept = np.polyfit(log_sizes, log_stds, 1)
        phi = math.exp(intercept)
        eta = slope

        # 计算所有 Z-score 的卡方拟合度
        all_z = (raws - mu * sizes) / (phi * sizes ** eta + 1e-10)
        all_z = all_z[np.isfinite(all_z)]

        if len(all_z) < 100:
            continue

        hist, edges = np.histogram(all_z, bins=50, density=True)
        bin_centers = (edges[:-1] + edges[1:]) / 2

        expected = _gumbel_pdf(bin_centers)
        expected = expected / expected.sum() * hist.sum()

        chi2_stat = np.sum((hist - expected) ** 2 / (expected + 1e-10))

        if chi2_stat < best_chi2:
            best_chi2 = chi2_stat
            best_fit = {"mu": round(mu, 6), "phi": round(phi, 6), "eta": round(eta, 4)}
            best_ts = ts

    if best_fit is None or best_ts is None:
        if verbose:
            print("WARNING: Calibration failed, using default parameters")
        return dict(DEFAULT_FIT_PARAMS)

    best_fit["TS"] = round(best_ts, 2)
    best_fit["chi2"] = round(float(best_chi2), 4)

    if verbose:
        print("Best fit:")
        print(f"  TS  = {best_fit['TS']}")
        print(f"  μ   = {best_fit['mu']}")
        print(f"  φ   = {best_fit['phi']}")
        print(f"  η   = {best_fit['eta']}")
        print(f"  χ²  = {best_fit['chi2']}")

    return best_fit


def quick_calibrate(
    target_fps: Dict[str, List[DataStructs.ExplicitBitVect]],
    sample_size: int = 50000,
    n_repeats: int = 10,
    seed: int = 42,
    verbose: bool = True,
) -> Dict:
    """
    快速校准（用于开发测试，~10 分钟）。

    从 target_fps 中随机抽 sample_size 个指纹作为背景池，
    使用稀疏的参数网格进行校准。
    """
    # 构建指纹池
    all_fps = []
    for fps in target_fps.values():
        all_fps.extend(fps)

    if len(all_fps) > sample_size:
        random.seed(seed)
        all_fps = random.sample(all_fps, sample_size)

    # 稀疏 TS 网格
    ts_range = [round(i * 0.05, 2) for i in range(0, 20)]

    # 稀疏集合大小网格
    set_size_range = [10, 30, 50, 100, 200, 500]

    return calibrate_background(
        all_fps,
        ts_range=ts_range,
        set_size_range=set_size_range,
        n_repeats=n_repeats,
        seed=seed,
        verbose=verbose,
    )


def save_fit_params(fit_params: Dict, output_path: str):
    """保存校准参数为 JSON"""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(fit_params, f, indent=2, ensure_ascii=False)
    print(f"Fit params saved to: {output_path}")


def load_fit_params(input_path: str) -> Dict:
    """从 JSON 加载校准参数"""
    if not os.path.exists(input_path):
        print(f"Fit params not found: {input_path}, using defaults")
        return dict(DEFAULT_FIT_PARAMS)
    with open(input_path, "r", encoding="utf-8") as f:
        return json.load(f)
