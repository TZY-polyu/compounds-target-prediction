"""
local_sea — 本地 SEA+TC 靶标预测引擎

基于 UCSF Shoichet Lab 的 Similarity Ensemble Approach (SEA) 算法，
使用 RDKit + ChEMBL 在本地计算靶标预测，无需网络依赖。

模块:
    data_extract: ChEMBL SQLite 数据提取与清洗
    fingerprints: 分子指纹计算 (ECFP4 2048-bit)
    calibration:  背景统计模型校准 (EVD 极值分布)
    predictor:    SEA+TC 预测引擎
"""

from .predictor import predict_targets

__version__ = "0.1.0"
