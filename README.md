# Compounds Target Prediction · 化合物靶标预测

Predict protein targets for small-molecule compounds using **SEA+TC** (Similarity Ensemble Approach with Tanimoto Coefficient enhancement) — the UCSF Shoichet Laboratory method.

基于 **SEA+TC** 算法预测小分子化合物可能结合的蛋白靶点 — UCSF Shoichet 实验室方法。

---

## What It Does · 功能

Given any valid SMILES string, this tool predicts which protein targets the compound is likely to bind. It compares the compound's ECFP4 fingerprint against 4,309 pre-built target-ligand fingerprint sets.

输入任意有效 SMILES，预测该化合物可能结合的蛋白靶点。引擎将化合物的 ECFP4 指纹与 4,309 个预构建靶点指纹集合进行统计检验。

| | |
|---|---|
| 🚫 **No training required** · 无需训练 | Works on novel, virtual, or known compounds |
| ⚡ **~170ms per query** · 单次查询约170ms | Local inference, no network needed |
| 🎯 **4,309 targets** · 覆盖4309个靶点 | Powered by ChEMBL bioactivity data |
| 🔬 **SEA+TC dual threshold** · 双阈值 | Improved accuracy (Irwin et al. 2018) |

## Quick Start · 快速开始

```bash
pip install rdkit numpy scipy

# Single compound · 单个化合物
python scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" --pvalue 0.05 --top-n 5
```

## Output Example · 输出示例

```json
{
  "target_key": "CHEMBL5847",
  "target_name": "AKR1C2",
  "gene_symbol": "AKR1C2",
  "description": "Aldo-keto reductase family 1 member C2 (Homo sapiens)",
  "pvalue": 8.29e-24,
  "maxtc": 0.4722,
  "probability": 0.7361
}
```

## Data · 数据

The pre-built fingerprint database (`target_fps.pkl`, ~114MB) is **not included** in this repository. To generate it, download ChEMBL SQLite and follow the rebuild guide.

预构建指纹数据库（`target_fps.pkl`，约114MB）**未包含**在本仓库。请下载 ChEMBL SQLite 后按重建指南操作。

→ [Rebuild Instructions · 重建指南](SKILL.md#regenerating-the-fingerprint-database)

## References · 参考文献

- Keiser et al. (2007) *Nature Biotechnology* 25(2), 197–206
- Keiser et al. (2009) *Nature* 462, 175–181
- Irwin et al. (2018) *J. Chem. Inf. Model.* 58(7)
- [SEA Wiki (Shoichet Lab)](https://wiki.docking.org/index.php/Category:SEA)
