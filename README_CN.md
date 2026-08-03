# 化合物靶标预测

基于 **SEA+TC**（相似性集成方法 + Tanimoto 系数增强）算法，预测小分子化合物可能结合的蛋白靶点 — UCSF Shoichet 实验室方法。

## 功能

输入任意有效 SMILES 字符串，预测该化合物可能结合的蛋白靶点。引擎计算化合物的 ECFP4 指纹，与 4,309 个预构建的靶点-配体指纹集合进行极值分布统计检验。

- **无需训练** — 适用于已知、虚拟或全新化合物
- **单次查询约 170ms** — 本地推理，无需联网
- **覆盖 4,309 个靶点** — 基于 ChEMBL 生物活性数据
- **SEA+TC 双阈值** — 相比传统 SEA 精度更高 (Irwin et al. 2018)

## 快速开始

```bash
pip install rdkit numpy scipy

# 单个化合物预测
python scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" --pvalue 0.05 --top-n 5
```

## 输出示例

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

## 数据

预构建的指纹数据库（`target_fps.pkl`，约 114MB）**未包含**在本仓库中。生成方法：

1. 从 [ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/) 下载 ChEMBL SQLite
2. 按照 [SKILL.md](SKILL.md#regenerating-the-fingerprint-database) 中的重建步骤操作

## 参考文献

- Keiser et al. (2007) *Nature Biotechnology* 25(2), 197–206
- Keiser et al. (2009) *Nature* 462, 175–181
- Irwin et al. (2018) *J. Chem. Inf. Model.* 58(7)
- [SEA Wiki (Shoichet Lab)](https://wiki.docking.org/index.php/Category:SEA)
