# 数据说明

本文件夹 (`compounds-targets`) 存放 `compounds-target-prediction` skill 运行所需的预构建数据文件。该 skill 使用 **SEA+TC**（Similarity Ensemble Approach + Tanimoto Coefficient）算法，从小分子化合物的 SMILES 预测其可能结合的蛋白靶点。

---

## 文件清单

| 文件 | 大小 | 说明 |
|------|------|------|
| `target_fps.pkl` | ~114 MB | 预构建的靶点指纹数据库 |
| `target_info.json` | ~3.5 MB | 靶点元数据（基因符号、物种、蛋白描述等） |
| `fit_params.json` | 极小 | SEA 背景模型校准参数 |

---

## 各文件作用

### 1. `target_fps.pkl`

- **内容**：4315 个蛋白靶点的预计算 ECFP4 指纹集合（2048-bit Morgan fingerprint）。
- **来源**：基于 ChEMBL 生物活性数据，为每个靶点聚合其已知配体的指纹。
- **作用**：预测时，将输入化合物的 ECFP4 指纹与这些靶点指纹集合进行集合层面的相似性比较（SEA 统计检验）。
- **特点**：
  - 预构建完成，无需联网或重新训练。
  - 每个靶点对应一个配体指纹集合，用于计算与查询化合物的相似性。

### 2. `target_info.json`

- **内容**：靶点 ID 到人类可读注释的映射，字段包括：
  - `target_name`：基因符号（如 `PTGS2`、`MAOB`），优先从 ChEMBL `component_synonyms` 的 `GENE_SYMBOL` 提取；若无则退回蛋白名称或 ChEMBL ID。
  - `description`：完整蛋白名称与物种（如 `Prostaglandin G/H synthase 2 (Homo sapiens)`）。
  - `organism`：物种（如 `Homo sapiens`）。
  - `target_type`：靶点类型（如 `SINGLE PROTEIN`）。
- **作用**：为预测结果提供人类可读的基因符号和描述，替代原始的 ChEMBL ID。
- **注意**：元数据条目数（约 18,552）可能略多于指纹数据库中的实际靶点数，因为仅当靶点存在于指纹数据库中时才会被使用。

### 3. `fit_params.json`

- **内容**：SEA 背景模型的极值分布（Extreme Value Distribution, EVD）校准参数。
- **示例**：
  ```json
  {
    "TS": 0.35,
    "mu": 0.000207
  }
  ```
- **作用**：
  - `TS`：相似性阈值（Tanimoto  cutoff），仅当指纹相似度超过该阈值时纳入 SEA 统计。
  - `mu`：EVD 分布的位置参数，用于将原始相似性得分转换为 P-value。
- **特点**：一次性校准完成，日常预测无需重新计算。

---

## Top-N 筛选规则

预测流程中对结果的筛选与排序逻辑如下：

### 1. 阳性判定（SEA+TC 双阈值）

一个化合物-靶点对被判定为阳性预测，当满足以下任一条件：

```
P-value < cutoff   OR   MaxTc >= 0.4
```

- **P-value**：SEA 集合相似性的统计显著性，越小表示越显著。
- **MaxTc**：查询化合物与靶点已知配体之间的最大 Tanimoto 系数，取值 0–1，越大表示结构越相似。
- 默认参数：`--pvalue 1.0`，`--maxtc 0.4`。

### 2. 排序规则

通过阳性判定后，所有靶点按以下规则排序：

1. **probability 降序**（综合置信度，越大越靠前）。
2. **probability 相同时，按 pvalue 升序**（P-value 越小越靠前）。

### 3. Top-N 截取

- 默认 `--top-n 5`，仅返回排序后的前 5 个靶点。
- 设置 `--top-n 0` 返回全部阳性预测结果。
- 输出字段包括：`target_key`（ChEMBL ID）、`target_name`（基因符号）、`gene_symbol`（与 `target_name` 同值，用于下游 `target-intersect` 匹配）、`description`、`pvalue`、`maxtc`、`probability`。

### 4. 输出示例

```json
{
  "method": "sea+tc",
  "total_compounds": 1,
  "results": [
    {
      "smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "method": "sea+tc",
      "total_predictions": 70,
      "targets": [
        {
          "target_key": "CHEMBL5847",
          "target_name": "AKR1C2",
          "gene_symbol": "AKR1C2",
          "description": "Aldo-keto reductase family 1 member C2 (Homo sapiens)",
          "pvalue": 8.29e-24,
          "maxtc": 0.4722,
          "probability": 0.7361
        }
      ]
    }
  ]
}
```

---

## 使用方式

skill 入口脚本：

```bash
python /path/to/opencode_AISNPT/.opencode/skills/compounds-target-prediction/scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --pvalue 0.05 \
  --top-n 5
```

脚本默认通过相对路径 `../compounds-targets-data/` 读取本文件夹中的数据文件。
