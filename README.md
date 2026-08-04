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

## Project Background · 项目背景

This tool was born out of practical necessity while screening large compound sets:

本项目源于实际筛选大量化合物时的需求：

1. **SwissTargetPrediction first · 最初使用 SwissTargetPrediction** — This site **bans IPs for rapid bulk access**. We first tried to work around it locally and got our local IP banned; after deploying to a server, running it directly without knowing this got the **server IP banned too**. Note: keeping requests ≥2s apart avoids the ban, but this dramatically increases runtime for batch screening. SwissTargetPrediction 会**因短时间内大量访问而封禁 IP**。最初在本地尝试绕过限制，本地 IP 被封；之后部署到服务器，因事先不知情直接运行，服务器 IP 也被封禁。若将访问间隔控制在 2 秒/次可避免封禁，但会大幅增加批量筛选的运行时间。
2. **SEA web server next · 改用 SEA 网站** — The new SEA website (sea.bkslab.org) has anti-crawling measures, and the old version is slow (~90 seconds per query), making batch screening impractical. SEA 新版网站有反爬措施，旧版查询约 90 秒一个，批量筛选不可行。
3. **Local deployment · 最终本地化部署** — Downloaded the full ChEMBL SQLite database and deployed SEA+TC locally. No rate limits, ~170ms per query, fully reproducible. 下载完整 ChEMBL 数据库，本地部署 SEA+TC。无速率限制，单次约 170ms，结果完全可复现。

> **Important · 重要**：This is a prediction tool. Results may differ from official web servers (SwissTargetPrediction, SEA, etc.) due to database versions, cutoffs, and algorithm parameters. **All predictions must be confirmed by experiments.** 这是预测工具，结果可能与官方网站（SwissTargetPrediction、SEA 等）不完全一致——数据库版本、阈值、算法参数不同都会导致差异。**一切预测结果以实验为准。**

## Setup · 环境配置

### 1. Install dependencies · 安装依赖

```bash
pip install rdkit numpy scipy
```

### 2. Extract pre-built data · 解压预构建数据（推荐，开箱即用）

The pre-built fingerprint database is bundled as `compounds-targets-data.tar.gz` (55 MB) in this repo. Extract it and start predicting immediately — **no 5 GB ChEMBL download needed**:

指纹数据库已随仓库打包在 `compounds-targets-data.tar.gz`（55 MB）。解压即可直接预测——**无需下载 5 GB 的 ChEMBL 数据库**：

```bash
tar -xzf compounds-targets-data.tar.gz
# → compounds-targets-data/target_fps.pkl (~114 MB)
# → compounds-targets-data/target_info.json (~3.5 MB)
# → compounds-targets-data/fit_params.json
```

> The `.pkl` and `target_info.json` files are not tracked directly in git (too large); they are shipped inside the archive. Do not re-download the 29 GB ChEMBL database unless you intend to rebuild the fingerprint database yourself. 原始 pkl 和 target_info.json 因体积较大不直接入库，而是打包在压缩包中。除非要自行重建指纹库，否则无需下载 29 GB 的 ChEMBL 数据库。

### 3. Optional: Rebuild from ChEMBL · 可选：从 ChEMBL 自行重建

For users who want the latest database or custom filtering, download ChEMBL and rebuild (requires ~29 GB disk space, ~1-2 h calibration):

需要最新数据或自定义筛选的用户，可下载 ChEMBL 自行重建（需约 29 GB 磁盘，校准约 1-2 小时）：

```bash
# Download ChEMBL SQLite (~5.4 GB compressed)
python scripts/download_chembl.py
# → Downloads, extracts, auto-generates target_info.json

# Then extract targets → fingerprints → calibrate
python scripts/local_sea/data_extract.py chembl.db --output compounds-targets-data/
python scripts/local_sea/fingerprints.py ...
```

→ Full guide · 完整指南：[SKILL.md](SKILL.md#regenerating-the-fingerprint-database)

## Quick Start · 快速开始

```bash
# Single compound · 单个化合物
python scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" --pvalue 0.05 --top-n 5

# Output → result/pred_<hash>_<ts>.json
```

Results are saved to the `result/` directory in JSON format.

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

| File · 文件 | Size · 大小 | Source · 来源 |
|---|---|---|
| `compounds-targets-data.tar.gz` | 55 MB | **Bundled in repo** · 随仓库提供（解压即用） |
| `target_fps.pkl` | ~114 MB | Inside archive · 在压缩包内 |
| `target_info.json` | ~3.5 MB | Inside archive · 在压缩包内 |
| `fit_params.json` | <1 KB | In repo + archive · 仓库及压缩包内 |
| `chembl_*.db` | ~29 GB | Optional · 可选：`scripts/download_chembl.py` |
| `result/*.json` | varies | Prediction output · 预测结果输出 |

**Default path (recommended) · 默认路径（推荐）**：extract the bundled archive, no large downloads. 解压随仓库的压缩包即可，无需大文件下载。

**Advanced path · 进阶路径**：rebuild from the full ChEMBL database via the download script. 通过下载脚本从完整 ChEMBL 数据库重建。

## Limitations · 局限性

| Limitation · 局限 | Detail · 说明 |
|---|---|
| 🔬 **Requires experimental validation** · 需要实验验证 | Predictions are statistical hypotheses, not confirmed binding. All results should be verified by biochemical assays. 预测结果是统计假设，非确认结合，需生化实验验证。 |
| 📚 **Bounded by ChEMBL coverage** · 受 ChEMBL 覆盖范围限制 | Only targets with known ligands in ChEMBL (~4,300 targets) can be predicted. Novel or poorly studied targets are invisible to this method. 只能预测 ChEMBL 中有已知配体的靶点。 |
| 🧪 **Similarity bias** · 相似性偏差 | SEA relies on chemical similarity. Compounds with entirely novel scaffolds may miss true targets that share no similar known ligands. 基于化学相似性，全新骨架可能漏掉真实靶点。 |
| 🎯 **No affinity prediction** · 无法预测亲和力 | The method predicts *whether* a compound may bind, not *how strongly* (no Kd/IC50). 只预测是否可能结合，不预测结合强度。 |
| ⚠️ **Promiscuous compounds** · 泛结合化合物 | Compounds like quercetin may return 200+ hits. High prediction count often indicates non-specific binding, not genuine polypharmacology. 泛结合化合物可能产生大量假阳性。 |
| 💾 **Large disk requirement** · 磁盘需求大 | ChEMBL SQLite requires ~29 GB after extraction. Fingerprint database adds ~114 MB. 数据库解压后约 29 GB，指纹文件约 114 MB。 |

## References · 参考文献

- Keiser et al. (2007) *Nature Biotechnology* 25(2), 197–206
- Keiser et al. (2009) *Nature* 462, 175–181
- Irwin et al. (2018) *J. Chem. Inf. Model.* 58(7)
- [SEA Wiki (Shoichet Lab)](https://wiki.docking.org/index.php/Category:SEA)
