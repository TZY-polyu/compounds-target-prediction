# 化合物靶标预测

> **🌐 [English README](./README.md)**

基于 **SEA+TC**（相似性集成方法 + Tanimoto 系数增强）算法预测小分子化合物可能结合的蛋白靶点 — UCSF Shoichet 实验室方法。

---

## 核心优势

| | |
|---|---|
| ⚡ **快** | 单次查询约 170ms |
| 💾 **轻量** | 指纹库仅约 114 MB，任何普通电脑都能运行 |
| 🚫 **无需训练、无需联网** | 适用于全新/虚拟/已知化合物，完全离线 |
| 🎯 **覆盖 4,309 个靶点** | 基于 ChEMBL 生物活性数据 |
| 🔬 **SEA+TC 双阈值** | 相比经典 SEA 精度更高（Irwin et al. 2018） |
| 📦 **无需大文件下载** | 解压随仓库附带的 55 MB 数据包即可使用，无需下载 29 GB 数据库 |

## 快速开始

```bash
# 1. 安装依赖
pip install rdkit numpy scipy

# 2. 解压预构建数据（只需一次）
tar -xzf compounds-targets-data.tar.gz

# 3. 单个化合物预测
python scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" --pvalue 0.05 --top-n 5

# 输出 → result/<smiles>_<时间戳>.json
```

结果以 JSON 格式保存到 `result/` 目录。

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

## 与在线查询工具对比

我们使用开发期间实际抓取的 SwissTargetPrediction（STP）批量结果（`data/comparison/swiss_results.csv`），用 5 个代表性化合物（阿司匹林、咖啡因、多巴胺、多巴胺类似物、对乙酰氨基酚）与本工具做了真实对比：

| 化合物 | STP 靶点 | 本地靶点 | 重叠 | STP 侧重叠率 | 本地侧重叠率 |
|---|---|---|---|---|---|
| 阿司匹林 | 100 | 70 | 30 | 30% | 43% |
| 咖啡因 | 100 | 40 | 19 | 19% | 48% |
| 多巴胺类似物 | 100 | 100 | 54 | 54% | 54% |
| 多巴胺 | 100 | 51 | 28 | 28% | 55% |
| 对乙酰氨基酚 | 100 | 192 | 42 | 42% | 22% |

两个工具的重叠率约为 **20–55%**，各自都有对方未命中的额外靶点——差异源于数据库版本、打分函数和阈值不同。

### 已知药靶验证

| 药物 | 已知靶点 | STP | 本地 |
|---|---|---|---|
| 阿司匹林 | PTGS2 (CHEMBL230) | ✓ | ✓ |
| 阿司匹林 | PTGS1 (CHEMBL221) | ✓ | ✗ |
| 咖啡因 | ADORA2A (CHEMBL251) | ✓ | ✓ |
| 咖啡因 | ADORA1 (CHEMBL253) | ✓ | ✗ |
| 多巴胺 | DRD2 (CHEMBL217) | ✓ | ✓ |
| 多巴胺 | DRD3 (CHEMBL234) | ✓ | ✓ |
| 对乙酰氨基酚 | PTGS2 (CHEMBL230) | ✗ | ✓ |

两个工具都能命中大部分公认的药靶，但**都不完整**——各有遗漏。

### 关于准确率

靶点预测**没有统一的准确率数字**——它高度依赖化合物本身和阈值选择。参考经验值：

| 场景 | 预期一致率 |
|---|---|
| 与在线工具 top-5 对比 | ~40–60% |
| 公认药靶（如 阿司匹林→PTGS2） | 两工具通常都能命中 |
| 高置信预测（pvalue < 1e-5） | 大多有生物学合理性 |
| 完整预测列表 | 重叠率低（打分机制不同） |

> **诚实预期**：与在线工具完整列表重叠率约 **30–50%**，top 命中约 **40–60%**。这是正常现象——不同工具互补而非互替。建议多工具交叉验证，预测仅作假设。

## 环境配置

### 1. 安装依赖

```bash
pip install rdkit numpy scipy
```

### 2. 解压预构建数据（推荐，开箱即用）

指纹数据库已随仓库打包在 `compounds-targets-data.tar.gz`（55 MB）。解压即可直接预测——**无需下载 5 GB 的 ChEMBL 数据库**：

```bash
tar -xzf compounds-targets-data.tar.gz
# → compounds-targets-data/target_fps.pkl（约 114 MB）
# → compounds-targets-data/target_info.json（约 3.5 MB）
# → compounds-targets-data/fit_params.json
```

> 原始 pkl 和 target_info.json 因体积较大不直接入库，而是打包在压缩包中。除非要自行重建指纹库，否则无需下载 29 GB 的 ChEMBL 数据库。

### 3. 可选：从 ChEMBL 自行重建

需要最新数据或自定义筛选的用户，可下载 ChEMBL 自行重建（需约 29 GB 磁盘，校准约 1-2 小时）：

```bash
# 下载 ChEMBL SQLite（压缩包约 5.4 GB）
python scripts/download_chembl.py
# → 下载、解压、自动生成 target_info.json

# 然后提取靶点 → 计算指纹 → 校准
python scripts/local_sea/data_extract.py chembl.db --output compounds-targets-data/
python scripts/local_sea/fingerprints.py ...
```

→ 完整指南：[SKILL.md](SKILL.md#regenerating-the-fingerprint-database)

## 数据

| 文件 | 大小 | 来源 |
|---|---|---|
| `compounds-targets-data.tar.gz` | 55 MB | **随仓库提供**（解压即用） |
| `target_fps.pkl` | ~114 MB | 在压缩包内 |
| `target_info.json` | ~3.5 MB | 在压缩包内 |
| `fit_params.json` | <1 KB | 仓库及压缩包内 |
| `chembl_*.db` | ~29 GB | 可选：`scripts/download_chembl.py` |
| `result/*.json` | 不定 | 预测结果输出 |

**默认路径（推荐）**：解压随仓库的压缩包即可，无需大文件下载。

**进阶路径**：通过下载脚本从完整 ChEMBL 数据库重建。

## 局限性

| 局限 | 说明 |
|---|---|
| 🔬 **需要实验验证** | 预测结果是统计假设，非确认结合，需生化实验验证。 |
| 📚 **受 ChEMBL 覆盖范围限制** | 只能预测 ChEMBL 中有已知配体的靶点（约 4,300 个）。全新或研究较少的靶点无法预测。 |
| 🧪 **相似性偏差** | 基于化学相似性，全新骨架可能漏掉真实靶点。 |
| 🎯 **无法预测亲和力** | 只预测是否可能结合，不预测结合强度（无 Kd/IC50）。 |
| ⚠️ **泛结合化合物** | 泛结合化合物可能产生大量假阳性。 |
| 💾 **磁盘需求大** | ChEMBL 数据库解压后约 29 GB，指纹文件约 114 MB。 |

## 项目背景

本项目源于实际筛选大量化合物时的需求：

1. **最初使用 SwissTargetPrediction** — 该网站会**因短时间内大量访问而封禁 IP**。最初在本地尝试绕过限制，本地 IP 被封；之后部署到服务器，因事先不知情直接运行，服务器 IP 也被封禁。若将访问间隔控制在 2 秒/次可避免封禁，但会大幅增加批量筛选的运行时间。
2. **改用 SEA 网站** — SEA 新版网站（sea.bkslab.org）有反爬措施，旧版查询约 90 秒一个，批量筛选不可行。
3. **最终本地化部署** — 下载完整 ChEMBL 数据库，本地部署 SEA+TC。无速率限制，单次约 170ms，结果完全可复现。

> **重要**：这是预测工具，结果可能与官方网站（SwissTargetPrediction、SEA 等）不完全一致——数据库版本、阈值、算法参数不同都会导致差异。**一切预测结果以实验为准。**

## 参考文献

- Keiser et al. (2007) *Nature Biotechnology* 25(2), 197–206
- Keiser et al. (2009) *Nature* 462, 175–181
- Irwin et al. (2018) *J. Chem. Inf. Model.* 58(7)
- [SEA Wiki (Shoichet Lab)](https://wiki.docking.org/index.php/Category:SEA)

---

2026.6.30-2026.8.7实习期间完成，没工资，只有一张证明，想着留一个纪念从而上传，没指望会有人看到，但毕竟这是我本人第一个正儿八经的项目，虽然没啥含金量。如果项目有问题，可以及时联系，我会尽快修改
