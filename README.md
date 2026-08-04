# Compounds Target Prediction

> **🌐 [中文版 (Chinese README)](./README_CN.md)**

Predict protein targets for small-molecule compounds using **SEA+TC** (Similarity Ensemble Approach with Tanimoto Coefficient enhancement) — the UCSF Shoichet Laboratory method.

---

## Key Features

| | |
|---|---|
| ⚡ **Fast** | ~170ms per query |
| 💾 **Lightweight** | Fingerprint DB only ~114 MB — runs on any ordinary computer |
| 🚫 **No training / no network** | Works on novel, virtual, or known compounds, fully offline |
| 🎯 **4,309 targets** | Powered by ChEMBL bioactivity data |
| 🔬 **SEA+TC dual threshold** | Improved accuracy over classic SEA (Irwin et al. 2018) |
| 📦 **Zero large downloads** | Extract the bundled 55 MB archive and start — no 29 GB ChEMBL database needed |

## Quick Start

```bash
# 1. Install dependencies
pip install rdkit numpy

# 2. Extract pre-built data (one-time)
tar -xzf compounds-targets-data.tar.gz

# 3. Single compound prediction
python scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" --pvalue 0.05 --top-n 5

# Output → result/<smiles>_<ts>.json
```

Results are saved to the `result/` directory in JSON format.

## Output Example

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

## Comparison with Web Servers

We compared this local engine against **SwissTargetPrediction (STP)** on 5 representative compounds (aspirin, caffeine, dopamine, dopamine analog, paracetamol), using the actual STP batch results captured during development (`data/comparison/swiss_results.csv`):

| Compound | STP targets | Local targets | Overlap | STP-side rate | Local-side rate |
|---|---|---|---|---|---|
| Aspirin | 100 | 70 | 30 | 30% | 43% |
| Caffeine | 100 | 40 | 19 | 19% | 48% |
| Dopamine analog | 100 | 100 | 54 | 54% | 54% |
| Dopamine | 100 | 51 | 28 | 28% | 55% |
| Paracetamol | 100 | 192 | 42 | 42% | 22% |

The two tools overlap on roughly **20–55%** of predictions. Both tools return distinct extra targets the other misses — differences come from database versions, scoring functions, and thresholds.

### Ground-truth check

| Drug | Known target | STP | Local |
|---|---|---|---|
| Aspirin | PTGS2 (CHEMBL230) | ✓ | ✓ |
| Aspirin | PTGS1 (CHEMBL221) | ✓ | ✗ |
| Caffeine | ADORA2A (CHEMBL251) | ✓ | ✓ |
| Caffeine | ADORA1 (CHEMBL253) | ✓ | ✗ |
| Dopamine | DRD2 (CHEMBL217) | ✓ | ✓ |
| Dopamine | DRD3 (CHEMBL234) | ✓ | ✓ |
| Paracetamol | PTGS2 (CHEMBL230) | ✗ | ✓ |

Both tools correctly identified most well-established drug targets, but **neither is complete** — each misses some known targets the other catches.

### About accuracy

There is **no single universal accuracy number** for target prediction — it depends heavily on the compound and the cutoff. Typical expectations:

| Scenario | Expected agreement |
|---|---|
| Top-5 hits vs. web servers | ~40–60% agreement |
| Well-known drug targets (e.g. aspirin→PTGS2) | Both tools usually hit |
| High-confidence calls (pvalue < 1e-5) | Most are biologically plausible |
| Entire prediction lists | Low overlap (different scoring) |

> **Honest expectation**: expect **~30–50% overlap** with web servers on full lists, and **~40–60% on top hits**. This is normal — different tools are complementary, not interchangeable. Use multiple tools and treat predictions as hypotheses.

## Setup

### 1. Install dependencies

```bash
pip install rdkit numpy
```

### 2. Extract pre-built data (recommended)

The pre-built fingerprint database is bundled as `compounds-targets-data.tar.gz` (55 MB) in this repo. Extract it and start predicting immediately — **no 5 GB ChEMBL download needed**:

```bash
tar -xzf compounds-targets-data.tar.gz
# → compounds-targets-data/target_fps.pkl (~114 MB)
# → compounds-targets-data/target_info.json (~3.5 MB)
# → compounds-targets-data/fit_params.json
```

> The `.pkl` and `target_info.json` files are not tracked directly in git (too large); they are shipped inside the archive. Do not re-download the 29 GB ChEMBL database unless you intend to rebuild the fingerprint database yourself.

### 3. Optional: Rebuild from ChEMBL

For users who want the latest database or custom filtering, download ChEMBL and rebuild (requires ~29 GB disk space, ~1-2 h calibration):

```bash
# Download ChEMBL SQLite (~5.4 GB compressed)
python scripts/download_chembl.py
# → Downloads, extracts, auto-generates target_info.json

# Then extract targets → fingerprints → calibrate
python scripts/local_sea/data_extract.py chembl.db \
  --output compounds-targets-data/target_ligands.pkl \
  --to-fingerprints
# → target_ligands.pkl + target_fps.pkl
```

→ Full guide: [SKILL.md](SKILL.md#regenerating-the-fingerprint-database)

## Data

| File | Size | Source |
|---|---|---|
| `compounds-targets-data.tar.gz` | 55 MB | **Bundled in repo** (extract & use) |
| `target_fps.pkl` | ~114 MB | Inside archive |
| `target_info.json` | ~3.5 MB | Inside archive |
| `fit_params.json` | <1 KB | In repo + archive |
| `chembl_*.db` | ~29 GB | Optional: `scripts/download_chembl.py` |
| `result/*.json` | varies | Prediction output |

**Default path (recommended)**: extract the bundled archive, no large downloads.

**Advanced path**: rebuild from the full ChEMBL database via the download script.

## Limitations

| Limitation | Detail |
|---|---|
| 🔬 **Requires experimental validation** | Predictions are statistical hypotheses, not confirmed binding. All results should be verified by biochemical assays. |
| 📚 **Bounded by ChEMBL coverage** | Only targets with known ligands in ChEMBL (~4,300 targets) can be predicted. Novel or poorly studied targets are invisible to this method. |
| 🧪 **Similarity bias** | SEA relies on chemical similarity. Compounds with entirely novel scaffolds may miss true targets that share no similar known ligands. |
| 🎯 **No affinity prediction** | The method predicts *whether* a compound may bind, not *how strongly* (no Kd/IC50). |
| ⚠️ **Promiscuous compounds** | Compounds like quercetin may return 200+ hits. High prediction count often indicates non-specific binding, not genuine polypharmacology. |
| 💾 **Large disk requirement** | ChEMBL SQLite requires ~29 GB after extraction. Fingerprint database adds ~114 MB. |

## Project Background

This tool was born out of practical necessity while screening large compound sets:

1. **SwissTargetPrediction first** — This site **bans IPs for rapid bulk access**. We first tried to work around it locally and got our local IP banned; after deploying to a server, running it directly without knowing this got the **server IP banned too**. Note: keeping requests ≥2s apart avoids the ban, but this dramatically increases runtime for batch screening.
2. **SEA web server next** — The new SEA website (sea.bkslab.org) has anti-crawling measures, and the old version is slow (~90 seconds per query), making batch screening impractical.
3. **Local deployment** — Downloaded the full ChEMBL SQLite database and deployed SEA+TC locally. No rate limits, ~170ms per query, fully reproducible.

> **Important**: This is a prediction tool. Results may differ from official web servers (SwissTargetPrediction, SEA, etc.) due to database versions, cutoffs, and algorithm parameters. **All predictions must be confirmed by experiments.**

## References

- Keiser et al. (2007) *Nature Biotechnology* 25(2), 197–206
- Keiser et al. (2009) *Nature* 462, 175–181
- Irwin et al. (2018) *J. Chem. Inf. Model.* 58(7)
- [SEA Wiki (Shoichet Lab)](https://wiki.docking.org/index.php/Category:SEA)

---

*Completed during my internship from June 30 to August 7, 2026. No salary — just a certificate. I uploaded this as a keepsake, not expecting anyone to see it, but it is my first real project after all, however modest. If you find any issues, feel free to reach out and I will fix them as soon as possible.*
