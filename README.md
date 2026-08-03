# Compounds Target Prediction

Predict protein targets for small-molecule compounds using **SEA+TC** (Similarity Ensemble Approach with Tanimoto Coefficient enhancement) — the UCSF Shoichet Laboratory method.

## What It Does

Given any valid SMILES string, this tool predicts which protein targets the compound is likely to bind. It compares the compound's ECFP4 fingerprint against 4,309 pre-built target-ligand fingerprint sets using a statistical test based on extreme value distribution.

- **No training required** — works on novel, virtual, or known compounds
- **~170ms per query** — local inference, no network needed
- **4,309 targets** — powered by ChEMBL bioactivity data
- **SEA+TC dual threshold** — improved accuracy over traditional SEA (Irwin et al. 2018)

## Quick Start

```bash
pip install rdkit numpy scipy

# Single compound
python scripts/compounds_target_pred.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" --pvalue 0.05 --top-n 5
```

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

## Data

The pre-built fingerprint database (`target_fps.pkl`, ~114MB) is **not included** in this repository. To generate it:

1. Download ChEMBL SQLite from [ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/)
2. Follow the rebuild instructions in [SKILL.md](SKILL.md#regenerating-the-fingerprint-database)

## References

- Keiser et al. (2007) *Nature Biotechnology* 25(2), 197–206
- Keiser et al. (2009) *Nature* 462, 175–181
- Irwin et al. (2018) *J. Chem. Inf. Model.* 58(7)
- [SEA Wiki (Shoichet Lab)](https://wiki.docking.org/index.php/Category:SEA)
