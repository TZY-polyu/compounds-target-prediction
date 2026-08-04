#!/usr/bin/env python3
"""
data_extract — 从 ChEMBL SQLite 提取靶点-配体数据
====================================================
从 ChEMBL SQLite 数据库中提取"靶点 → 已知配体 SMILES 集合"，
输出供 fingerprints.py 使用（compute_fp_matrix 输入格式）。

筛选逻辑（复现 compounds-target-prediction 预构建指纹库）:
    1. 只保留 target_type = 'SINGLE PROTEIN'（单一蛋白靶点）
    2. 只保留有活性的配体记录: pchembl_value >= --pchembl（默认 5.0，即 IC50 <= 10 uM）
    3. 同一靶点的配体 SMILES 去重
    4. 靶点至少要有 --min-ligands 个配体（默认 1）

用法:
    python scripts/local_sea/data_extract.py chembl_34.db --output target_ligands.pkl
    python scripts/local_sea/data_extract.py chembl_34.db --output target_ligands.json --human-only
"""

import argparse
import json
import pickle
import sqlite3
import sys
from pathlib import Path

try:
    from .fingerprints import compute_fp_matrix, save_fingerprints
except ImportError:
    from fingerprints import compute_fp_matrix, save_fingerprints


EXTRACT_SQL = """
SELECT
    td.chembl_id                       AS target_id,
    td.target_type                     AS target_type,
    cs.canonical_smiles                AS smiles,
    a.pchembl_value                    AS pchembl_value
FROM activities a
JOIN assays ass        ON a.assay_id = ass.assay_id
JOIN target_dictionary td ON ass.tid = td.tid
JOIN compound_structures cs ON a.molregno = cs.molregno
WHERE td.target_type = 'SINGLE PROTEIN'
  AND cs.canonical_smiles IS NOT NULL
  AND a.pchembl_value IS NOT NULL
"""


def extract_target_ligands(
    db_path: str,
    pchembl_cutoff: float = 5.0,
    min_ligands: int = 1,
    human_only: bool = False,
) -> dict:
    """从 ChEMBL SQLite 提取 {target_id: [{"smiles": ...}, ...]}."""
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    sql = EXTRACT_SQL
    if human_only:
        # 只保留 Homo sapiens 相关: 在 component_synonyms 中标记的物种
        sql = """
        SELECT DISTINCT
            td.chembl_id AS target_id,
            td.target_type AS target_type,
            cs.canonical_smiles AS smiles,
            a.pchembl_value AS pchembl_value
        FROM activities a
        JOIN assays ass ON a.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        JOIN compound_structures cs ON a.molregno = cs.molregno
        JOIN target_components tc ON tc.tid = td.tid
        JOIN component_synonyms comp ON comp.component_id = tc.component_id
        WHERE td.target_type = 'SINGLE PROTEIN'
          AND cs.canonical_smiles IS NOT NULL
          AND a.pchembl_value IS NOT NULL
          AND comp.syn_type = 'ORGANISM'
          AND UPPER(comp.component_synonym) LIKE '%HOMO SAPIENS%'
        """

    print(f"Querying {db_path} ...", file=sys.stderr)
    print(f"  pchembl_cutoff = {pchembl_cutoff}, min_ligands = {min_ligands}, "
          f"human_only = {human_only}", file=sys.stderr)

    cur.execute(sql)

    # 聚合: target_id -> 有序 SMILES 集合（按 pchembl_value 降序，保证活性高的优先）
    targets: dict = {}
    for row in cur.fetchall():
        pv = row["pchembl_value"]
        if pv is None or pv < pchembl_cutoff:
            continue
        tid = row["target_id"]
        smi = row["smiles"]
        if tid not in targets:
            targets[tid] = {"ligands": {}, "n_high_confidence": 0}
        # 用 SMILES 作 key 去重; 记录最高 pchembl_value
        prev = targets[tid]["ligands"].get(smi)
        if prev is None or pv > prev:
            targets[tid]["ligands"][smi] = pv
            if pv >= 6.0:  # IC50 <= 1uM 视为高置信
                targets[tid]["n_high_confidence"] += 1

    conn.close()

    # 组装输出结构, 过滤 min_ligands
    result: dict = {}
    for tid, data in targets.items():
        if len(data["ligands"]) < min_ligands:
            continue
        result[tid] = [{"smiles": s} for s in sorted(
            data["ligands"], key=lambda s: -data["ligands"][s]
        )]

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Extract target-ligand data from ChEMBL SQLite"
    )
    parser.add_argument("db", help="Path to chembl_*.db SQLite file")
    parser.add_argument("--output", "-o", default="target_ligands.pkl",
                        help="Output file (.pkl or .json, default: target_ligands.pkl)")
    parser.add_argument("--pchembl", type=float, default=5.0,
                        help="pchembl_value cutoff (default: 5.0, i.e. IC50<=10uM)")
    parser.add_argument("--min-ligands", type=int, default=1,
                        help="Min ligands per target (default: 1)")
    parser.add_argument("--human-only", action="store_true",
                        help="Only keep Homo sapiens targets")
    parser.add_argument("--to-fingerprints", action="store_true",
                        help="Also compute ECFP4 fingerprints and save as target_fps.pkl")
    args = parser.parse_args()

    if not Path(args.db).exists():
        print(f"Error: database not found: {args.db}", file=sys.stderr)
        sys.exit(1)

    result = extract_target_ligands(
        args.db,
        pchembl_cutoff=args.pchembl,
        min_ligands=args.min_ligands,
        human_only=args.human_only,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix == ".json":
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=False)
    else:
        with open(output_path, "wb") as f:
            pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)

    n_targets = len(result)
    n_ligands = sum(len(v) for v in result.values())
    print(f"\nExtracted {n_targets} targets, {n_ligands} ligands total")
    print(f"Saved → {output_path}")

    if args.to_fingerprints:
        print("\nComputing ECFP4 fingerprints ...")
        target_fps, stats = compute_fp_matrix(result)
        print(f"  {stats['n_targets']} targets, "
              f"{stats['total_success']} fingerprints computed, "
              f"{stats['total_failed']} failed")
        fp_path = output_path.parent / "target_fps.pkl"
        save_fingerprints(target_fps, str(fp_path))


if __name__ == "__main__":
    main()
