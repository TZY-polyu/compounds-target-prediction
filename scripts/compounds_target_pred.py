#!/usr/bin/env python3
"""
Compounds Target Prediction Skill
==================================
Predict protein targets for small-molecule compounds using SEA+TC
(Similarity Ensemble Approach with Tanimoto Coefficient enhancement).

UCSF Shoichet Lab algorithm, implemented locally with RDKit + ChEMBL.
ECFP4 2048-bit fingerprints, Extreme Value Distribution P-values.
~170ms per query. No network required.
"""

import argparse
import json
import hashlib
import os
import sys
from datetime import datetime, timezone
from typing import List, Dict, Optional
from pathlib import Path


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def find_scripts_dir() -> Path:
    """Return the scripts/ directory."""
    return Path(__file__).resolve().parent


def get_cache_path(prefix: str, smiles: str, param: float, cache_dir: str) -> Path:
    h = hashlib.md5(f"{smiles}_{param}".encode()).hexdigest()[:12]
    return Path(cache_dir).resolve() / f"{prefix}_{h}.json"


# ── Backend ───────────────────────────────────────────────────────────────

_FPS = None
_FIT = None
_TARGET_INFO = None


def init(fp_path: str, fit_path: str, target_info_path: Optional[str] = None):
    """Load fingerprint database and calibration parameters."""
    global _FPS, _FIT, _TARGET_INFO

    from local_sea.fingerprints import load_fingerprints
    from local_sea.calibration import load_fit_params

    print(f"Loading fingerprint database: {fp_path}", file=sys.stderr)
    _FPS = load_fingerprints(fp_path)
    print(f"  {len(_FPS)} targets loaded", file=sys.stderr)

    print(f"Loading calibration params: {fit_path}", file=sys.stderr)
    _FIT = load_fit_params(fit_path)
    if _FIT:
        print(f"  TS={_FIT.get('TS', 'N/A')}, μ={_FIT.get('mu', 'N/A')}", file=sys.stderr)

    if target_info_path and os.path.exists(target_info_path):
        with open(target_info_path) as f:
            _TARGET_INFO = json.load(f)
        print(f"  Target info loaded: {len(_TARGET_INFO)} entries", file=sys.stderr)

    print(file=sys.stderr)


def predict_one(smiles: str, pvalue_cutoff: float, cache_dir: str,
                maxtc_cutoff: float = 0.4) -> List[Dict]:
    """Predict targets for a single SMILES string."""
    h = hashlib.md5(f"{smiles}_{pvalue_cutoff}_{maxtc_cutoff}".encode()).hexdigest()[:12]
    cache_path = Path(cache_dir).resolve() / f"sea_{h}.json"
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    from local_sea.predictor import predict_targets

    targets = predict_targets(
        smiles,
        _FPS,
        fit_params=_FIT,
        pvalue_cutoff=pvalue_cutoff,
        maxtc_cutoff=maxtc_cutoff,
    )

    # Merge target metadata
    if _TARGET_INFO:
        for t in targets:
            info = _TARGET_INFO.get(t["target_id"], {})
            t["target_name"] = info.get("target_name", t["target_id"])
            t["description"] = info.get("description", info.get("target_name", ""))
    else:
        for t in targets:
            t["target_name"] = t.get("target_name", t["target_id"])
            t["description"] = t.get("description", "")

    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(targets, f, indent=2, ensure_ascii=False)

    return targets


def normalize(t: Dict) -> Dict:
    """Normalize to standard output schema."""
    target_name = t.get("target_name", "")
    return {
        "target_key": t.get("target_id", ""),
        "target_name": target_name,
        "gene_symbol": target_name,
        "description": t.get("description", ""),
        "pvalue": t.get("pvalue", None),
        "maxtc": t.get("maxtc", None),
        "probability": t.get("probability", 0),
    }


# ── CLI ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compounds Target Prediction — local SEA+TC engine"
    )
    parser.add_argument("--smiles", help="Single SMILES string")
    parser.add_argument("--input-file", help="JSON file with a compounds array for batch mode")
    parser.add_argument("--pvalue", type=float, default=0.05,
                        help="P-value cutoff (default: 0.05)")
    parser.add_argument("--maxtc", type=float, default=0.4,
                        help="MaxTc cutoff for SEA+TC (default: 0.4)")
    parser.add_argument("--cache-dir", default=".cache",
                        help="Cache directory (default: .cache)")
    parser.add_argument("--output", help="Output JSON file path (default: result/pred_<hash>_<ts>.json or result/batch_<ts>.json)")
    parser.add_argument("--input-manifest", help="Path for input manifest JSON")
    parser.add_argument("--top-n", type=int, default=50,
                        help="Return top N predictions (default: 50, use 0 for all)")
    parser.add_argument("--fingerprints", default="../compounds-targets-data/target_fps.pkl",
                        help="Fingerprint database path")
    parser.add_argument("--fit-params", default="../compounds-targets-data/fit_params.json",
                        help="Calibration parameters path")
    parser.add_argument("--target-info", default="../compounds-targets-data/target_info.json",
                        help="Target metadata JSON (default: compounds-targets-data/target_info.json)")

    args = parser.parse_args()

    started_at = utc_now_iso()

    scripts_dir = find_scripts_dir()
    fp_path = (scripts_dir / args.fingerprints).resolve()
    fit_path = (scripts_dir / args.fit_params).resolve()
    target_info_path = (scripts_dir / args.target_info).resolve() if args.target_info else None

    init(str(fp_path), str(fit_path),
         str(target_info_path) if target_info_path else None)

    os.makedirs(args.cache_dir, exist_ok=True)

    # Collect SMILES
    smiles_list = []
    if args.smiles:
        smiles_list.append(args.smiles)
    elif args.input_file:
        with open(args.input_file) as f:
            data = json.load(f)
            for comp in data.get("compounds", data.get("results", [])):
                smi = comp.get("smiles", comp.get("precursor_smiles", ""))
                if smi:
                    smiles_list.append(smi)
    else:
        print("Error: provide --smiles or --input-file", file=sys.stderr)
        sys.exit(1)

    if args.input_manifest:
        input_manifest = {
            "started_at": started_at,
            "method": "sea+tc",
            "pvalue_cutoff": args.pvalue,
            "maxtc_cutoff": args.maxtc,
            "top_n": args.top_n,
            "compounds": [{"index": i, "smiles": s} for i, s in enumerate(smiles_list, start=1)],
        }
        with open(args.input_manifest, "w", encoding="utf-8") as f:
            json.dump(input_manifest, f, indent=2, ensure_ascii=False)
        print(f"Input manifest → {args.input_manifest}")

    all_results = []
    errors = []
    for smi in smiles_list:
        try:
            targets = predict_one(smi, args.pvalue, args.cache_dir, args.maxtc)
            targets.sort(key=lambda x: (-x.get("probability", 0), x.get("pvalue", 1)))
            total_count = len(targets)
            if args.top_n and args.top_n > 0:
                targets = targets[:args.top_n]

            all_results.append({
                "smiles": smi,
                "method": "sea+tc",
                "total_predictions": total_count,
                "targets": [normalize(t) for t in targets],
            })
        except Exception as e:
            errors.append({
                "smiles": smi,
                "error": str(e),
            })

    finished_at = utc_now_iso()
    try:
        t0 = datetime.fromisoformat(started_at)
        t1 = datetime.fromisoformat(finished_at)
        duration_seconds = round((t1 - t0).total_seconds(), 6)
    except Exception:
        duration_seconds = -1.0

    total_preds = sum(r.get("total_predictions", 0) for r in all_results)

    out = json.dumps({
        "method": "sea+tc",
        "total_compounds": len(smiles_list),
        "total_predictions": total_preds,
        "valid_count": len(all_results),
        "invalid_count": len(errors),
        "started_at": started_at,
        "finished_at": finished_at,
        "duration_seconds": duration_seconds,
        "results": all_results,
        "errors": errors,
    }, indent=2, ensure_ascii=False)

    output_path = args.output
    if not output_path:
        os.makedirs("result", exist_ok=True)
        if len(smiles_list) == 1 and args.smiles:
            smi_hash = hashlib.md5(smiles_list[0].encode()).hexdigest()[:8]
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = f"result/pred_{smi_hash}_{ts}.json"
        else:
            ts = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = f"result/batch_{ts}.json"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(out)
    print(f"Output → {output_path}")


if __name__ == "__main__":
    main()
