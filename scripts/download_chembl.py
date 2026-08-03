#!/usr/bin/env python3
"""
ChEMBL Database Downloader + Target Info Generator
===================================================
Automatically download ChEMBL SQLite, extract, and generate target_info.json.

Usage:
    python scripts/download_chembl.py
    python scripts/download_chembl.py --output /path/to/dir
    python scripts/download_chembl.py --yes   # skip confirmation

The script will:
    1. Check if chembl_*.db already exists in the output directory
    2. Ask for confirmation before downloading (~5.4GB compressed, ~29GB extracted)
    3. Download from the official ChEMBL FTP mirror
    4. Extract the .tar.gz archive
    5. Auto-generate target_info.json with gene symbols from ChEMBL
"""

import argparse
import json
import os
import sqlite3
import sys
import tarfile
import urllib.request
from pathlib import Path
from urllib.parse import urlparse


# ── Constants ──────────────────────────────────────────────────────────────

CHEMBL_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
CHEMBL_DB_NAME = "chembl_34_sqlite.tar.gz"  # Update when new version released
EXPECTED_SIZE_GB = 5.4


# ── Helpers ────────────────────────────────────────────────────────────────

def format_size(bytes_val: int) -> str:
    """Format bytes to human-readable string."""
    for unit in ("B", "KB", "MB", "GB"):
        if bytes_val < 1024:
            return f"{bytes_val:.1f} {unit}"
        bytes_val /= 1024
    return f"{bytes_val:.1f} TB"


def progress_hook(block_num: int, block_size: int, total_size: int):
    """Show download progress bar."""
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(downloaded / total_size * 100, 100)
        bar_len = 30
        filled = int(bar_len * pct / 100)
        bar = "█" * filled + "░" * (bar_len - filled)
        print(
            f"\r  [{bar}] {pct:5.1f}%  "
            f"{format_size(downloaded)} / {format_size(total_size)}",
            end="",
            flush=True,
        )
    else:
        print(f"\r  Downloaded: {format_size(downloaded)}", end="", flush=True)


def find_existing_db(output_dir: Path) -> list[Path]:
    """Find existing ChEMBL SQLite files in the output directory."""
    return sorted(output_dir.glob("chembl_*.db"))


def find_latest_db(output_dir: Path) -> Path | None:
    """Return the most recent ChEMBL SQLite file, or None."""
    dbs = find_existing_db(output_dir)
    return dbs[-1] if dbs else None


def extract_tar_gz(archive_path: Path, output_dir: Path) -> Path:
    """Extract .tar.gz and return the .db file path."""
    print(f"\nExtracting {archive_path.name} ...")
    with tarfile.open(archive_path, "r:gz") as tar:
        members = tar.getmembers()
        db_member = None
        for m in members:
            if m.name.endswith(".db"):
                db_member = m
                break

        if db_member:
            size_gb = db_member.size / (1024**3)
            print(f"  Database file: {db_member.name} ({size_gb:.1f} GB)")
            print("  Extracting (this may take a few minutes) ...")
            tar.extract(db_member, path=output_dir, filter="data")
            db_path = output_dir / db_member.name
            print(f"  Extracted → {db_path}")
            return db_path
        else:
            print("  Extracting all files ...")
            tar.extractall(path=output_dir, filter="data")
            db_files = sorted(output_dir.glob("chembl_*.db"))
            if db_files:
                print(f"  Extracted → {db_files[0]}")
                return db_files[0]
            else:
                raise RuntimeError("No .db file found in archive")


# ── Target Info Generator ──────────────────────────────────────────────────

def generate_target_info(db_path: Path, output_dir: Path) -> Path:
    """
    Generate target_info.json from ChEMBL SQLite.
    Extracts gene symbols via target_components + component_synonyms.
    """
    info_path = output_dir / "target_info.json"

    print(f"\nGenerating target metadata from {db_path.name} ...")
    conn = sqlite3.connect(str(db_path))
    cur = conn.cursor()

    # Step 1: Get all targets
    cur.execute(
        "SELECT chembl_id, pref_name, organism, target_type, tid "
        "FROM target_dictionary"
    )
    targets = cur.fetchall()
    print(f"  {len(targets)} targets found in target_dictionary")

    # Step 2: For each target, find GENE_SYMBOL from component_synonyms
    info = {}
    for i, (chembl_id, pref_name, organism, target_type, tid) in enumerate(targets):
        cur.execute(
            "SELECT DISTINCT cs.component_synonym "
            "FROM target_components tc "
            "JOIN component_synonyms cs ON tc.component_id = cs.component_id "
            "WHERE tc.tid = ? AND cs.syn_type = 'GENE_SYMBOL'",
            (tid,),
        )
        symbols = [r[0] for r in cur.fetchall()]

        gene_str = ";".join(symbols) if symbols else None
        target_name = gene_str if gene_str else (pref_name or chembl_id)
        description = (
            f"{pref_name} ({organism})"
            if pref_name and organism
            else (pref_name or gene_str or chembl_id)
        )

        info[chembl_id] = {
            "target_name": target_name,
            "description": description,
            "organism": organism or "",
            "target_type": target_type or "",
        }

        if (i + 1) % 5000 == 0:
            print(f"  Progress: {i + 1}/{len(targets)} targets processed")

    conn.close()

    with open(info_path, "w", encoding="utf-8") as f:
        json.dump(info, f, ensure_ascii=False, indent=2)

    file_size = os.path.getsize(info_path)
    print(f"  Saved {len(info)} entries → {info_path} ({format_size(file_size)})")
    return info_path


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Download ChEMBL SQLite and generate target metadata."
    )
    parser.add_argument(
        "--output", "-o",
        default="compounds-targets-data",
        help="Output directory (default: compounds-targets-data)",
    )
    parser.add_argument(
        "--yes", "-y",
        action="store_true",
        help="Skip confirmation prompt",
    )
    parser.add_argument(
        "--url",
        help=f"Custom download URL (default: {CHEMBL_BASE_URL}{CHEMBL_DB_NAME})",
    )
    parser.add_argument(
        "--skip-target-info",
        action="store_true",
        help="Skip target_info.json generation",
    )
    args = parser.parse_args()

    output_dir = Path(args.output).resolve()
    os.makedirs(output_dir, exist_ok=True)

    download_url = args.url or f"{CHEMBL_BASE_URL}{CHEMBL_DB_NAME}"
    archive_name = os.path.basename(urlparse(download_url).path)
    archive_path = output_dir / archive_name

    # ── Step 1: Check if .db already exists ────────────────────────────────
    db_path = find_latest_db(output_dir)
    if db_path:
        size = os.path.getsize(db_path)
        print(f"✅ ChEMBL database already exists: {db_path.name} ({format_size(size)})")

        # Check if target_info.json needs (re)generation
        info_path = output_dir / "target_info.json"
        if not args.skip_target_info:
            if info_path.exists():
                print(f"✅ target_info.json already exists ({format_size(os.path.getsize(info_path))})")
                print("\nAll ready. To re-download, delete the existing files first.")
            else:
                print("\ntarget_info.json not found, generating ...")
                generate_target_info(db_path, output_dir)
                print("\n✅ All ready.")
        else:
            print("\nSkipping target_info.json (--skip-target-info)")
        return

    # ── Step 2: Check if archive already downloaded ────────────────────────
    if archive_path.exists():
        size = os.path.getsize(archive_path)
        print(f"📦 Archive already downloaded: {archive_path.name} ({format_size(size)})")
        print("   Skipping download, will extract directly.")
    else:
        # ── Step 3: Confirm ────────────────────────────────────────────────
        if not args.yes:
            print("=" * 60)
            print("  ChEMBL Database Download")
            print("=" * 60)
            print()
            print(f"  URL:      {download_url}")
            print(f"  Size:     ~{EXPECTED_SIZE_GB} GB (compressed)")
            print(f"  After:    ~29 GB (extracted)")
            print(f"  Save to:  {output_dir}/")
            print()
            print("  The database is required for target prediction.")
            print("  target_info.json will be auto-generated after extraction.")
            print()

            answer = input("  Download now? [Y/n] ").strip().lower()
            if answer and answer not in ("y", "yes"):
                print("\n  Cancelled. Run this script again when ready.")
                sys.exit(0)

        # ── Step 4: Download ───────────────────────────────────────────────
        print(f"\n⬇  Downloading {archive_name} ...")
        print(f"   From: {download_url}")
        print()

        try:
            urllib.request.urlretrieve(
                download_url,
                archive_path,
                reporthook=progress_hook,
            )
            print()
            final_size = os.path.getsize(archive_path)
            print(f"\n✅ Download complete: {format_size(final_size)}")
        except Exception as e:
            print(f"\n❌ Download failed: {e}")
            if archive_path.exists():
                archive_path.unlink()
            sys.exit(1)

    # ── Step 5: Extract ────────────────────────────────────────────────────
    try:
        db_path = extract_tar_gz(archive_path, output_dir)
    except Exception as e:
        print(f"\n❌ Extraction failed: {e}")
        sys.exit(1)

    # ── Step 6: Generate target_info.json ──────────────────────────────────
    if not args.skip_target_info:
        generate_target_info(db_path, output_dir)

    # ── Step 7: Verify ─────────────────────────────────────────────────────
    print(f"\n✅ ChEMBL database ready: {db_path.name} ({format_size(os.path.getsize(db_path))})")
    print(f"Next step: generate fingerprint database → see SKILL.md")


if __name__ == "__main__":
    main()
