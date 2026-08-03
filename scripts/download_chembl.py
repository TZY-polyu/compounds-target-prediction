#!/usr/bin/env python3
"""
ChEMBL Database Downloader
==========================
Automatically download and extract the latest ChEMBL SQLite database.

Usage:
    python scripts/download_chembl.py
    python scripts/download_chembl.py --output /path/to/dir
    python scripts/download_chembl.py --yes   # skip confirmation

The script will:
    1. Check if chembl_*.db already exists in the output directory
    2. Ask for confirmation before downloading (~5.4GB compressed, ~29GB extracted)
    3. Download from the official ChEMBL FTP mirror
    4. Extract the .tar.gz archive
"""

import argparse
import os
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


def find_existing_db(output_dir: Path) -> list[str]:
    """Find existing ChEMBL SQLite files in the output directory."""
    existing = sorted(output_dir.glob("chembl_*.db"))
    return [str(p) for p in existing]


def find_existing_archive(output_dir: Path) -> list[str]:
    """Find existing downloaded archives."""
    existing = sorted(output_dir.glob("chembl_*_sqlite.tar.gz"))
    return [str(p) for p in existing]


def extract_tar_gz(archive_path: Path, output_dir: Path) -> Path:
    """Extract .tar.gz and return the .db file path."""
    print(f"\nExtracting {archive_path.name} ...")
    with tarfile.open(archive_path, "r:gz") as tar:
        members = tar.getmembers()
        # Find the .db file
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
            # No .db found, extract all
            print("  Extracting all files ...")
            tar.extractall(path=output_dir, filter="data")
            db_files = sorted(output_dir.glob("chembl_*.db"))
            if db_files:
                print(f"  Extracted → {db_files[0]}")
                return db_files[0]
            else:
                raise RuntimeError("No .db file found in archive")


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Download the latest ChEMBL SQLite database for target prediction."
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
    args = parser.parse_args()

    output_dir = Path(args.output).resolve()
    os.makedirs(output_dir, exist_ok=True)

    download_url = args.url or f"{CHEMBL_BASE_URL}{CHEMBL_DB_NAME}"
    archive_name = os.path.basename(urlparse(download_url).path)
    archive_path = output_dir / archive_name

    # ── Step 1: Check if .db already exists ────────────────────────────────
    existing_dbs = find_existing_db(output_dir)
    if existing_dbs:
        print("✅ ChEMBL database already exists:")
        for db in existing_dbs:
            size = os.path.getsize(db)
            print(f"   {db} ({format_size(size)})")
        print("\nSkipping download. To re-download, delete the existing .db file first.")
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
            print("  This download is required for target prediction.")
            print("  Without it, the SEA+TC engine cannot run.")
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
            print()  # newline after progress bar
            final_size = os.path.getsize(archive_path)
            print(f"\n✅ Download complete: {format_size(final_size)}")
        except Exception as e:
            print(f"\n❌ Download failed: {e}")
            # Clean up partial download
            if archive_path.exists():
                archive_path.unlink()
            sys.exit(1)

    # ── Step 5: Extract ────────────────────────────────────────────────────
    try:
        db_path = extract_tar_gz(archive_path, output_dir)
    except Exception as e:
        print(f"\n❌ Extraction failed: {e}")
        sys.exit(1)

    # ── Step 6: Verify ─────────────────────────────────────────────────────
    if db_path.exists():
        db_size = os.path.getsize(db_path)
        print(f"\n✅ ChEMBL database ready: {db_path.name} ({format_size(db_size)})")
        print(f"\nNext step: generate fingerprint database → see SKILL.md")
    else:
        print(f"\n❌ Error: .db file not found after extraction")
        sys.exit(1)


if __name__ == "__main__":
    main()
