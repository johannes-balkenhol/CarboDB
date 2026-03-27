#!/usr/bin/env python3
"""
00_setup_project.py
===================
CarboxyDB — Run this ONCE before anything else.

Creates the full project folder structure and checks that required
external tools are installed (HMMER3, BLAST+, CD-HIT, Python packages).

Usage:
    python 00_setup_project.py
    python 00_setup_project.py --skip-tool-check
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# ── Folder tree ────────────────────────────────────────────────────────────────
# Every directory the pipeline needs, created up front so no script ever has to
# mkdir on the fly (except timestamped run dirs).
DIRS = [
    # Raw downloads — never modified after creation
    "data/raw/brenda",
    "data/raw/uniprot/swissprot",
    "data/raw/uniprot/trembl",
    "data/raw/uniprot/negatives",

    # Intermediate processing
    "data/interim",

    # Primary merged dataset (single source of truth)
    "data/primary",

    # Per-layer feature files (TSV, one row per CDB_ID)
    "data/features/composition",
    "data/features/domains",
    "data/features/motifs",
    "data/features/blast",
    "data/features/esm2",
    "data/features/meme",       # PENDING: filled by MEME subproject

    # ML-ready splits
    "data/ml",

    # Benchmark outputs
    "data/benchmark",

    # SHAP outputs (Figure 3)
    "data/shap",

    # HMMER / BLAST reference databases
    "data/dbs/pfam",
    "data/dbs/blast",
    "data/dbs/prosite",

    # Trained model artefacts
    "models",

    # SQLite database
    "database",

    # Scripts (already exists if you cloned the repo)
    "scripts",

    # Logs — one file per script run
    "logs",
]

# ── External tool requirements ─────────────────────────────────────────────────
TOOLS = {
    "hmmscan":      "HMMER3  — sudo apt install hmmer  OR  conda install -c bioconda hmmer",
    "blastp":       "BLAST+  — sudo apt install ncbi-blast+  OR  conda install -c bioconda blast",
    "makeblastdb":  "BLAST+  — (same package as blastp)",
    "cd-hit":       "CD-HIT  — sudo apt install cd-hit  OR  conda install -c bioconda cd-hit",
    "fimo":         "MEME suite (optional, Layer D) — conda install -c bioconda meme",
}

# ── Python package requirements ────────────────────────────────────────────────
PY_PACKAGES = [
    "requests", "tqdm", "pandas", "numpy", "Bio",
    "xgboost", "sklearn", "zeep",   # zeep for BRENDA SOAP API
]
PY_OPTIONAL = ["shap", "torch", "esm"]   # ESM-2 embeddings


# ── Helpers ────────────────────────────────────────────────────────────────────

def _ok(msg):  print(f"  \033[32m✓\033[0m  {msg}")
def _warn(msg): print(f"  \033[33m⚠\033[0m  {msg}")
def _err(msg):  print(f"  \033[31m✗\033[0m  {msg}")


def create_dirs():
    print("\n── Creating folder structure ──")
    root = Path(__file__).resolve().parent
    for d in DIRS:
        p = root / d
        p.mkdir(parents=True, exist_ok=True)
    _ok(f"All {len(DIRS)} directories ready under {root}/")
    return root


def check_tools():
    print("\n── Checking external tools ──")
    missing = []
    for tool, hint in TOOLS.items():
        optional = "optional" in hint.lower()
        if shutil.which(tool):
            _ok(f"{tool}")
        elif optional:
            _warn(f"{tool} not found (optional)  →  {hint}")
        else:
            _err(f"{tool} not found  →  {hint}")
            missing.append(tool)
    if missing:
        print(f"\n  REQUIRED tools missing: {missing}")
        print("  Install them before running scripts 04 (domains) and 06 (CD-HIT).")
    else:
        print("  All required tools found.")
    return missing


def check_python():
    print("\n── Checking Python packages ──")
    missing = []
    for pkg in PY_PACKAGES:
        import importlib
        try:
            importlib.import_module(pkg.replace("-", "_"))
            _ok(pkg)
        except ImportError:
            _err(f"{pkg} missing  →  pip install {pkg}")
            missing.append(pkg)
    for pkg in PY_OPTIONAL:
        import importlib
        try:
            importlib.import_module(pkg)
            _ok(f"{pkg} (optional)")
        except ImportError:
            _warn(f"{pkg} not installed (optional — needed for ESM-2 embeddings)")
    if missing:
        print(f"\n  Install missing packages:  pip install {' '.join(missing)}")
    return missing


def check_config(root: Path):
    print("\n── Checking config.py ──")
    cfg = root / "config.py"
    if cfg.exists():
        _ok("config.py found")
    else:
        _err("config.py not found — copy config.py to the project root")


def write_gitignore(root: Path):
    gi = root / ".gitignore"
    if gi.exists():
        return
    gi.write_text(
        "# CarboxyDB — auto-generated\n"
        "data/raw/\n"
        "data/primary/\n"
        "data/features/\n"
        "data/ml/\n"
        "data/dbs/\n"
        "database/*.sqlite\n"
        "models/*.pkl\n"
        "models/*.json\n"
        "logs/\n"
        "__pycache__/\n"
        "*.pyc\n"
        ".env\n"
        "credentials.json\n"
    )
    _ok(".gitignore written (raw data and models excluded from git)")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-tool-check", action="store_true")
    args = ap.parse_args()

    print("=" * 70)
    print("CarboxyDB — Project Setup")
    print("=" * 70)

    root = create_dirs()
    write_gitignore(root)
    check_config(root)

    tool_issues = []
    if not args.skip_tool_check:
        tool_issues = check_tools()

    py_issues = check_python()

    print("\n" + "=" * 70)
    if not tool_issues and not py_issues:
        print("✓ Setup complete — ready to run the pipeline.")
        print("\n  Next step: python scripts/01_brenda_download.py")
    else:
        print("⚠ Setup complete with warnings — fix issues above before running.")
    print("=" * 70)


if __name__ == "__main__":
    main()
