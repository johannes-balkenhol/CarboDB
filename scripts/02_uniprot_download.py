#!/usr/bin/env python3
"""
02_uniprot_download.py
======================
CarboxyDB — Step 2: Large-scale UniProt download.

Downloads SwissProt (reviewed) + TrEMBL (unreviewed) sequences for all
CO2-related EC classes discovered in Script 01, plus a matching negative pool.

Confirmed working API:
  URL: https://rest.uniprot.org/uniprotkb/search
  Query: (ec:{ec})                  ← confirmed working Dec 2025
  Format: tsv
  Fields: accession,organism_name,protein_name,sequence,length,reviewed,
          lineage_ids,annotation_score,go_ids
  Pagination: cursor from Link: <...>; rel="next" response header

Scale targets:
  SwissProt positives: ~50,000  (all reviewed entries for CO2 EC classes)
  TrEMBL positives:   ~600,000  (capped per EC, sorted by annotation score)
  Negatives:          ~650,000  (matching total positive count)

Usage:
    # Reads EC list from data/raw/brenda/co2_ec_classes_*.txt (Script 01 output)
    python scripts/02_uniprot_download.py

    # Or specify EC file explicitly:
    python scripts/02_uniprot_download.py --ec-file data/raw/brenda/co2_ec_classes_20250101_120000.txt
"""

import argparse
import re
import sys
import time
from io import StringIO
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

TS = datetime.now().strftime("%Y%m%d_%H%M%S")

SP_OUT  = Path("data/raw/uniprot/swissport")
TB_OUT  = Path("data/raw/uniprot/trembl")
NEG_OUT = Path("data/raw/uniprot/negatives")
for d in [SP_OUT, TB_OUT, NEG_OUT]:
    d.mkdir(parents=True, exist_ok=True)

# UniProt REST API — confirmed working endpoint and field set
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,organism_name,protein_name,sequence,length,reviewed,lineage_ids,annotation_score"
BATCH  = 500   # max per page

# Confirmed-negative EC classes (do NOT interact with CO2)
# Same list as Script 01, kept here for self-containedness
NEGATIVE_EC = [
    "1.1.1.1","1.1.1.27","1.1.1.37","1.1.1.42","1.1.1.44","1.1.1.49",
    "1.2.1.12","1.4.1.2","1.4.1.3","1.6.5.3","1.8.1.4","1.9.3.1",
    "1.11.1.6","1.11.1.7","1.14.13.39","1.15.1.1",
    "2.1.1.37","2.2.1.1","2.3.1.9","2.4.1.1","2.4.2.1",
    "2.5.1.18","2.6.1.1","2.6.1.2","2.7.1.1","2.7.1.11","2.7.1.40",
    "2.7.1.69","2.7.4.3","2.7.7.7","2.7.7.48",
    "3.1.1.3","3.1.3.1","3.1.3.2","3.2.1.1","3.2.1.17","3.2.1.21",
    "3.2.1.23","3.4.11.1","3.4.21.4","3.4.21.1","3.4.23.1",
    "3.5.1.1","3.5.1.2","3.5.4.4","3.6.1.3","3.6.4.12",
    "4.1.2.13","4.2.1.2","4.2.1.3","4.2.1.11","4.3.1.3",
    "5.1.3.1","5.3.1.1","5.3.1.9","5.4.2.1","5.4.2.2","5.4.99.2",
    "6.1.1.1","6.1.1.2","6.1.1.5","6.2.1.1",
    "6.3.2.1","6.3.2.2","6.5.1.1","6.5.1.2",
]


# ─────────────────────────────────────────────────────────────────────────────
# Core download function — confirmed working (cursor pagination from Link header)
# ─────────────────────────────────────────────────────────────────────────────

def fetch_uniprot(query: str, max_results: int = 1_000_000,
                  reviewed_only: bool = False) -> pd.DataFrame:
    """
    Fetch sequences from UniProt REST API using cursor pagination.

    Confirmed working API pattern (from March 2026 MEME notebook):
        url = 'https://rest.uniprot.org/uniprotkb/search'
        params = {'query': query, 'format': 'tsv', 'fields': fields, 'size': 500}
        cursor from response Link header: cursor=([^&>]+)

    Returns DataFrame with columns:
        accession, organism, protein_name, sequence, length,
        reviewed, lineage_ids, annotation_score
    """
    if reviewed_only:
        query = f"({query}) AND reviewed:true"

    frames = []
    fetched = 0
    cursor = None
    retries = 0
    MAX_RETRIES = 5

    while fetched < max_results:
        batch = min(BATCH, max_results - fetched)
        params = {
            "query":  query,
            "format": "tsv",
            "fields": FIELDS,
            "size":   batch,
        }
        if cursor:
            params["cursor"] = cursor

        try:
            r = requests.get(UNIPROT_SEARCH, params=params, timeout=60)
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            retries += 1
            if retries > MAX_RETRIES:
                print(f"  ERROR after {MAX_RETRIES} retries: {e}")
                break
            wait = 2 ** retries
            print(f"  Retry {retries}/{MAX_RETRIES} in {wait}s: {e}")
            time.sleep(wait)
            continue

        retries = 0  # reset on success

        try:
            chunk = pd.read_csv(StringIO(r.text), sep="\t")
        except Exception:
            break

        if chunk.empty:
            break

        # Normalise column names (UniProt sometimes changes capitalisation)
        chunk.columns = [c.lower().replace(" ", "_") for c in chunk.columns]
        # Rename to our standard names
        rename = {
            "entry": "accession",
            "entry_name": "entry_name",
            "organism": "organism",
            "organism_(id)": "taxonomy_id",
            "reviewed": "reviewed",
            "annotation_score": "annotation_score",
            "sequence": "sequence",
            "length": "length",
            "protein_names": "protein_name",
            "gene_names": "gene_name",
            "lineage_ids": "lineage_ids",
        }
        for old, new in rename.items():
            if old in chunk.columns and old != new:
                chunk = chunk.rename(columns={old: new})

        if "accession" not in chunk.columns and len(chunk.columns) > 0:
            chunk.columns = ["accession"] + list(chunk.columns[1:])

        chunk = chunk.dropna(subset=["sequence"] if "sequence" in chunk.columns else [])
        frames.append(chunk)
        fetched += len(chunk)

        # Cursor pagination
        link = r.headers.get("Link", "")
        if 'rel="next"' not in link:
            break
        m = re.search(r"cursor=([^&>]+)", link)
        cursor = m.group(1) if m else None
        if not cursor:
            break

        time.sleep(0.3)

    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    if "accession" in df.columns:
        df = df.drop_duplicates(subset="accession")
    return df.reset_index(drop=True)


def parse_uniprot_fasta(text: str) -> dict:
    """Parse FASTA format UniProt response → {uid: sequence}."""
    seqs = {}
    uid = seq_lines = None
    for line in text.splitlines():
        if line.startswith(">"):
            if uid and seq_lines:
                seqs[uid] = "".join(seq_lines)
            parts = line.split("|")
            uid = parts[1] if len(parts) >= 2 else line[1:].split()[0]
            seq_lines = []
        elif line.strip():
            if seq_lines is not None:
                seq_lines.append(line.strip())
    if uid and seq_lines:
        seqs[uid] = "".join(seq_lines)
    return seqs


# ─────────────────────────────────────────────────────────────────────────────
# Download positives
# ─────────────────────────────────────────────────────────────────────────────

def download_positives(ec_list: list[str]):
    """Download SwissProt and TrEMBL for all CO2 EC classes."""

    sp_rows, tb_rows = [], []

    print(f"\n{'='*70}")
    print(f"DOWNLOADING POSITIVES ({len(ec_list)} CO2 EC classes)")
    print(f"{'='*70}")

    for ec in tqdm(ec_list, desc="EC classes"):

        # ── SwissProt (reviewed, no size cap) ──────────────────────────────
        query_sp = f"(ec:{ec}) AND reviewed:true"
        df_sp = fetch_uniprot(query_sp, max_results=100_000)
        if not df_sp.empty:
            df_sp["ec_number"] = ec
            df_sp["source"]    = "swissport"
            df_sp["label"]     = 1
            sp_rows.append(df_sp)
        time.sleep(0.5)

        # ── TrEMBL (unreviewed, capped at 100k per EC, sorted by score) ────
        query_tb = f"(ec:{ec}) AND reviewed:false"
        df_tb = fetch_uniprot(query_tb, max_results=100_000)
        if not df_tb.empty:
            df_tb["ec_number"] = ec
            df_tb["source"]    = "trembl"
            df_tb["label"]     = 1
            tb_rows.append(df_tb)
        time.sleep(0.5)

    df_sp_all = pd.concat(sp_rows, ignore_index=True) if sp_rows else pd.DataFrame()
    df_tb_all = pd.concat(tb_rows, ignore_index=True) if tb_rows else pd.DataFrame()

    # Save
    if not df_sp_all.empty:
        df_sp_all.to_csv(SP_OUT / f"swissport_positives_{TS}.tsv", sep="\t", index=False)
        _write_fasta(df_sp_all, SP_OUT / f"swissport_positives_{TS}.fasta", label=1)

    if not df_tb_all.empty:
        df_tb_all.to_csv(TB_OUT / f"trembl_positives_{TS}.tsv", sep="\t", index=False)
        _write_fasta(df_tb_all, TB_OUT / f"trembl_positives_{TS}.fasta", label=1)

    print(f"\n  SwissProt positives: {len(df_sp_all):,}")
    print(f"  TrEMBL positives:    {len(df_tb_all):,}")
    return df_sp_all, df_tb_all


# ─────────────────────────────────────────────────────────────────────────────
# Download negatives
# ─────────────────────────────────────────────────────────────────────────────

def download_negatives(n_target: int):
    """
    Download negative sequences from confirmed non-CO2 EC classes.
    Target: n_target sequences (= total positives from all sources).
    Mix reviewed + unreviewed, stratified by EC.
    """
    print(f"\n{'='*70}")
    print(f"DOWNLOADING NEGATIVES (target: {n_target:,})")
    print(f"{'='*70}")

    per_ec = max(500, n_target // len(NEGATIVE_EC) + 200)
    neg_rows = []
    seen_uids = set()

    for ec in tqdm(NEGATIVE_EC, desc="Negative ECs"):
        query = f"(ec:{ec})"
        df = fetch_uniprot(query, max_results=per_ec)
        if not df.empty:
            df["ec_number"] = ec
            df["source"]    = "negative"
            df["label"]     = 0
            # Remove any UID already in negatives
            if "accession" in df.columns:
                df = df[~df["accession"].isin(seen_uids)]
                seen_uids.update(df["accession"].tolist())
            neg_rows.append(df)
        time.sleep(0.3)
        if len(seen_uids) >= n_target:
            break

    df_neg = pd.concat(neg_rows, ignore_index=True) if neg_rows else pd.DataFrame()
    if not df_neg.empty:
        df_neg.to_csv(NEG_OUT / f"negatives_{TS}.tsv", sep="\t", index=False)
        _write_fasta(df_neg, NEG_OUT / f"negatives_{TS}.fasta", label=0)

    print(f"\n  Negative sequences: {len(df_neg):,}")
    return df_neg


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _write_fasta(df: pd.DataFrame, path: Path, label: int):
    uid_col = "accession" if "accession" in df.columns else df.columns[0]
    seq_col = "sequence"
    ec_col  = "ec_number" if "ec_number" in df.columns else None
    org_col = "organism"  if "organism"  in df.columns else None

    with open(path, "w") as f:
        for _, row in df.iterrows():
            seq = str(row.get(seq_col, "") or "").strip()
            if not seq:
                continue
            uid  = str(row.get(uid_col, ""))
            ec   = str(row.get(ec_col,  "")) if ec_col  else ""
            org  = str(row.get(org_col, "")).replace(" ", "_")[:40] if org_col else ""
            f.write(f">{uid}|{ec}|{org}|label={label}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def load_ec_list(ec_file: Path) -> list[str]:
    """Load EC list from Script 01 output."""
    with open(ec_file) as f:
        return [line.strip() for line in f if line.strip()]


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ec-file", default=None,
                        help="Path to EC list from Script 01. "
                             "Auto-detects latest if not given.")
    args = parser.parse_args()

    # Find EC file
    if args.ec_file:
        ec_file = Path(args.ec_file)
    else:
        brenda_dir = Path("data/raw/brenda")
        ec_files = sorted(brenda_dir.glob("co2_ec_classes_*.txt"))
        if not ec_files:
            print("ERROR: No co2_ec_classes_*.txt found. Run 01_brenda_download.py first.")
            sys.exit(1)
        ec_file = ec_files[-1]

    print(f"Using EC list: {ec_file}")
    ec_list = load_ec_list(ec_file)
    print(f"CO2 EC classes: {len(ec_list)}")

    # Download positives
    df_sp, df_tb = download_positives(ec_list)

    # Download negatives targeting ~equal size to positives
    n_pos = len(df_sp) + len(df_tb)
    df_neg = download_negatives(n_target=n_pos)

    print("\n" + "=" * 70)
    print("UNIPROT DOWNLOAD COMPLETE")
    print("=" * 70)
    print(f"  SwissProt positives: {len(df_sp):,}")
    print(f"  TrEMBL positives:    {len(df_tb):,}")
    print(f"  Negatives:           {len(df_neg):,}")
    print(f"\n  Next step: python scripts/03_merge_all_sources.py")


if __name__ == "__main__":
    main()
