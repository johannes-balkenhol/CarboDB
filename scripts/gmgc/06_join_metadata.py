#!/usr/bin/env python3
"""
Join the v5 predictions with GMGC sample metadata and taxonomy.

Steps:
  1. Strip _1 suffix from prediction seq_ids (transeq added it)
  2. Look up each unigene's source sample via GMGC10.meta.tsv.gz
     (gives us "original_name" which maps to sample_id in sample.meta)
     NOTE: a unigene's representative is from ONE sample (the assembly source)
  3. Join sample metadata (GPS, country, etc.)
  4. Join taxonomy assignments

Reads:  data/gmgc/results/predictions_pilot.parquet
        data/gmgc/downloads/GMGC10.meta.tsv.gz
        data/gmgc/downloads/GMGC10.sample.meta.tsv.gz
        data/gmgc/downloads/GMGC10.taxonomy.tsv.gz
Writes: data/gmgc/results/predictions_joined.parquet
"""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

ROOT     = Path.home() / "Projects_shared/CarboDB_v3"
PRED_IN  = ROOT / "data/gmgc/results/predictions_pilot.parquet"
META_F   = ROOT / "data/gmgc/downloads/GMGC10.meta.tsv.gz"
SAMPLE_F = ROOT / "data/gmgc/downloads/GMGC10.sample.meta.tsv.gz"
TAXON_F  = ROOT / "data/gmgc/downloads/GMGC10.taxonomy.tsv.gz"
OUT_FILE = ROOT / "data/gmgc/results/predictions_joined.parquet"

# ── 1. LOAD PREDICTIONS, STRIP _N FROM SEQ_ID ────────────────────────
print(f"Loading predictions from {PRED_IN}", file=sys.stderr)
df = pd.read_parquet(PRED_IN)
print(f"  {len(df):,} predictions", file=sys.stderr)

# Strip transeq's frame suffix: GMGC10.xxx.UNKNOWN_1 -> GMGC10.xxx.UNKNOWN
df["gmgc_id"] = df["seq_id"].str.replace(r"_\d+$", "", regex=True)
print(f"  Example seq_id -> gmgc_id: {df['seq_id'].iloc[0]} -> {df['gmgc_id'].iloc[0]}",
      file=sys.stderr)

# ── 2. SAMPLE LOOKUP (gmgc_id -> original_name -> sample_id) ─────────
# GMGC10.meta.tsv.gz columns: new_names, original_name, sizes, is_complete
# We need original_name → strip ORF suffix to get sample identifier
# e.g., original_name "TARA_ASW_039_C_4_60_5_scaffold10_2" → sample "TARA_ASW_039"

# Pre-filter to just the IDs we need
target_ids = set(df["gmgc_id"].unique())
print(f"  Looking up source samples for {len(target_ids):,} unique unigenes", file=sys.stderr)

meta_rows = []
with gzip.open(META_F, "rt") as f:
    next(f)  # header
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        new_name, orig_name, size, is_complete = parts[:4]
        if new_name in target_ids:
            meta_rows.append((new_name, orig_name, int(size), is_complete == "True"))

meta_df = pd.DataFrame(meta_rows, columns=["gmgc_id", "original_name", "size_nt", "is_complete"])
print(f"  {len(meta_df):,} unigenes matched in meta.tsv.gz", file=sys.stderr)

# Extract sample_id from original_name (the source assembly's sample name)
# Several patterns exist; the most common: TARA_ASW_039_..., ExDog10_..., etc.
# Sample names in sample.meta.tsv.gz are like "ExDog1" or "TARA_039_DCM"
# A reasonable heuristic: take the first N tokens before _scaffold or _NODE.
def extract_sample(orig):
    if not isinstance(orig, str):
        return None
    # Strip everything from _scaffold onwards
    for suffix_marker in ["_scaffold", "_NODE", "_contig", "_k119_", "_k127_"]:
        idx = orig.find(suffix_marker)
        if idx > 0:
            return orig[:idx]
    # Fallback: take everything up to last _<number>
    m = re.match(r"^(.+?)_\d+$", orig)
    return m.group(1) if m else orig

meta_df["sample_guess"] = meta_df["original_name"].apply(extract_sample)
print(f"  Example: {meta_df['original_name'].iloc[0]} -> {meta_df['sample_guess'].iloc[0]}",
      file=sys.stderr)

# ── 3. JOIN SAMPLE METADATA ──────────────────────────────────────────
print(f"Loading sample metadata", file=sys.stderr)
sample_df = pd.read_csv(SAMPLE_F, sep="\t")
print(f"  {len(sample_df):,} samples", file=sys.stderr)

# Try matching sample_guess to either sample_id or original_name (the public name)
meta_df = meta_df.merge(
    sample_df, left_on="sample_guess", right_on="original_name",
    suffixes=("", "_sample"), how="left",
)

n_matched = meta_df["sample_id"].notna().sum()
print(f"  {n_matched:,} / {len(meta_df):,} unigenes matched to a sample "
      f"({100*n_matched/len(meta_df):.1f}%)", file=sys.stderr)

# ── 4. JOIN TAXONOMY (one assignment per gmgc_id, can be 'no rank') ──
print(f"Loading taxonomy (this is large, takes a couple of minutes)", file=sys.stderr)
taxa_rows = []
with gzip.open(TAXON_F, "rt") as f:
    next(f)  # header
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        gene, name, rank, taxid = parts[:4]
        if gene in target_ids:
            taxa_rows.append((gene, name, rank, taxid))
taxa_df = pd.DataFrame(taxa_rows, columns=["gmgc_id", "tax_name", "tax_rank", "tax_id"])
print(f"  {len(taxa_df):,} taxonomy assignments matched", file=sys.stderr)

# ── 5. MERGE EVERYTHING ──────────────────────────────────────────────
out = df.merge(meta_df.drop(columns=["original_name_sample"], errors="ignore"),
               on="gmgc_id", how="left")
out = out.merge(taxa_df, on="gmgc_id", how="left")

print(f"\nFinal joined predictions: {len(out):,}", file=sys.stderr)
print(f"  With GPS:       {out['latitude'].notna().sum():,}", file=sys.stderr)
print(f"  With habitat:   {out['habitat'].notna().sum():,}", file=sys.stderr)
print(f"  With taxonomy:  {out['tax_name'].notna().sum():,}", file=sys.stderr)

out.to_parquet(OUT_FILE, index=False)
print(f"\nWritten {OUT_FILE}", file=sys.stderr)
