#!/usr/bin/env python3
"""Join FULL marine v5 predictions with sample metadata and taxonomy."""
import gzip, re, sys
from pathlib import Path
import pandas as pd

ROOT     = Path.home() / "Projects_shared/CarboDB_v3"
PRED_IN  = ROOT / "data/gmgc/results/predictions_full.parquet"
META_F   = ROOT / "data/gmgc/downloads/GMGC10.meta.tsv.gz"
SAMPLE_F = ROOT / "data/gmgc/downloads/GMGC10.sample.meta.tsv.gz"
TAXON_F  = ROOT / "data/gmgc/downloads/GMGC10.taxonomy.tsv.gz"
OUT_FILE = ROOT / "data/gmgc/results/predictions_joined_full.parquet"

print(f"Loading {PRED_IN}", file=sys.stderr)
df = pd.read_parquet(PRED_IN)
print(f"  {len(df):,} predictions", file=sys.stderr)

df["gmgc_id"] = df["seq_id"].str.replace(r"_\d+$", "", regex=True)
target_ids = set(df["gmgc_id"].unique())
print(f"Looking up {len(target_ids):,} unique unigenes", file=sys.stderr)

meta_rows = []
with gzip.open(META_F, "rt") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        if parts[0] in target_ids:
            meta_rows.append((parts[0], parts[1], int(parts[2]), parts[3]=="True"))
meta_df = pd.DataFrame(meta_rows, columns=["gmgc_id","original_name","size_nt","is_complete"])
print(f"  Matched in meta: {len(meta_df):,}", file=sys.stderr)

def extract_sample(orig):
    if not isinstance(orig, str): return None
    for sm in ["_scaffold", "_NODE", "_contig", "_k119_", "_k127_"]:
        idx = orig.find(sm)
        if idx > 0: return orig[:idx]
    m = re.match(r"^(.+?)_\d+$", orig)
    return m.group(1) if m else orig

meta_df["sample_guess"] = meta_df["original_name"].apply(extract_sample)

print("Loading sample metadata", file=sys.stderr)
sample_df = pd.read_csv(SAMPLE_F, sep="\t")
meta_df = meta_df.merge(sample_df, left_on="sample_guess", right_on="original_name",
                         suffixes=("","_sample"), how="left")
n_matched = meta_df["sample_id"].notna().sum()
print(f"  Sample-matched: {n_matched:,} / {len(meta_df):,} "
      f"({100*n_matched/len(meta_df):.1f}%)", file=sys.stderr)

print("Loading taxonomy (large file)", file=sys.stderr)
taxa_rows = []
with gzip.open(TAXON_F, "rt") as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        if parts[0] in target_ids:
            taxa_rows.append((parts[0], parts[1], parts[2], parts[3]))
taxa_df = pd.DataFrame(taxa_rows, columns=["gmgc_id","tax_name","tax_rank","tax_id"])
print(f"  Taxonomy matched: {len(taxa_df):,}", file=sys.stderr)

out = df.merge(meta_df.drop(columns=["original_name_sample"], errors="ignore"),
                on="gmgc_id", how="left")
out = out.merge(taxa_df, on="gmgc_id", how="left")

print(f"\nFinal:  {len(out):,} predictions", file=sys.stderr)
print(f"  GPS:       {out['latitude'].notna().sum():,}", file=sys.stderr)
print(f"  Habitat:   {out['habitat'].notna().sum():,}", file=sys.stderr)
print(f"  Taxonomy:  {out['tax_name'].notna().sum():,}", file=sys.stderr)
out.to_parquet(OUT_FILE, index=False)
print(f"\nWrote {OUT_FILE}", file=sys.stderr)
