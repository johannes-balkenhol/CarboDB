#!/usr/bin/env python3
"""
02_prepare_training_data.py
===========================
CarboxyDB — Step 2: Build balanced training / validation / test datasets
from the raw BRENDA output produced by 01_brenda_download.py.

What this does:
  A. Loads positives TSV + negatives TSV from data/raw/brenda/
  B. Validates sequences (AA alphabet, length 50–5000)
  C. Adds Km values to positive sequences that have them
  D. Balances: positives with Km form the gold-standard core;
     all positives are kept; negatives are sampled to ~5× positives
  E. Deduplicates at 100% identity
  F. Writes three outputs:
       data/interim/training_data.tsv  — balanced training set with all metadata
       data/interim/training_data.fasta
       data/interim/km_training.tsv    — only rows with experimental Km (for regression)

After this you run:
  03_annotate_features.py   → extract all sequence features
  04_train_models.py        → train XGBoost v3 and v5

Usage:
    python scripts/02_prepare_training_data.py
    # (auto-detects latest BRENDA files in data/raw/brenda/)
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from datetime import datetime

TS = datetime.now().strftime("%Y%m%d_%H%M%S")

RAW_DIR  = Path("data/raw/brenda")
OUT_DIR  = Path("data/interim")
OUT_DIR.mkdir(parents=True, exist_ok=True)

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 50, 5000


# ── Helpers ──────────────────────────────────────────────────────────────────

def latest(pattern: str) -> Path:
    """Find the most recent file matching a glob pattern."""
    files = sorted(RAW_DIR.glob(pattern))
    if not files:
        print(f"ERROR: No file matching {RAW_DIR / pattern}")
        sys.exit(1)
    return files[-1]

def is_valid(seq: str) -> bool:
    seq = str(seq).upper().strip()
    if len(seq) < MIN_LEN or len(seq) > MAX_LEN:
        return False
    chars = set(seq)
    return len(chars - VALID_AA) == 0 and len(chars) >= 5

def clean_seq(seq: str) -> str:
    return str(seq).upper().strip()


# ── Load raw files ────────────────────────────────────────────────────────────

def load_raw():
    pos_file = latest("brenda_positives_*.tsv")
    neg_file = latest("brenda_negatives_*.tsv")
    km_file  = latest("brenda_km_curated_*.tsv")

    print(f"Loading positives:  {pos_file.name}")
    print(f"Loading negatives:  {neg_file.name}")
    print(f"Loading curated Km: {km_file.name}")

    df_pos = pd.read_csv(pos_file, sep="\t", dtype=str).fillna("")
    df_neg = pd.read_csv(neg_file, sep="\t", dtype=str).fillna("")
    df_km  = pd.read_csv(km_file,  sep="\t", dtype=str).fillna("")

    print(f"\nRaw counts:")
    print(f"  Positives: {len(df_pos):,}")
    print(f"  Negatives: {len(df_neg):,}")
    print(f"  Km (curated): {len(df_km):,}")

    return df_pos, df_neg, df_km


# ── Validate ─────────────────────────────────────────────────────────────────

def validate(df: pd.DataFrame, name: str) -> pd.DataFrame:
    before = len(df)
    df = df.copy()
    df["sequence"] = df["sequence"].apply(clean_seq)
    mask = df["sequence"].apply(is_valid)
    df = df[mask].reset_index(drop=True)
    print(f"  {name}: {before:,} → {len(df):,} after validation "
          f"({before - len(df):,} removed)")
    return df


# ── Deduplicate ───────────────────────────────────────────────────────────────

def dedup(df: pd.DataFrame, name: str) -> pd.DataFrame:
    before = len(df)
    df = df.drop_duplicates(subset="uniprot_id").reset_index(drop=True)
    print(f"  {name}: {before:,} → {len(df):,} after dedup on UniProt ID")
    return df


# ── Remove cross-contamination ────────────────────────────────────────────────

def remove_overlap(df_pos: pd.DataFrame, df_neg: pd.DataFrame):
    """Remove any UniProt ID that appears in both sets."""
    pos_uids = set(df_pos["uniprot_id"])
    before = len(df_neg)
    df_neg = df_neg[~df_neg["uniprot_id"].isin(pos_uids)].reset_index(drop=True)
    removed = before - len(df_neg)
    if removed:
        print(f"  Removed {removed:,} negatives that also appear in positives")
    return df_neg


# ── Attach Km to positives ────────────────────────────────────────────────────

def attach_km(df_pos: pd.DataFrame, df_km: pd.DataFrame) -> pd.DataFrame:
    df_km = df_km[["uniprot_id", "km_best_mM"]].copy()
    df_km["km_best_mM"] = pd.to_numeric(df_km["km_best_mM"], errors="coerce")
    df_km = df_km[df_km["km_best_mM"] > 0]

    df_pos = df_pos.merge(df_km, on="uniprot_id", how="left")

    n_km = df_pos["km_best_mM"].notna().sum()
    print(f"  Positives with experimental Km: {n_km:,} / {len(df_pos):,}")

    # log10 Km in mM (for regression)
    df_pos["km_log10_mM"] = np.where(
        df_pos["km_best_mM"].notna() & (df_pos["km_best_mM"] > 0),
        np.log10(df_pos["km_best_mM"].astype(float)),
        np.nan,
    )
    return df_pos


# ── Sample negatives ──────────────────────────────────────────────────────────

def sample_negatives(df_pos: pd.DataFrame, df_neg: pd.DataFrame,
                     ratio: float = 5.0) -> pd.DataFrame:
    """
    Sample negatives to ratio × n_positives.
    Keep all if negatives are fewer than target.
    Sample stratified by EC class if possible.
    """
    target = int(len(df_pos) * ratio)
    if len(df_neg) <= target:
        print(f"  Keeping all {len(df_neg):,} negatives (< target {target:,})")
        return df_neg
    # Stratified: proportional to EC class size
    sampled = (
        df_neg.groupby("ec_number", group_keys=False)
        .apply(lambda g: g.sample(
            min(len(g), max(1, int(len(g) / len(df_neg) * target))),
            random_state=42
        ))
        .reset_index(drop=True)
    )
    # Top up if rounding left us short
    if len(sampled) < target:
        extra = df_neg[~df_neg["uniprot_id"].isin(sampled["uniprot_id"])]
        extra = extra.sample(min(len(extra), target - len(sampled)), random_state=42)
        sampled = pd.concat([sampled, extra], ignore_index=True)
    print(f"  Sampled {len(sampled):,} negatives from {len(df_neg):,} "
          f"(ratio {len(sampled)/len(df_pos):.1f}×)")
    return sampled


# ── Write outputs ─────────────────────────────────────────────────────────────

def write_fasta(df: pd.DataFrame, path: Path):
    with open(path, "w") as f:
        for _, row in df.iterrows():
            label = int(row.get("label", -1))
            km_tag = (f"|km={row['km_best_mM']:.6f}mM"
                      if pd.notna(row.get("km_best_mM")) else "")
            f.write(
                f">{row['uniprot_id']}|{row['ec_number']}"
                f"|label={label}{km_tag}\n"
            )
            seq = str(row["sequence"])
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Step 02 — Prepare Training Data")
    print("=" * 70)

    df_pos_raw, df_neg_raw, df_km = load_raw()

    print("\n── Validate sequences ──")
    df_pos = validate(df_pos_raw, "Positives")
    df_neg = validate(df_neg_raw, "Negatives")

    print("\n── Deduplicate ──")
    df_pos = dedup(df_pos, "Positives")
    df_neg = dedup(df_neg, "Negatives")

    print("\n── Remove overlap ──")
    df_neg = remove_overlap(df_pos, df_neg)

    print("\n── Attach Km values ──")
    df_pos = attach_km(df_pos, df_km)

    print("\n── Sample negatives (5× positives) ──")
    df_neg_sampled = sample_negatives(df_pos, df_neg, ratio=5.0)

    # Assign labels
    df_pos["label"] = 1
    df_neg_sampled["label"] = 0

    # Combined training set
    KEEP_COLS = [
        "uniprot_id", "ec_number", "enzyme_name", "organism",
        "sequence", "length", "label", "source",
        "km_best_mM", "km_log10_mM",
    ]
    for col in KEEP_COLS:
        if col not in df_pos.columns:
            df_pos[col] = None
        if col not in df_neg_sampled.columns:
            df_neg_sampled[col] = None

    df_train = pd.concat(
        [df_pos[KEEP_COLS], df_neg_sampled[KEEP_COLS]],
        ignore_index=True
    ).sample(frac=1, random_state=42).reset_index(drop=True)   # shuffle

    # ── Save ──────────────────────────────────────────────────────────────────
    train_tsv = OUT_DIR / "training_data.tsv"
    df_train.to_csv(train_tsv, sep="\t", index=False)
    print(f"\n✓ training_data.tsv:    {len(df_train):,} rows")

    write_fasta(df_train, OUT_DIR / "training_data.fasta")
    print(f"✓ training_data.fasta:  {len(df_train):,} sequences")

    # Km-only subset (for regression training)
    df_km_train = df_train[df_train["km_best_mM"].notna()].copy()
    df_km_train.to_csv(OUT_DIR / "km_training.tsv", sep="\t", index=False)
    print(f"✓ km_training.tsv:      {len(df_km_train):,} rows with experimental Km")

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("TRAINING DATA SUMMARY")
    print("=" * 70)
    print(f"  Total rows:          {len(df_train):,}")
    print(f"  Positives:           {(df_train['label']==1).sum():,}")
    print(f"  Negatives:           {(df_train['label']==0).sum():,}")
    print(f"  With Km (positives): {df_km_train['label'].eq(1).sum():,}")
    print(f"\n  Positive EC classes:")
    pos_only = df_train[df_train["label"] == 1]
    for ec, cnt in pos_only["ec_number"].value_counts().items():
        km_n = pos_only[pos_only["ec_number"] == ec]["km_best_mM"].notna().sum()
        print(f"    {ec}: {cnt:>5,} seqs  ({km_n} with Km)")
    print(f"\n  Km range (positives): "
          f"{df_km_train['km_best_mM'].min():.4f} – "
          f"{df_km_train['km_best_mM'].max():.2f} mM")
    print(f"\n  Output → {OUT_DIR}/")
    print("\n  Next step: python scripts/03_annotate_features.py")


if __name__ == "__main__":
    main()
