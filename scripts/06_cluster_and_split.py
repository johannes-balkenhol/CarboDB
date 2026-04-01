#!/usr/bin/env python3
"""
06_cluster_and_split.py
=======================
CarboDB — Step 06: Sequence clustering + stratified train/val/test split.

Two-level split strategy:
  1. CD-HIT 40% identity clustering → prevents sequence leakage between splits
  2. Stratified assignment → ensures EC class + Km range + taxonomy diversity

Creates three datasets:
  - Binary classification: label=1 vs label=0 (excludes label=2)
  - EC class prediction:   label=1 only, tier 1+2 evidence, 40 EC classes
  - Km regression:         label=1 with BRENDA experimental Km (2,971 seqs)

Output:
  data/splits/cluster_assignments.tsv  — cdb_id → cluster_id
  data/splits/split_binary.tsv         — cdb_id, label, split (train/val/test)
  data/splits/split_ec.tsv             — cdb_id, ec_number, split
  data/splits/split_km.tsv             — cdb_id, km_log10_mM, split
  data/splits/split_summary.json       — counts per split per task

Usage:
  python scripts/06_cluster_and_split.py
  python scripts/06_cluster_and_split.py --identity 0.4  # default
  python scripts/06_cluster_and_split.py --dry-run
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("06_cluster_and_split")

SPLIT_DIR  = PATHS.PRIMARY.parent / "splits"
FASTA_ALL  = PATHS.PRIMARY / "master.fasta"

# Split ratios
TRAIN_RATIO = 0.80
VAL_RATIO   = 0.10
TEST_RATIO  = 0.10

# Km bins for stratification (log10 mM)
KM_BINS = [-4, -1, 0, 2]   # <0.1mM | 0.1-1mM | >1mM
KM_LABELS = ["low", "medium", "high"]


# ═══════════════════════════════════════════════════════════════════════════
# Step 1: Run CD-HIT
# ═══════════════════════════════════════════════════════════════════════════

def run_cdhit(fasta_path: Path, identity: float, threads: int, dry_run: bool) -> Path:
    """Run CD-HIT on all sequences. Returns path to cluster file."""
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)
    out_prefix = SPLIT_DIR / f"cdhit_{int(identity*100)}"
    clstr_file = Path(str(out_prefix) + ".clstr")

    if clstr_file.exists():
        log.info("CD-HIT cluster file exists — skipping: %s", clstr_file)
        return clstr_file

    # Choose word length based on identity
    # CD-HIT recommendation: id>0.7→5, id>0.6→4, id>0.5→3, id>0.4→2
    word_len = 2 if identity <= 0.4 else 3 if identity <= 0.5 else 4

    cmd = [
        "cd-hit",
        "-i", str(fasta_path),
        "-o", str(out_prefix),
        "-c", str(identity),
        "-n", str(word_len),
        "-T", str(threads),
        "-M", "80000",        # 80GB RAM
        "-d", "0",            # full sequence name in output
        "-g", "1",            # accurate mode
        "-aS", "0.8",         # alignment coverage for shorter seq
        "-sc", "1",           # sort clusters by size
    ]

    log.info("Running CD-HIT: identity=%.2f, word_len=%d", identity, word_len)
    log.info("Command: %s", " ".join(cmd))

    if dry_run:
        log.info("DRY RUN — skipping CD-HIT execution")
        return clstr_file

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.error("CD-HIT failed: %s", result.stderr)
        raise RuntimeError("CD-HIT failed")

    log.info("CD-HIT complete — cluster file: %s", clstr_file)
    return clstr_file


# ═══════════════════════════════════════════════════════════════════════════
# Step 2: Parse CD-HIT cluster file
# ═══════════════════════════════════════════════════════════════════════════

def parse_cdhit_clusters(clstr_file: Path) -> pd.DataFrame:
    """Parse CD-HIT .clstr file → DataFrame with cdb_id, cluster_id, is_representative."""
    log.info("Parsing CD-HIT clusters from %s", clstr_file)

    if not clstr_file.exists():
        log.warning("Cluster file not found — creating dummy clusters")
        return pd.DataFrame(columns=["cdb_id", "cluster_id", "is_representative"])

    rows = []
    cluster_id = -1

    with open(clstr_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[1])
            elif line:
                # Format: 0   476aa, >CDB000001... *
                # or:     1   476aa, >CDB000002... at 95.23%
                is_rep = line.endswith("*")
                parts = line.split(">")
                if len(parts) < 2:
                    continue
                cdb_id = parts[1].split("...")[0].split("|")[0].strip()
                rows.append({
                    "cdb_id":           cdb_id,
                    "cluster_id":       cluster_id,
                    "is_representative": is_rep,
                })

    df = pd.DataFrame(rows)
    log.info("  Parsed %d sequences in %d clusters",
             len(df), df["cluster_id"].nunique())
    log.info("  Cluster size distribution:")
    sizes = df.groupby("cluster_id").size()
    log.info("    median=%d  mean=%.1f  max=%d  singletons=%d",
             sizes.median(), sizes.mean(), sizes.max(), (sizes == 1).sum())

    return df


# ═══════════════════════════════════════════════════════════════════════════
# Step 3: Load master + assign strata
# ═══════════════════════════════════════════════════════════════════════════

def load_master_with_strata() -> pd.DataFrame:
    """Load master.tsv and assign stratification keys."""
    log.info("Loading master.tsv...")
    df = pd.read_csv(PATHS.MASTER_TSV, sep="\t", dtype=str,
                     usecols=["cdb_id", "uniprot_id", "ec_number", "label",
                               "source", "organism", "km_best_mM", "km_log10_mM"])

    df["label"]       = df["label"].astype(int)
    df["km_log10_mM"] = pd.to_numeric(df["km_log10_mM"], errors="coerce")
    df["km_best_mM"]  = pd.to_numeric(df["km_best_mM"],  errors="coerce")

    # Taxonomy group: extract kingdom from organism name (simple heuristic)
    df["kingdom"] = df["organism"].apply(_infer_kingdom)

    # Km bin for stratification
    df["km_bin"] = pd.cut(
        df["km_log10_mM"],
        bins=KM_BINS,
        labels=KM_LABELS,
        include_lowest=True
    ).astype(str)
    df["km_bin"] = df["km_bin"].replace("nan", "none")

    # Stratification key for binary task: ec_number + kingdom
    df["strat_binary"] = df["ec_number"].fillna("unk") + "|" + df["kingdom"]

    # Stratification key for EC task: ec_number
    df["strat_ec"] = df["ec_number"].fillna("unk")

    # Stratification key for Km task: ec_number + km_bin
    df["strat_km"] = df["ec_number"].fillna("unk") + "|" + df["km_bin"]

    log.info("  Loaded %d sequences", len(df))
    log.info("  label=0: %d  label=1: %d  label=2: %d",
             (df["label"]==0).sum(), (df["label"]==1).sum(), (df["label"]==2).sum())
    log.info("  With Km: %d", df["km_best_mM"].notna().sum())

    return df


def _infer_kingdom(organism: str) -> str:
    """Simple organism → kingdom mapping."""
    if not isinstance(organism, str):
        return "unknown"
    org = organism.lower()
    if any(x in org for x in ["homo", "mus ", "rat ", "bos ", "sus ", "gallus"]):
        return "vertebrate"
    if any(x in org for x in ["arabidopsis", "oryza", "zea ", "solanum", "nicotiana",
                                "spinacia", "triticum", "glycine"]):
        return "plant"
    if any(x in org for x in ["saccharomyces", "aspergillus", "candida", "neurospora"]):
        return "fungi"
    if any(x in org for x in ["escherichia", "bacillus", "pseudomonas", "streptomyces",
                                "mycobacterium", "staphylococcus", "salmonella"]):
        return "bacteria"
    if any(x in org for x in ["methan", "sulfolobus", "pyrococcus", "halobacter",
                                "archaeo", "thermo"]):
        return "archaea"
    return "other"


# ═══════════════════════════════════════════════════════════════════════════
# Step 4: Cluster-based stratified split
# ═══════════════════════════════════════════════════════════════════════════

def cluster_stratified_split(
    df: pd.DataFrame,
    clusters: pd.DataFrame,
    strat_col: str,
    seed: int = 42
) -> pd.Series:
    """
    Assign train/val/test splits at the CLUSTER level.
    All members of a cluster go to the same split.
    Stratified by strat_col to ensure diversity.
    Returns Series of split assignments indexed by cdb_id.
    """
    # Merge clusters with master
    merged = df.merge(clusters[["cdb_id", "cluster_id"]], on="cdb_id", how="left")

    # For sequences without cluster assignment, assign to own cluster
    no_cluster = merged["cluster_id"].isna()
    if no_cluster.sum() > 0:
        log.warning("  %d sequences not in cluster file — assigning singleton clusters",
                    no_cluster.sum())
        max_id = merged["cluster_id"].max() if not merged["cluster_id"].isna().all() else 0
        merged.loc[no_cluster, "cluster_id"] = (
            range(int(max_id) + 1, int(max_id) + 1 + no_cluster.sum())
        )
    merged["cluster_id"] = merged["cluster_id"].astype(int)

    # Get cluster-level representative strat label (majority vote within cluster)
    cluster_strat = (
        merged.groupby("cluster_id")[strat_col]
        .agg(lambda x: x.value_counts().index[0])
        .reset_index()
    )

    # Assign clusters to splits
    rng = np.random.default_rng(seed)
    cluster_strat = cluster_strat.sample(frac=1, random_state=seed).reset_index(drop=True)

    # Stratified split at cluster level
    split_assignments = {}
    strat_groups = cluster_strat.groupby(strat_col)

    for strat, group in strat_groups:
        n = len(group)
        n_test = max(1, int(n * TEST_RATIO))
        n_val  = max(1, int(n * VAL_RATIO))
        n_train = n - n_test - n_val

        idx = group["cluster_id"].values
        rng.shuffle(idx)

        for c in idx[:n_train]:
            split_assignments[c] = "train"
        for c in idx[n_train:n_train+n_val]:
            split_assignments[c] = "val"
        for c in idx[n_train+n_val:]:
            split_assignments[c] = "test"

    # Map back to sequences
    merged["split"] = merged["cluster_id"].map(split_assignments).fillna("train")

    return merged.set_index("cdb_id")["split"]


# ═══════════════════════════════════════════════════════════════════════════
# Step 5: Build and save splits
# ═══════════════════════════════════════════════════════════════════════════

def build_splits(df: pd.DataFrame, clusters: pd.DataFrame, dry_run: bool):
    SPLIT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Binary classification split ──────────────────────────────────────
    log.info("── Binary classification split ──")
    binary_df = df[df["label"].isin([0, 1])].copy()
    binary_splits = cluster_stratified_split(binary_df, clusters, "strat_binary")
    binary_df["split"] = binary_df["cdb_id"].map(binary_splits)

    log.info("  Train: %d  Val: %d  Test: %d",
             (binary_df["split"]=="train").sum(),
             (binary_df["split"]=="val").sum(),
             (binary_df["split"]=="test").sum())

    # Check balance in test set
    test_binary = binary_df[binary_df["split"] == "test"]
    log.info("  Test positives: %d  Test negatives: %d",
             (test_binary["label"]==1).sum(), (test_binary["label"]==0).sum())

    # ── EC class prediction split ────────────────────────────────────────
    log.info("── EC class prediction split ──")
    ec_df = df[(df["label"] == 1)].copy()
    # Filter to tier 1+2 evidence only (BRENDA + SwissProt)
    tier12_sources = {"brenda", "swissprot"}
    ec_df = ec_df[ec_df["source"].isin(tier12_sources)]
    log.info("  EC dataset: %d sequences (tier 1+2 only)", len(ec_df))

    # Check EC class coverage
    ec_counts = ec_df["ec_number"].value_counts()
    log.info("  EC classes with >10 seqs: %d", (ec_counts > 10).sum())
    log.info("  EC classes with >100 seqs: %d", (ec_counts > 100).sum())

    ec_splits = cluster_stratified_split(ec_df, clusters, "strat_ec")
    ec_df["split"] = ec_df["cdb_id"].map(ec_splits)

    log.info("  Train: %d  Val: %d  Test: %d",
             (ec_df["split"]=="train").sum(),
             (ec_df["split"]=="val").sum(),
             (ec_df["split"]=="test").sum())

    # ── Km regression split ──────────────────────────────────────────────
    log.info("── Km regression split ──")
    km_df = df[
        (df["label"] == 1) &
        (df["km_best_mM"].notna()) &
        (df["source"] == "brenda")
    ].copy()
    log.info("  Km dataset: %d sequences", len(km_df))

    # Log Km distribution
    log.info("  Km log10 range: %.2f to %.2f",
             km_df["km_log10_mM"].min(), km_df["km_log10_mM"].max())
    log.info("  Km bins: %s",
             km_df["km_bin"].value_counts().to_dict())

    km_splits = cluster_stratified_split(km_df, clusters, "strat_km")
    km_df["split"] = km_df["cdb_id"].map(km_splits)

    log.info("  Train: %d  Val: %d  Test: %d",
             (km_df["split"]=="train").sum(),
             (km_df["split"]=="val").sum(),
             (km_df["split"]=="test").sum())

    # Check Km distribution in each split
    for split in ["train", "val", "test"]:
        sub = km_df[km_df["split"] == split]
        log.info("  %s — mean log10Km=%.2f  std=%.2f  n=%d",
                 split, sub["km_log10_mM"].mean(), sub["km_log10_mM"].std(), len(sub))

    if dry_run:
        log.info("DRY RUN — not writing split files")
        return

    # ── Save split files ─────────────────────────────────────────────────
    binary_df[["cdb_id", "uniprot_id", "ec_number", "label",
               "source", "kingdom", "split"]].to_csv(
        SPLIT_DIR / "split_binary.tsv", sep="\t", index=False
    )
    log.info("Saved: %s", SPLIT_DIR / "split_binary.tsv")

    ec_df[["cdb_id", "uniprot_id", "ec_number", "label",
           "source", "kingdom", "split"]].to_csv(
        SPLIT_DIR / "split_ec.tsv", sep="\t", index=False
    )
    log.info("Saved: %s", SPLIT_DIR / "split_ec.tsv")

    km_df[["cdb_id", "uniprot_id", "ec_number", "km_best_mM",
           "km_log10_mM", "km_bin", "split"]].to_csv(
        SPLIT_DIR / "split_km.tsv", sep="\t", index=False
    )
    log.info("Saved: %s", SPLIT_DIR / "split_km.tsv")

    # ── Save summary ─────────────────────────────────────────────────────
    summary = {
        "created_at": TS,
        "clustering_identity": 0.4,
        "n_clusters": int(clusters["cluster_id"].nunique()),
        "binary": {
            "total": len(binary_df),
            "train": int((binary_df["split"]=="train").sum()),
            "val":   int((binary_df["split"]=="val").sum()),
            "test":  int((binary_df["split"]=="test").sum()),
            "train_pos": int(((binary_df["split"]=="train") & (binary_df["label"]==1)).sum()),
            "train_neg": int(((binary_df["split"]=="train") & (binary_df["label"]==0)).sum()),
            "test_pos":  int(((binary_df["split"]=="test") & (binary_df["label"]==1)).sum()),
            "test_neg":  int(((binary_df["split"]=="test") & (binary_df["label"]==0)).sum()),
        },
        "ec_prediction": {
            "total":      len(ec_df),
            "n_ec_classes": int(ec_df["ec_number"].nunique()),
            "train": int((ec_df["split"]=="train").sum()),
            "val":   int((ec_df["split"]=="val").sum()),
            "test":  int((ec_df["split"]=="test").sum()),
        },
        "km_regression": {
            "total":  len(km_df),
            "train":  int((km_df["split"]=="train").sum()),
            "val":    int((km_df["split"]=="val").sum()),
            "test":   int((km_df["split"]=="test").sum()),
            "km_range_mM": [
                float(km_df["km_best_mM"].min()),
                float(km_df["km_best_mM"].max())
            ],
            "km_log10_std": float(km_df["km_log10_mM"].std()),
        },
    }

    with open(SPLIT_DIR / "split_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    log.info("Saved: %s", SPLIT_DIR / "split_summary.json")

    return summary


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity",  type=float, default=0.4,
                    help="CD-HIT sequence identity threshold (default: 0.4)")
    ap.add_argument("--threads",   type=int,   default=80)
    ap.add_argument("--dry-run",   action="store_true")
    ap.add_argument("--skip-cdhit", action="store_true",
                    help="Skip CD-HIT if cluster file already exists")
    args = ap.parse_args()

    log.info("Step 06: Clustering + split  (identity=%.2f)", args.identity)

    # Check CD-HIT available
    if not shutil.which("cd-hit"):
        log.error("cd-hit not found — install via: conda install -c bioconda cd-hit")
        sys.exit(1)

    # Run CD-HIT
    clstr_file = run_cdhit(FASTA_ALL, args.identity, args.threads, args.dry_run)

    # Parse clusters
    clusters = parse_cdhit_clusters(clstr_file)

    if not args.dry_run:
        clusters.to_csv(SPLIT_DIR / "cluster_assignments.tsv", sep="\t", index=False)
        log.info("Saved cluster assignments: %s",
                 SPLIT_DIR / "cluster_assignments.tsv")

    # Load master
    df = load_master_with_strata()

    # Build splits
    summary = build_splits(df, clusters, args.dry_run)

    if summary:
        log.info("\n" + "="*60)
        log.info("SPLIT SUMMARY")
        log.info("="*60)
        log.info("Binary:  train=%d  val=%d  test=%d",
                 summary["binary"]["train"],
                 summary["binary"]["val"],
                 summary["binary"]["test"])
        log.info("EC pred: train=%d  val=%d  test=%d",
                 summary["ec_prediction"]["train"],
                 summary["ec_prediction"]["val"],
                 summary["ec_prediction"]["test"])
        log.info("Km regr: train=%d  val=%d  test=%d",
                 summary["km_regression"]["train"],
                 summary["km_regression"]["val"],
                 summary["km_regression"]["test"])

    log.info("Done. Next: python scripts/07_build_feature_matrix.py")


if __name__ == "__main__":
    main()
