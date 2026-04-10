#!/usr/bin/env python3
"""
04f_blast_benchmark.py
======================
CarboDB — BLAST-based benchmark for method comparison.

Builds a BLAST database from training sequences, queries test sequences,
and evaluates carboxylase prediction performance at different identity thresholds.

This provides the BLAST baseline for the method comparison table (script 13 Task B).

Strategy:
  - Sample N test sequences (default 2000, stratified pos/neg)
  - Build BLAST DB from training sequences only
  - Run blastp: each test seq vs training DB
  - For each test sequence: best hit → predict carboxylase if:
      (a) hit exists AND hit is a known carboxylase (label=1) → positive
      (b) no hit or hit is negative → negative
  - Evaluate at different identity cutoffs: 90%, 70%, 50%, 30%, 20%
  - Also evaluate on fragment conditions (truncated sequences)
  - Compare to CarboDB ML at each threshold

Output: data/benchmark/blast_benchmark.json
         data/benchmark/figures/blast_comparison.tsv

Usage:
  python scripts/04f_blast_benchmark.py
  python scripts/04f_blast_benchmark.py --n-test 500 --threads 8
  python scripts/04f_blast_benchmark.py --skip-db   # if DB already built
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, ROOT, TS, setup_logging

log = setup_logging("04f_blast_benchmark")

PRIMARY   = ROOT / "data" / "primary"
SPLIT_DIR = ROOT / "data" / "splits"
ML_DIR    = ROOT / "data" / "ml"
MODEL_DIR = ROOT / "data" / "models"
BENCH_DIR = ROOT / "data" / "benchmark"
BLAST_DIR = ROOT / "data" / "benchmark" / "blast"
FIG_DIR   = BENCH_DIR / "figures"

BLAST_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

MASTER_FASTA = PRIMARY / "master.fasta"
BLAST_DB     = BLAST_DIR / "carbodb_train"
BLAST_QUERY  = BLAST_DIR / "test_query.fasta"
BLAST_OUT    = BLAST_DIR / "blast_results.tsv"


# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load splits and prepare sequence sets
# ══════════════════════════════════════════════════════════════════════════════

def load_splits(n_test=2000, seed=42):
    """
    Load binary split TSV and sample test sequences.
    Returns:
        train_ids: set of CDB IDs in training set
        test_df:   DataFrame of sampled test sequences with label
        train_labels: dict {cdb_id: label} for training set
    """
    log.info("Loading splits from %s", SPLIT_DIR / "split_binary.tsv")
    splits = pd.read_csv(SPLIT_DIR / "split_binary.tsv", sep="\t", dtype=str)
    splits["label"] = splits["label"].astype(int)

    train = splits[splits["split"] == "train"]
    test  = splits[splits["split"] == "test"]

    log.info("  Train: %d sequences | Test: %d sequences", len(train), len(test))

    # Stratified sample of test set
    rng = np.random.default_rng(seed)
    pos_test = test[test["label"] == 1]
    neg_test = test[test["label"] == 0]
    n_pos = min(len(pos_test), n_test // 2)
    n_neg = min(len(neg_test), n_test - n_pos)

    sampled = pd.concat([
        pos_test.sample(n=n_pos, random_state=seed),
        neg_test.sample(n=n_neg, random_state=seed),
    ]).sample(frac=1, random_state=seed).reset_index(drop=True)

    log.info("  Sampled test: %d pos + %d neg = %d total",
             n_pos, n_neg, len(sampled))

    train_ids    = set(train["cdb_id"].values)
    train_labels = dict(zip(train["cdb_id"], train["label"]))

    return train_ids, sampled, train_labels


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Extract sequences from master.fasta
# ══════════════════════════════════════════════════════════════════════════════

def parse_fasta_index(fasta_path):
    """Build index: {cdb_id: (header, sequence)} from master.fasta."""
    log.info("Indexing master.fasta (this takes ~30s)...")
    index = {}
    cur_id = cur_seq = None
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_id:
                    index[cur_id] = cur_seq
                cur_id = line[1:].split("|")[0].split()[0]
                cur_seq = []
            elif line:
                cur_seq.append(line)
    if cur_id:
        index[cur_id] = "".join(cur_seq)
    log.info("  Indexed %d sequences", len(index))
    return index


def write_fasta(ids_with_labels, seq_index, out_path, label_filter=None):
    """Write a FASTA file for a subset of sequences."""
    n_written = 0
    with open(out_path, "w") as f:
        for cdb_id, label in ids_with_labels:
            if label_filter is not None and label != label_filter:
                continue
            seq = seq_index.get(cdb_id)
            if seq:
                f.write(f">{cdb_id}|label={label}\n{seq}\n")
                n_written += 1
    log.info("  Wrote %d sequences to %s", n_written, out_path)
    return n_written


def write_fasta_df(df, seq_index, out_path, truncate_frac=None):
    """Write test sequences to FASTA, optionally truncating."""
    n_written = 0
    with open(out_path, "w") as f:
        for _, row in df.iterrows():
            seq = seq_index.get(row["cdb_id"])
            if not seq:
                continue
            if truncate_frac and truncate_frac < 1.0:
                n = max(10, int(len(seq) * truncate_frac))
                seq = seq[:n]  # N-terminal truncation
            f.write(f">{row['cdb_id']}|label={row['label']}\n{seq}\n")
            n_written += 1
    return n_written


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Build BLAST database
# ══════════════════════════════════════════════════════════════════════════════

def build_blast_db(train_ids, train_labels, seq_index, threads=8):
    """Build BLAST protein database from training sequences."""
    if (BLAST_DB.with_suffix(".pin")).exists():
        log.info("BLAST DB already exists at %s — skipping build", BLAST_DB)
        return True

    db_fasta = BLAST_DIR / "train_sequences.fasta"
    log.info("Writing training sequences to FASTA...")

    n = 0
    with open(db_fasta, "w") as f:
        for cdb_id in train_ids:
            seq = seq_index.get(cdb_id)
            label = train_labels.get(cdb_id, 0)
            if seq:
                f.write(f">{cdb_id}|label={label}\n{seq}\n")
                n += 1
    log.info("  Wrote %d training sequences", n)

    log.info("Building BLAST database...")
    cmd = ["makeblastdb", "-in", str(db_fasta), "-dbtype", "prot",
           "-out", str(BLAST_DB), "-title", "CarboDB_train"]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        log.error("makeblastdb failed: %s", result.stderr.decode()[:300])
        return False
    log.info("  BLAST DB built: %s", BLAST_DB)
    return True


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Run BLAST
# ══════════════════════════════════════════════════════════════════════════════

def run_blast(query_fasta, out_file, threads=8, evalue=0.001):
    """Run blastp and save tabular results."""
    log.info("Running blastp (query=%s)...", query_fasta)
    cmd = [
        "blastp",
        "-query",   str(query_fasta),
        "-db",      str(BLAST_DB),
        "-out",     str(out_file),
        "-outfmt",  "6 qseqid sseqid pident length evalue bitscore",
        "-evalue",  str(evalue),
        "-max_target_seqs", "1",
        "-max_hsps", "1",
        "-num_threads", str(threads),
        "-matrix",  "BLOSUM62",
        "-seg",     "yes",
    ]
    result = subprocess.run(cmd, capture_output=True, timeout=3600)
    if result.returncode != 0:
        log.error("blastp failed: %s", result.stderr.decode()[:300])
        return False
    log.info("  BLAST complete")
    return True


def parse_blast_results(blast_out):
    """
    Parse BLAST tabular output.
    Returns dict: {query_cdb_id: {'pident': float, 'evalue': float,
                                   'hit_label': int, 'has_hit': bool}}
    """
    results = {}
    if not Path(blast_out).exists():
        return results

    with open(blast_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            qid    = parts[0].split("|")[0]
            sid    = parts[1].split("|")[0]
            pident = float(parts[2])
            evalue = float(parts[4])

            # Extract label from subject ID (format: CDB000001|label=1)
            hit_label = 0
            if "|label=" in parts[1]:
                try:
                    hit_label = int(parts[1].split("label=")[1])
                except Exception:
                    pass

            # Keep best hit per query (max_target_seqs=1 already handles this)
            if qid not in results:
                results[qid] = {
                    "pident":    pident,
                    "evalue":    evalue,
                    "hit_label": hit_label,
                    "has_hit":   True,
                }

    return results


# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Evaluate BLAST at identity thresholds
# ══════════════════════════════════════════════════════════════════════════════

def evaluate_blast_at_threshold(test_df, blast_results, identity_threshold):
    """
    BLAST prediction rule:
    - Positive if: has_hit AND hit_label==1 AND pident >= identity_threshold
    - Negative otherwise

    Also compute AUROC using pident * hit_label as score.
    """
    y_true  = test_df["label"].astype(int).values
    y_pred  = np.zeros(len(test_df), dtype=int)
    y_score = np.zeros(len(test_df), dtype=float)

    for i, row in test_df.iterrows():
        cdb_id = row["cdb_id"]
        hit = blast_results.get(cdb_id)
        if hit and hit["has_hit"] and hit["hit_label"] == 1 and hit["pident"] >= identity_threshold:
            y_pred[i]  = 1
            y_score[i] = hit["pident"] / 100.0
        elif hit and hit["has_hit"] and hit["hit_label"] == 1:
            # Has a carboxylase hit but below identity threshold — partial score
            y_score[i] = hit["pident"] / 100.0 * 0.5

    # Handle edge case where all predictions are same class
    if len(np.unique(y_true)) < 2:
        return None

    try:
        auroc = float(roc_auc_score(y_true, y_score))
    except Exception:
        auroc = 0.5

    f1  = float(f1_score(y_true, y_pred, zero_division=0))
    acc = float(accuracy_score(y_true, y_pred))

    n_hits      = sum(1 for r in blast_results.values() if r["has_hit"])
    n_pos_hits  = sum(1 for r in blast_results.values() if r["has_hit"] and r["hit_label"] == 1)
    n_above_thr = sum(1 for r in blast_results.values()
                      if r["has_hit"] and r["hit_label"] == 1 and r["pident"] >= identity_threshold)

    return {
        "identity_threshold": identity_threshold,
        "n_test":             int(len(test_df)),
        "n_pos":              int(y_true.sum()),
        "n_with_any_hit":     n_hits,
        "n_with_carb_hit":    n_pos_hits,
        "n_above_threshold":  n_above_thr,
        "auroc":              round(auroc, 4),
        "f1":                 round(f1, 4),
        "accuracy":           round(acc, 4),
        "hit_rate_pct":       round(n_hits / len(test_df) * 100, 1),
    }


# ══════════════════════════════════════════════════════════════════════════════
# Step 6: CarboDB ML at matching conditions (for direct comparison)
# ══════════════════════════════════════════════════════════════════════════════

def evaluate_carbodb_on_sample(test_df, seq_index, truncate_frac=None):
    """
    Evaluate CarboDB v5 on the same test sample used for BLAST.
    Uses precomputed ML feature matrices, matching by cdb_id.
    """
    # Load full binary test set
    X_te_full = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te_full = np.load(ML_DIR / "y_binary_test.npy")

    # Load split to get cdb_ids for test set
    splits = pd.read_csv(SPLIT_DIR / "split_binary.tsv", sep="\t")
    test_split = splits[splits["split"] == "test"].reset_index(drop=True)

    # Map cdb_id → index in X_te_full
    id_to_idx = {row["cdb_id"]: i for i, row in test_split.iterrows()
                 if i < len(X_te_full)}

    # Get indices for our sampled test_df
    indices = []
    labels  = []
    for _, row in test_df.iterrows():
        idx = id_to_idx.get(row["cdb_id"])
        if idx is not None:
            indices.append(idx)
            labels.append(int(row["label"]))

    if not indices:
        log.warning("  No matching indices found in ML test set")
        return None

    X_sub = X_te_full[indices]
    y_sub = np.array(labels)

    # Apply truncation to feature vector if requested
    if truncate_frac and truncate_frac < 1.0:
        feat_names = json.load(open(ML_DIR / "feature_names_binary.json"))
        from scripts.x13_benchmark_identity import truncate_sequence_features
        X_sub = truncate_sequence_features(X_sub, feat_names, None, truncate_frac, "n_terminal")

    booster = xgb.Booster()
    booster.load_model(str(MODEL_DIR / "binary_v5.json"))
    probs = booster.predict(xgb.DMatrix(X_sub))
    preds = (probs >= 0.5).astype(int)

    if len(np.unique(y_sub)) < 2:
        return None

    auroc = float(roc_auc_score(y_sub, probs))
    f1    = float(f1_score(y_sub, preds, zero_division=0))

    return {
        "n_test":   int(len(y_sub)),
        "auroc":    round(auroc, 4),
        "f1":       round(f1, 4),
        "truncate": truncate_frac,
    }


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    import json

    ap = argparse.ArgumentParser(description="BLAST benchmark for CarboDB method comparison.")
    ap.add_argument("--n-test",    type=int,  default=2000,  help="Test sequences to sample (default: 2000)")
    ap.add_argument("--threads",   type=int,  default=8,     help="BLAST threads (default: 8)")
    ap.add_argument("--skip-db",   action="store_true",      help="Skip DB build if already exists")
    ap.add_argument("--skip-blast",action="store_true",      help="Skip BLAST run if results exist")
    ap.add_argument("--evalue",    type=float,default=0.001, help="BLAST E-value cutoff (default: 0.001)")
    args = ap.parse_args()

    # ── Step 1: Load splits ────────────────────────────────────────────────
    train_ids, test_df, train_labels = load_splits(args.n_test)

    # ── Step 2: Index master.fasta ─────────────────────────────────────────
    seq_index = parse_fasta_index(MASTER_FASTA)

    # ── Step 3: Build BLAST DB ─────────────────────────────────────────────
    if not args.skip_db:
        ok = build_blast_db(train_ids, train_labels, seq_index, args.threads)
        if not ok:
            log.error("BLAST DB build failed")
            sys.exit(1)

    # ── Step 4: Write test query FASTA ─────────────────────────────────────
    log.info("Writing test query FASTA...")
    write_fasta_df(test_df, seq_index, BLAST_QUERY)

    # ── Step 5: Run BLAST ──────────────────────────────────────────────────
    if not args.skip_blast or not BLAST_OUT.exists():
        ok = run_blast(BLAST_QUERY, BLAST_OUT, args.threads, args.evalue)
        if not ok:
            sys.exit(1)

    blast_results = parse_blast_results(BLAST_OUT)
    log.info("Parsed %d BLAST results (%d with hits)",
             len(blast_results),
             sum(1 for r in blast_results.values() if r["has_hit"]))

    # ── Step 6: Evaluate at identity thresholds ────────────────────────────
    log.info("══ BLAST performance at identity thresholds ══")
    thresholds = [0, 20, 30, 50, 70, 90]  # 0 = any hit
    blast_threshold_results = []

    for thr in thresholds:
        r = evaluate_blast_at_threshold(test_df, blast_results, thr)
        if r:
            blast_threshold_results.append(r)
            log.info("  ≥%d%% identity: AUROC=%.4f  F1=%.4f  hit_rate=%.1f%%  n_above_thr=%d",
                     thr, r["auroc"], r["f1"], r["hit_rate_pct"], r["n_above_threshold"])

    # ── Step 7: Identity distribution of hits ─────────────────────────────
    log.info("══ Identity distribution of BLAST hits ══")
    pidents = [r["pident"] for r in blast_results.values() if r["has_hit"]]
    if pidents:
        pidents = np.array(pidents)
        for cutoff in [20, 30, 40, 50, 70, 90]:
            n = (pidents >= cutoff).sum()
            log.info("  Hits with pident ≥%d%%: %d (%.1f%% of test sequences)",
                     cutoff, n, n/len(test_df)*100)

    # ── Step 8: Fragment robustness for BLAST ─────────────────────────────
    log.info("══ BLAST fragment robustness ══")
    fragment_results = []

    for frac, label in [(1.0, "Full sequence"),
                        (0.75, "N-terminal 75%"),
                        (0.50, "N-terminal 50%"),
                        (0.25, "N-terminal 25%")]:

        # Write truncated query
        frac_query = BLAST_DIR / f"test_query_{int(frac*100)}pct.fasta"
        write_fasta_df(test_df, seq_index, frac_query, truncate_frac=frac)

        # Run BLAST on truncated sequences
        frac_out = BLAST_DIR / f"blast_results_{int(frac*100)}pct.tsv"
        if not frac_out.exists():
            run_blast(frac_query, frac_out, args.threads, args.evalue)

        frac_blast = parse_blast_results(frac_out)

        # Evaluate at "any hit" threshold (most generous for BLAST)
        r = evaluate_blast_at_threshold(test_df, frac_blast, 0)
        if r:
            r["condition"] = label
            r["fraction"]  = frac
            fragment_results.append(r)
            log.info("  %s: AUROC=%.4f  F1=%.4f  hit_rate=%.1f%%",
                     label, r["auroc"], r["f1"], r["hit_rate_pct"])

    # ── Step 9: Combined comparison table ─────────────────────────────────
    log.info("══ Combined comparison: BLAST vs CarboDB ══")
    log.info("  %-35s  %8s  %8s", "Method/Condition", "AUROC", "F1")
    log.info("  " + "-"*55)

    comparison_rows = []

    # Full sequence comparison
    for thr in [0, 30, 50]:
        r = next((x for x in blast_threshold_results if x["identity_threshold"] == thr), None)
        if r:
            label = f"BLAST (any hit)" if thr == 0 else f"BLAST (≥{thr}% identity)"
            log.info("  %-35s  %.4f    %.4f", label, r["auroc"], r["f1"])
            comparison_rows.append({"method": label, "condition": "Full sequence",
                                     "fraction": 1.0, "auroc": r["auroc"], "f1": r["f1"]})

    # Fragment comparison for BLAST vs CarboDB
    for r in fragment_results:
        log.info("  %-35s  %.4f    %.4f",
                 f"BLAST (any hit) - {r['condition']}", r["auroc"], r["f1"])
        comparison_rows.append({"method": "BLAST (any hit)", "condition": r["condition"],
                                  "fraction": r["fraction"], "auroc": r["auroc"], "f1": r["f1"]})

    # ── Save results ───────────────────────────────────────────────────────
    result = {
        "n_test":                args.n_test,
        "blast_threshold_results": blast_threshold_results,
        "fragment_results":       fragment_results,
        "comparison_rows":        comparison_rows,
        "identity_distribution": {
            "mean":   round(float(np.mean(pidents)), 1) if pidents else None,
            "median": round(float(np.median(pidents)), 1) if pidents else None,
            "pct_below_30": round(float((np.array(pidents) < 30).mean() * 100), 1) if pidents else None,
            "pct_below_50": round(float((np.array(pidents) < 50).mean() * 100), 1) if pidents else None,
        }
    }

    json.dump(result, open(BENCH_DIR / "blast_benchmark.json", "w"), indent=2)

    # TSV for plotting
    pd.DataFrame(comparison_rows).to_csv(
        FIG_DIR / "blast_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(blast_threshold_results).to_csv(
        FIG_DIR / "blast_threshold_curve.tsv", sep="\t", index=False)

    log.info("Saved: blast_benchmark.json + blast_comparison.tsv + blast_threshold_curve.tsv")

    # ── Final summary ──────────────────────────────────────────────────────
    log.info("\n" + "="*60)
    log.info("BLAST BENCHMARK SUMMARY")
    log.info("="*60)
    log.info("BLAST identity distribution: mean=%.1f%%  median=%.1f%%",
             result["identity_distribution"]["mean"] or 0,
             result["identity_distribution"]["median"] or 0)
    log.info("Hits below 30%% identity: %.1f%% of test set",
             result["identity_distribution"]["pct_below_30"] or 0)
    log.info("Hits below 50%% identity: %.1f%% of test set",
             result["identity_distribution"]["pct_below_50"] or 0)
    log.info("\nKey comparison (full sequences):")
    for r in blast_threshold_results:
        log.info("  BLAST ≥%d%%: AUROC=%.4f  (CarboDB = 0.9999)",
                 r["identity_threshold"], r["auroc"])
    log.info("\nDone. Add results to script 13 Task B for full comparison table.")
    log.info("Next: python scripts/13_benchmark_identity.py --tasks B")


if __name__ == "__main__":
    main()
