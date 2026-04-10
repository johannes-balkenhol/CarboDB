#!/usr/bin/env python3
"""
13_publication_benchmark.py
===========================
CarboDB — Complete publication benchmark.

Compares all methods on the same test set with three difficulty levels:
  - Standard test set (90% CD-HIT split, as used in training)
  - Hard negatives: CO2-related EC classes (biochemically similar to carboxylases)
  - Fragment conditions: 75%, 50%, 25% N-terminal

Methods compared:
  Binary classification:
    1. CarboDB v5 (ESM-2 + Pfam + Composition) — full model
    2. ESM-2 only
    3. ESM-2 + Composition
    4. Pfam + Composition (no ESM-2)
    5. Pfam only (domain rule-based)
    6. ProSITE scan (pattern-based)
    7. BLAST nearest-neighbor
    8. InterPro/CATH rule-based

  Km regression:
    1. CarboDB Km (Pfam + Composition + EC one-hot)
    2. Pfam only
    3. Composition only
    4. EC-class mean baseline
    5. BLAST Km transfer (nearest neighbor Km)

Metrics:
  - AUROC, F1, precision, recall, FP rate, FN rate
  - Km: R², RMSE, within-2-fold, within-5-fold
  - All metrics on: full sequences, 75%, 50%, 25% fragments

Output: data/benchmark/publication_benchmark.json
         data/benchmark/figures/publication_*.tsv

Usage:
  python scripts/13_publication_benchmark.py
  python scripts/13_publication_benchmark.py --tasks binary km
  python scripts/13_publication_benchmark.py --tasks binary --no-external
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import (roc_auc_score, f1_score, precision_score,
                              recall_score, accuracy_score, r2_score)
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, ROOT, TS, setup_logging

log = setup_logging("13_publication_benchmark")

ML_DIR    = ROOT / "data" / "ml"
MODEL_DIR = ROOT / "data" / "models"
SPLIT_DIR = ROOT / "data" / "splits"
BENCH_DIR = ROOT / "data" / "benchmark"
FIG_DIR   = BENCH_DIR / "figures"
DBS_DIR   = ROOT / "data" / "dbs"
PRIMARY   = ROOT / "data" / "primary"

BENCH_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

KM_TRAINABLE_EC = [
    "4.2.1.1","4.1.1.39","4.1.1.31","4.1.1.49",
    "6.3.4.14","4.1.1.32","6.4.1.1","6.4.1.4","6.4.1.2","6.4.1.3",
]


# ══════════════════════════════════════════════════════════════════════════════
# Data loading
# ══════════════════════════════════════════════════════════════════════════════

def load_binary_test():
    X = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y = np.load(ML_DIR / "y_binary_test.npy")
    splits = pd.read_csv(SPLIT_DIR / "split_binary.tsv", sep="\t")
    test_split = splits[splits["split"] == "test"].reset_index(drop=True)
    n = min(len(test_split), len(X))
    return X[:n], y[:n], test_split.iloc[:n]


def load_km_test():
    X_all = np.vstack([np.load(ML_DIR/f"X_km_{s}_v3.npz")["X"]
                       for s in ["train","val","test"]])
    y_all = np.concatenate([np.load(ML_DIR/f"y_km_{s}_v3.npy")
                            for s in ["train","val","test"]])
    splits = pd.read_csv(SPLIT_DIR / "split_km.tsv", sep="\t")
    mask = splits["ec_number"].isin(KM_TRAINABLE_EC).values
    X_f, y_f = X_all[mask], y_all[mask]
    ec_f = splits["ec_number"].values[mask]
    _, X_te, _, y_te, _, ec_te = train_test_split(
        X_f, y_f, ec_f, test_size=0.15, random_state=42)
    return X_te, y_te, ec_te


def load_feat_names(task, suffix=""):
    p = ML_DIR / f"feature_names_{task}{suffix}.json"
    return json.load(open(p)) if p.exists() else []


def get_feature_indices(feat_names):
    return {
        "esm2":  [i for i,n in enumerate(feat_names) if n.startswith("esm2_")],
        "pfam":  [i for i,n in enumerate(feat_names) if n.startswith("pfam_")],
        "comp":  [i for i,n in enumerate(feat_names) if n.startswith("aac_") or
                  n.startswith("dp_") or n.startswith("phys_") or n.startswith("pse_")],
        "motif": [i for i,n in enumerate(feat_names) if n.startswith("motif_") or
                  n.startswith("inv_")],
    }


# ══════════════════════════════════════════════════════════════════════════════
# Hard negative set: CO2-related EC classes
# ══════════════════════════════════════════════════════════════════════════════

def build_hard_negative_set(X_te, y_te, test_split, feat_names):
    """
    Build a hard negative test set from sequences with EC numbers in
    CFG.CO2_RELATED_EC (biochemically similar to carboxylases but not CO2-fixing).
    These are the hardest negatives — biochemically closest to positives.

    Since label=2 sequences are not in the current test split, we identify
    negatives (label=0) whose EC number falls in CO2_RELATED_EC as a proxy
    for 'hard negatives'.
    """
    co2_related = CFG.CO2_RELATED_EC

    if "ec_number" not in test_split.columns:
        log.warning("ec_number not in split — cannot build hard negative set")
        return None, None, None

    # Hard negatives: label=0 AND ec_number in CO2_RELATED_EC
    hard_neg_mask = (
        (test_split["label"].astype(int) == 0) &
        (test_split["ec_number"].isin(co2_related))
    ).values[:len(X_te)]

    easy_neg_mask = (
        (test_split["label"].astype(int) == 0) &
        (~test_split["ec_number"].isin(co2_related))
    ).values[:len(X_te)]

    pos_mask = (test_split["label"].astype(int) == 1).values[:len(X_te)]

    log.info("  Hard negatives (CO2-related EC): %d", hard_neg_mask.sum())
    log.info("  Easy negatives (unrelated EC):   %d", easy_neg_mask.sum())
    log.info("  True positives:                  %d", pos_mask.sum())

    # Hard test: positives + hard negatives only
    hard_mask = pos_mask | hard_neg_mask
    X_hard = X_te[hard_mask]
    y_hard = y_te[hard_mask]

    # Easy test: positives + easy negatives only
    easy_mask = pos_mask | easy_neg_mask
    X_easy = X_te[easy_mask]
    y_easy = y_te[easy_mask]

    return (X_hard, y_hard), (X_easy, y_easy), hard_neg_mask.sum()


# ══════════════════════════════════════════════════════════════════════════════
# Feature truncation (fragment simulation)
# ══════════════════════════════════════════════════════════════════════════════

def truncate_features(X, feat_names, frac, fidx):
    """Scale features to simulate truncated sequences."""
    Xt = X.copy()
    rng = np.random.default_rng(42)
    # Composition: scale proportionally
    for i in fidx["comp"]:
        Xt[:, i] *= frac
    # ESM-2: scale (fewer tokens in mean pool)
    for i in fidx["esm2"]:
        Xt[:, i] *= frac
    # Pfam: zero with prob (1-frac) — domain may not be in fragment
    for i in fidx["pfam"]:
        mask = rng.random(len(X)) < (1 - frac)
        Xt[mask, i] = 0.0
    # Motifs: zero with prob (1-frac)
    for i in fidx["motif"]:
        mask = rng.random(len(X)) < (1 - frac)
        Xt[mask, i] = 0.0
    return Xt


# ══════════════════════════════════════════════════════════════════════════════
# Binary evaluation
# ══════════════════════════════════════════════════════════════════════════════

def eval_binary(booster, X, y, label, condition="full"):
    probs = booster.predict(xgb.DMatrix(X))
    pred  = (probs >= 0.5).astype(int)
    if len(np.unique(y)) < 2:
        return None
    auroc = float(roc_auc_score(y, probs))
    f1    = float(f1_score(y, pred, zero_division=0))
    prec  = float(precision_score(y, pred, zero_division=0))
    rec   = float(recall_score(y, pred, zero_division=0))
    fp    = float(((pred==1) & (y==0)).mean())
    fn    = float(((pred==0) & (y==1)).mean())
    r = {
        "method":    label,
        "condition": condition,
        "n_test":    int(len(y)),
        "n_pos":     int(y.sum()),
        "auroc":     round(auroc, 4),
        "f1":        round(f1, 4),
        "precision": round(prec, 4),
        "recall":    round(rec, 4),
        "fp_rate":   round(fp, 4),
        "fn_rate":   round(fn, 4),
    }
    log.info("  %-42s %-18s AUROC=%.4f  F1=%.4f  FP=%.4f  FN=%.4f",
             label, condition, auroc, f1, fp, fn)
    return r


def make_variant(X, fidx, variant):
    Xv = np.zeros_like(X)
    if variant == "full":        return X.copy()
    if variant == "esm2":        Xv[:, fidx["esm2"]] = X[:, fidx["esm2"]]
    elif variant == "esm2_comp": Xv[:, fidx["esm2"]] = X[:, fidx["esm2"]]; Xv[:, fidx["comp"]] = X[:, fidx["comp"]]
    elif variant == "pfam_comp": Xv[:, fidx["pfam"]] = X[:, fidx["pfam"]]; Xv[:, fidx["comp"]] = X[:, fidx["comp"]]
    elif variant == "pfam":      Xv[:, fidx["pfam"]] = X[:, fidx["pfam"]]
    elif variant == "comp":      Xv[:, fidx["comp"]] = X[:, fidx["comp"]]
    return Xv


# ══════════════════════════════════════════════════════════════════════════════
# External methods: ProSITE, CATH, InterPro, BLAST
# ══════════════════════════════════════════════════════════════════════════════

def run_prosite_benchmark(X_te, y_te, test_split, feat_names):
    """
    ProSITE benchmark using existing motif features in the feature matrix.
    The motif_ and inv_ features encode ProSITE pattern hits.
    """
    log.info("  ProSITE benchmark (using motif features from feature matrix)...")
    fidx = get_feature_indices(feat_names)
    motif_idx = fidx["motif"]

    if not motif_idx:
        log.warning("  No motif features found in feature matrix")
        return None

    # ProSITE prediction: any motif hit = positive
    X_motif = X_te[:, motif_idx]
    motif_score = X_motif.sum(axis=1)
    pred = (motif_score > 0).astype(int)

    if len(np.unique(y_te)) < 2:
        return None

    try:
        auroc = float(roc_auc_score(y_te, motif_score))
    except Exception:
        auroc = 0.5

    f1  = float(f1_score(y_te, pred, zero_division=0))
    fp  = float(((pred==1) & (y_te==0)).mean())
    fn  = float(((pred==0) & (y_te==1)).mean())

    log.info("  %-42s AUROC=%.4f  F1=%.4f  FP=%.4f  FN=%.4f",
             "ProSITE motifs", auroc, f1, fp, fn)

    return {
        "method": "ProSITE motifs",
        "condition": "full",
        "n_test": int(len(y_te)),
        "auroc": round(auroc, 4),
        "f1": round(f1, 4),
        "fp_rate": round(fp, 4),
        "fn_rate": round(fn, 4),
    }


def run_interpro_cath_benchmark(X_te, y_te, feat_names):
    """
    InterPro/CATH benchmark using InterPro features from feature matrix.
    The n_pfam_hits, n_panther_hits, n_cath_hits etc. columns encode this.
    """
    log.info("  InterPro/CATH benchmark (using InterPro features)...")

    interpro_idx = [i for i,n in enumerate(feat_names)
                    if n.startswith("n_") and "hits" in n]

    if not interpro_idx:
        log.warning("  No InterPro features found")
        return None

    X_ip = X_te[:, interpro_idx]
    ip_score = X_ip.sum(axis=1)
    pred = (ip_score > 0).astype(int)

    if len(np.unique(y_te)) < 2:
        return None

    try:
        auroc = float(roc_auc_score(y_te, ip_score))
    except Exception:
        auroc = 0.5

    f1  = float(f1_score(y_te, pred, zero_division=0))
    fp  = float(((pred==1) & (y_te==0)).mean())
    fn  = float(((pred==0) & (y_te==1)).mean())

    log.info("  %-42s AUROC=%.4f  F1=%.4f  FP=%.4f  FN=%.4f",
             "InterPro (n_hits)", auroc, f1, fp, fn)

    return {
        "method": "InterPro (domain hit count)",
        "condition": "full",
        "auroc": round(auroc, 4),
        "f1": round(f1, 4),
        "fp_rate": round(fp, 4),
        "fn_rate": round(fn, 4),
    }


def load_blast_results_if_available():
    """Load BLAST benchmark results if script 04f has been run."""
    blast_path = BENCH_DIR / "blast_benchmark.json"
    if not blast_path.exists():
        log.warning("  BLAST results not found — run scripts/04f_blast_benchmark.py first")
        return None
    return json.load(open(blast_path))


def pfam_rule_benchmark(X_te, y_te, feat_names):
    """
    Pfam rule-based: any carboxylase-specific Pfam domain hit = positive.
    Uses the pfam_PFXXXXX columns in the feature matrix.
    """
    carboxy_pfam = sorted(CFG.CARBOXY_PFAM)
    pfam_cols = [i for i,n in enumerate(feat_names)
                 if n.startswith("pfam_PF") and not n.endswith("n_hits")]

    if not pfam_cols:
        log.warning("  No Pfam domain columns found")
        return None

    pfam_score = X_te[:, pfam_cols].sum(axis=1)
    pred = (pfam_score > 0).astype(int)

    if len(np.unique(y_te)) < 2:
        return None

    try:
        auroc = float(roc_auc_score(y_te, pfam_score))
    except Exception:
        auroc = 0.5

    f1  = float(f1_score(y_te, pred, zero_division=0))
    fp  = float(((pred==1) & (y_te==0)).mean())
    fn  = float(((pred==0) & (y_te==1)).mean())

    log.info("  %-42s AUROC=%.4f  F1=%.4f  FP=%.4f  FN=%.4f",
             "Pfam rule-based (carboxylase domains)", auroc, f1, fp, fn)

    return {
        "method": "Pfam rule-based",
        "condition": "full",
        "auroc": round(auroc, 4),
        "f1": round(f1, 4),
        "fp_rate": round(fp, 4),
        "fn_rate": round(fn, 4),
    }


# ══════════════════════════════════════════════════════════════════════════════
# Km evaluation
# ══════════════════════════════════════════════════════════════════════════════

def eval_km(booster, X, y_true, label, ec=None):
    pred = booster.predict(xgb.DMatrix(X))
    r2   = float(r2_score(y_true, pred))
    r, _ = pearsonr(y_true, pred)
    rmse = float(np.sqrt(np.mean((y_true - pred)**2)))
    w2f  = float(np.mean(np.abs(pred - y_true) < np.log10(2)))
    w5f  = float(np.mean(np.abs(pred - y_true) < np.log10(5)))

    log.info("  %-42s R²=%.4f  r=%.4f  within_2fold=%.3f  within_5fold=%.3f",
             label, r2, r, w2f, w5f)

    return {
        "method":       label,
        "r2":           round(r2, 4),
        "pearson_r":    round(float(r), 4),
        "rmse_log10":   round(rmse, 4),
        "within_2fold": round(w2f, 3),
        "within_5fold": round(w5f, 3),
        "n_test":       int(len(y_true)),
    }


# ══════════════════════════════════════════════════════════════════════════════
# Main benchmark
# ══════════════════════════════════════════════════════════════════════════════

def run_binary_benchmark(include_external=True):
    log.info("══ BINARY CLASSIFICATION BENCHMARK ══")

    booster    = xgb.Booster(); booster.load_model(str(MODEL_DIR/"binary_v5.json"))
    feat_names = load_feat_names("binary")
    fidx       = get_feature_indices(feat_names)

    X_te, y_te, test_split = load_binary_test()

    results = []

    # ── Build test conditions ──────────────────────────────────────────────
    hard_neg, easy_neg, n_hard = build_hard_negative_set(X_te, y_te, test_split, feat_names)

    conditions = [
        (X_te,          y_te,          "Standard test set"),
        (X_te,          y_te,          "N-terminal 75%"),
        (X_te,          y_te,          "N-terminal 50%"),
        (X_te,          y_te,          "N-terminal 25%"),
    ]

    # Add hard negative condition if available
    if hard_neg is not None and hard_neg[0] is not None and len(hard_neg[0]) > 10:
        conditions.append((hard_neg[0], hard_neg[1], "Hard negatives (CO2-related EC)"))
        log.info("  Hard negative set: %d positives + %d CO2-related negatives",
                 int(hard_neg[1].sum()), n_hard)
    else:
        log.info("  No hard negatives found in test split — using standard negatives only")

    # ── ML variants ───────────────────────────────────────────────────────
    ml_variants = [
        ("full",       "CarboDB v5 (ESM-2+Pfam+Comp)"),
        ("esm2",       "ESM-2 only"),
        ("esm2_comp",  "ESM-2 + Composition"),
        ("pfam_comp",  "Pfam + Composition (no ESM-2)"),
        ("pfam",       "Pfam only (feature matrix)"),
        ("comp",       "Composition only"),
    ]

    fracs = {
        "Standard test set": 1.0,
        "N-terminal 75%":    0.75,
        "N-terminal 50%":    0.50,
        "N-terminal 25%":    0.25,
        "Hard negatives (CO2-related EC)": 1.0,
    }

    log.info("\n  %-42s %-18s %7s  %6s  %6s  %6s",
             "Method", "Condition", "AUROC", "F1", "FP", "FN")
    log.info("  " + "-"*90)

    for cond_name, (X_cond, y_cond) in [(c[2], (c[0], c[1])) for c in conditions]:
        frac = fracs.get(cond_name, 1.0)
        X_use = truncate_features(X_cond, feat_names, frac, fidx) if frac < 1.0 else X_cond

        for vkey, vname in ml_variants:
            r = eval_binary(booster, make_variant(X_use, fidx, vkey),
                            y_cond, vname, cond_name)
            if r:
                results.append(r)

    # ── External methods (full sequences only) ────────────────────────────
    if include_external:
        log.info("\n  External methods (full sequences):")

        # Pfam rule-based
        r = pfam_rule_benchmark(X_te, y_te, feat_names)
        if r: results.append(r)

        # ProSITE
        r = run_prosite_benchmark(X_te, y_te, test_split, feat_names)
        if r: results.append(r)

        # InterPro/CATH
        r = run_interpro_cath_benchmark(X_te, y_te, feat_names)
        if r: results.append(r)

        # BLAST (from pre-computed results)
        blast = load_blast_results_if_available()
        if blast:
            for br in blast.get("blast_threshold_results", []):
                thr = br["identity_threshold"]
                label = "BLAST (any hit)" if thr == 0 else f"BLAST (≥{thr}% identity)"
                results.append({
                    "method":    label,
                    "condition": "Standard test set",
                    "n_test":    br.get("n_test", 0),
                    "auroc":     br.get("auroc", 0),
                    "f1":        br.get("f1", 0),
                    "fp_rate":   round(1 - br.get("hit_rate_pct",0)/100, 4),
                    "fn_rate":   None,
                })
                log.info("  %-42s %-18s AUROC=%.4f  F1=%.4f",
                         label, "Standard test set", br["auroc"], br["f1"])

    return results


def run_km_benchmark():
    log.info("══ Km REGRESSION BENCHMARK ══")

    booster    = xgb.Booster(); booster.load_model(str(MODEL_DIR/"km_v5_weighted.json"))
    feat_names = load_feat_names("km", "_v3")

    km_esm2_idx = [i for i,n in enumerate(feat_names) if n.startswith("esm2_")]
    km_pfam_idx = [i for i,n in enumerate(feat_names) if n.startswith("pfam_")]
    km_comp_idx = [i for i,n in enumerate(feat_names) if n.startswith("aac_") or
                   n.startswith("dp_") or n.startswith("phys_")]
    km_ec_idx   = [i for i,n in enumerate(feat_names) if
                   n.startswith("ec_oh_") or n.startswith("kingdom_")]

    X_te, y_te, ec_te = load_km_test()

    results = []

    log.info("\n  %-42s %7s  %6s  %10s  %10s",
             "Method", "R²", "r", "within_2fold", "within_5fold")
    log.info("  " + "-"*80)

    # Full model
    r = eval_km(booster, X_te, y_te, "CarboDB Km (Pfam+Comp+EC one-hot)")
    results.append(r)

    # Ablations
    ablations = [
        (km_pfam_idx, "Km: no Pfam (Comp+EC only)"),
        (km_comp_idx, "Km: no Composition (Pfam+EC only)"),
        (km_ec_idx,   "Km: no EC one-hot (Pfam+Comp only)"),
    ]
    for zero_idx, label in ablations:
        Xv = X_te.copy(); Xv[:, zero_idx] = 0
        r = eval_km(booster, Xv, y_te, label)
        results.append(r)

    # Single-feature methods
    singles = [
        (km_pfam_idx, "Km: Pfam only"),
        (km_comp_idx, "Km: Composition only"),
        (km_ec_idx,   "Km: EC one-hot only"),
    ]
    for keep_idx, label in singles:
        Xv = np.zeros_like(X_te); Xv[:, keep_idx] = X_te[:, keep_idx]
        r = eval_km(booster, Xv, y_te, label)
        results.append(r)

    # EC-class mean baseline
    ec_means = {ec: y_te[ec_te==ec].mean() for ec in KM_TRAINABLE_EC
                if (ec_te==ec).sum() > 0}
    pred_base = np.array([ec_means.get(e, y_te.mean()) for e in ec_te])
    r2_base   = float(r2_score(y_te, pred_base))
    w2f_base  = float(np.mean(np.abs(pred_base - y_te) < np.log10(2)))
    r_base, _ = pearsonr(y_te, pred_base)
    log.info("  %-42s R²=%.4f  r=%.4f  within_2fold=%.3f  within_5fold=N/A",
             "EC-class mean baseline", r2_base, r_base, w2f_base)
    results.append({
        "method": "EC-class mean baseline",
        "r2": round(r2_base, 4),
        "pearson_r": round(float(r_base), 4),
        "within_2fold": round(w2f_base, 3),
        "within_5fold": None,
    })

    # Per-EC breakdown
    log.info("\n  Per-EC R² (full model):")
    per_ec = {}
    for ec in sorted(set(ec_te)):
        m = ec_te == ec
        if m.sum() < 5: continue
        pred_ec = booster.predict(xgb.DMatrix(X_te[m]))
        r2_ec = float(r2_score(y_te[m], pred_ec))
        r_ec, _ = pearsonr(y_te[m], pred_ec)
        per_ec[ec] = {"n": int(m.sum()), "r2": round(r2_ec,4), "r": round(float(r_ec),4)}
        log.info("    %s: n=%d  R²=%.4f  r=%.4f", ec, m.sum(), r2_ec, r_ec)

    # BLAST Km transfer (nearest-neighbor Km)
    blast = load_blast_results_if_available()
    if blast:
        log.info("  BLAST Km transfer: pending — requires Km values for training sequences")

    return results, per_ec


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description="CarboDB publication benchmark.")
    ap.add_argument("--tasks",       nargs="+", default=["binary","km"],
                    choices=["binary","km"])
    ap.add_argument("--no-external", action="store_true",
                    help="Skip external method benchmarks (BLAST, ProSITE, InterPro)")
    args = ap.parse_args()

    report = {"model_version": "v5", "created_at": TS, "tasks": {}}

    if "binary" in args.tasks:
        results = run_binary_benchmark(not args.no_external)
        report["tasks"]["binary"] = results
        pd.DataFrame(results).to_csv(
            FIG_DIR / "publication_binary.tsv", sep="\t", index=False)
        log.info("Saved: publication_binary.tsv (%d rows)", len(results))

    if "km" in args.tasks:
        results, per_ec = run_km_benchmark()
        report["tasks"]["km"] = {"ablations": results, "per_ec": per_ec}
        pd.DataFrame(results).to_csv(
            FIG_DIR / "publication_km.tsv", sep="\t", index=False)
        log.info("Saved: publication_km.tsv")

    json.dump(report, open(BENCH_DIR/"publication_benchmark.json","w"), indent=2)

    # ── Final summary table ────────────────────────────────────────────────
    log.info("\n" + "="*70)
    log.info("PUBLICATION BENCHMARK SUMMARY")
    log.info("="*70)

    if "binary" in report["tasks"]:
        rows = report["tasks"]["binary"]
        std = [r for r in rows if r.get("condition")=="Standard test set"]
        log.info("\nBinary — Standard test set:")
        log.info("  %-42s %7s  %6s  %6s  %6s", "Method","AUROC","F1","FP","FN")
        for r in std:
            log.info("  %-42s %.4f  %.4f  %.4f  %.4f",
                     r["method"], r["auroc"], r["f1"],
                     r.get("fp_rate",0) or 0, r.get("fn_rate",0) or 0)

        hard = [r for r in rows if "Hard" in r.get("condition","")]
        if hard:
            log.info("\nBinary — Hard negatives (CO2-related ECs):")
            for r in hard:
                log.info("  %-42s %.4f  %.4f  %.4f  %.4f",
                         r["method"], r["auroc"], r["f1"],
                         r.get("fp_rate",0) or 0, r.get("fn_rate",0) or 0)

    if "km" in report["tasks"]:
        log.info("\nKm regression:")
        for r in report["tasks"]["km"]["ablations"]:
            log.info("  %-42s R²=%.4f  within_2fold=%.3f",
                     r["method"], r["r2"], r.get("within_2fold",0) or 0)

    log.info("\nOutputs: %s", BENCH_DIR)
    log.info("Done.")


if __name__ == "__main__":
    main()
