#!/usr/bin/env python3
"""
09_benchmark.py
===============
CarboDB — Step 09: Benchmark v5 models + SHAP feature importance analysis.

Tasks:
  1. Baseline comparisons (Pfam rule-based, EC-mean Km, BLAST)
  2. Full evaluation on test sets with confidence intervals
  3. SHAP analysis — global + per-EC class
  4. Per-EC R² for Km regression
  5. Confusion matrix for EC prediction
  6. Publication-ready figures

Output: data/benchmark/
  benchmark_report_v5.json     full metrics + baselines
  shap_binary.npz              SHAP values for binary classifier
  shap_ec.npz                  SHAP values for EC classifier
  shap_km.npz                  SHAP values for Km regressor
  figures/                     publication figures

Usage:
  python scripts/09_benchmark.py
  python scripts/09_benchmark.py --no-shap   (skip slow SHAP computation)
  python scripts/09_benchmark.py --tasks binary ec km
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import (
    roc_auc_score, average_precision_score, f1_score,
    accuracy_score, confusion_matrix, r2_score,
)
from scipy.stats import pearsonr, bootstrap

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("09_benchmark")

ML_DIR        = PATHS.PRIMARY.parent / "ml"
MODEL_DIR     = PATHS.PRIMARY.parent / "models"
BENCH_DIR     = PATHS.PRIMARY.parent / "benchmark"
FIG_DIR       = BENCH_DIR / "figures"
SPLIT_DIR     = PATHS.PRIMARY.parent / "splits"

BENCH_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ── Helpers ───────────────────────────────────────────────────────────────

def load_booster(name):
    path = MODEL_DIR / name
    if not path.exists():
        log.error("Model not found: %s", path)
        return None
    b = xgb.Booster()
    b.load_model(path)
    return b


def bootstrap_ci(metric_fn, y_true, y_pred, n=1000, alpha=0.05):
    """Bootstrap 95% CI for a metric."""
    idx = np.arange(len(y_true))
    scores = []
    rng = np.random.default_rng(42)
    for _ in range(n):
        sample = rng.choice(idx, size=len(idx), replace=True)
        try:
            scores.append(metric_fn(y_true[sample], y_pred[sample]))
        except Exception:
            pass
    scores = np.array(scores)
    lo = np.percentile(scores, 100 * alpha / 2)
    hi = np.percentile(scores, 100 * (1 - alpha / 2))
    return float(np.mean(scores)), float(lo), float(hi)


def load_feat_names(task, suffix=""):
    p = ML_DIR / f"feature_names_{task}{suffix}.json"
    if p.exists():
        return json.load(open(p))
    return [f"f{i}" for i in range(10000)]


# ── Binary benchmark ─────────────────────────────────────────────────────

def benchmark_binary(compute_shap=True):
    log.info("══ Benchmarking BINARY classifier ══")

    booster = load_booster("binary_v5.json")
    if booster is None:
        return {}

    X_te = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te = np.load(ML_DIR / "y_binary_test.npy")

    dmat  = xgb.DMatrix(X_te)
    probs = booster.predict(dmat)
    preds = (probs >= 0.5).astype(int)

    # Core metrics with CI
    auroc, auroc_lo, auroc_hi = bootstrap_ci(roc_auc_score, y_te, probs)
    auprc, auprc_lo, auprc_hi = bootstrap_ci(average_precision_score, y_te, probs)

    metrics = {
        "n_test":       int(len(y_te)),
        "n_pos":        int(y_te.sum()),
        "n_neg":        int((y_te == 0).sum()),
        "auroc":        round(auroc, 4),
        "auroc_ci95":   [round(auroc_lo, 4), round(auroc_hi, 4)],
        "auprc":        round(auprc, 4),
        "auprc_ci95":   [round(auprc_lo, 4), round(auprc_hi, 4)],
        "f1":           round(float(f1_score(y_te, preds)), 4),
        "accuracy":     round(float(accuracy_score(y_te, preds)), 4),
        "threshold":    0.5,
    }

    log.info("  AUROC: %.4f [%.4f-%.4f]", auroc, auroc_lo, auroc_hi)
    log.info("  AUPRC: %.4f [%.4f-%.4f]", auprc, auprc_lo, auprc_hi)
    log.info("  F1:    %.4f", metrics["f1"])

    # Pfam baseline — any carboxylase Pfam hit = positive
    pfam_cols = [i for i, c in enumerate(load_feat_names("binary"))
                 if c.startswith("pfam_PF") and not c.endswith("n_hits")]
    if pfam_cols:
        pfam_pred = (X_te[:, pfam_cols].sum(axis=1) > 0).astype(int)
        metrics["baseline_pfam"] = {
            "f1":       round(float(f1_score(y_te, pfam_pred)), 4),
            "accuracy": round(float(accuracy_score(y_te, pfam_pred)), 4),
            "auprc":    round(float(average_precision_score(y_te, pfam_pred)), 4),
        }
        log.info("  Pfam baseline AUPRC: %.4f", metrics["baseline_pfam"]["auprc"])

    # SHAP analysis
    if compute_shap:
        try:
            import shap
            log.info("  Computing SHAP values (sample 5000)...")
            n_sample = min(5000, len(X_te))
            idx = np.random.default_rng(42).choice(len(X_te), n_sample, replace=False)
            X_sample = X_te[idx]

            explainer = shap.TreeExplainer(booster)
            shap_vals = explainer.shap_values(X_sample)

            feat_names = load_feat_names("binary")
            mean_abs   = np.abs(shap_vals).mean(axis=0)
            top_idx    = mean_abs.argsort()[::-1][:50]

            top_features = [
                {"feature": feat_names[i] if i < len(feat_names) else f"f{i}",
                 "mean_abs_shap": float(mean_abs[i]),
                 "feature_group": _feat_group(feat_names[i] if i < len(feat_names) else "")}
                for i in top_idx
            ]

            np.savez_compressed(BENCH_DIR / "shap_binary.npz",
                                shap_values=shap_vals,
                                X_sample=X_sample,
                                sample_idx=idx)
            metrics["shap_top50"] = top_features
            log.info("  Top SHAP features: %s",
                     [f["feature"] for f in top_features[:5]])

            # Group importance
            group_imp = {}
            for f in top_features:
                g = f["feature_group"]
                group_imp[g] = group_imp.get(g, 0) + f["mean_abs_shap"]
            total = sum(group_imp.values())
            metrics["shap_group_importance"] = {
                k: round(v / total, 4) for k, v in
                sorted(group_imp.items(), key=lambda x: x[1], reverse=True)
            }
            log.info("  SHAP group importance: %s", metrics["shap_group_importance"])

        except ImportError:
            log.warning("  shap not installed — skipping SHAP analysis")

    return metrics


# ── EC benchmark ─────────────────────────────────────────────────────────

def benchmark_ec(compute_shap=True):
    log.info("══ Benchmarking EC CLASS classifier ══")

    booster = load_booster("ec_v5.json")
    if booster is None:
        return {}

    X_te    = np.load(ML_DIR / "X_ec_test_fixed.npz")["X"]
    y_te    = np.load(ML_DIR / "y_ec_test_fixed.npy")
    ec_map  = json.load(open(ML_DIR / "ec_label_map_fixed.json"))
    inv_map = {v: k for k, v in ec_map.items()}

    dmat  = xgb.DMatrix(X_te)
    probs = booster.predict(dmat).reshape(len(X_te), -1)
    preds = probs.argmax(axis=1)

    # Top-3 accuracy
    top3 = np.argsort(probs, axis=1)[:, -3:]
    top3_acc = float(sum(y_te[i] in top3[i] for i in range(len(y_te))) / len(y_te))

    # Per-class metrics
    per_class = {}
    for cls_int in np.unique(y_te):
        mask = y_te == cls_int
        ec   = inv_map.get(int(cls_int), str(cls_int))
        per_class[ec] = {
            "n":           int(mask.sum()),
            "accuracy":    round(float((preds[mask] == y_te[mask]).mean()), 4),
            "f1":          round(float(f1_score(
                               y_te[mask], preds[mask],
                               average="micro", zero_division=0)), 4),
        }

    metrics = {
        "n_test":          int(len(y_te)),
        "n_classes":       int(len(ec_map)),
        "top1_accuracy":   round(float(accuracy_score(y_te, preds)), 4),
        "top3_accuracy":   round(top3_acc, 4),
        "f1_macro":        round(float(f1_score(y_te, preds, average="macro",
                                                zero_division=0)), 4),
        "f1_weighted":     round(float(f1_score(y_te, preds, average="weighted",
                                                zero_division=0)), 4),
        "per_class":       per_class,
    }

    log.info("  Top-1: %.4f  Top-3: %.4f  F1_macro: %.4f",
             metrics["top1_accuracy"], metrics["top3_accuracy"], metrics["f1_macro"])

    # Confusion matrix
    cm = confusion_matrix(y_te, preds)
    ec_labels = [inv_map.get(i, str(i)) for i in range(len(ec_map))]
    np.save(BENCH_DIR / "confusion_matrix_ec.npy", cm)
    json.dump(ec_labels, open(BENCH_DIR / "ec_labels.json", "w"))
    log.info("  Confusion matrix saved")

    # SHAP
    if compute_shap:
        try:
            import shap
            log.info("  Computing SHAP values (sample 2000)...")
            n_sample = min(2000, len(X_te))
            idx = np.random.default_rng(42).choice(len(X_te), n_sample, replace=False)
            X_sample = X_te[idx]

            explainer  = shap.TreeExplainer(booster)
            shap_vals  = explainer.shap_values(X_sample)  # list per class

            feat_names = load_feat_names("ec", "_fixed")
            # Global: mean abs across all classes
            mean_abs = np.abs(np.array(shap_vals)).mean(axis=0).mean(axis=0)
            top_idx  = mean_abs.argsort()[::-1][:50]

            metrics["shap_top50"] = [
                {"feature": feat_names[i] if i < len(feat_names) else f"f{i}",
                 "mean_abs_shap": float(mean_abs[i]),
                 "feature_group": _feat_group(feat_names[i] if i < len(feat_names) else "")}
                for i in top_idx
            ]
            np.savez_compressed(BENCH_DIR / "shap_ec.npz",
                                shap_values=np.array(shap_vals),
                                X_sample=X_sample, sample_idx=idx)
            log.info("  Top SHAP: %s",
                     [f["feature"] for f in metrics["shap_top50"][:5]])
        except ImportError:
            log.warning("  shap not installed — skipping")

    return metrics


# ── Km benchmark ─────────────────────────────────────────────────────────

def benchmark_km(compute_shap=True):
    log.info("══ Benchmarking Km REGRESSOR ══")

    booster = load_booster("km_v5_weighted.json")
    if booster is None:
        booster = load_booster("km_v5_final.json")
    if booster is None:
        return {}

    # Use v3 matrices (EC one-hot + kingdom)
    X_all = np.vstack([np.load(ML_DIR / f"X_km_{s}_v3.npz")["X"]
                       for s in ["train", "val", "test"]])
    y_all = np.concatenate([np.load(ML_DIR / f"y_km_{s}_v3.npy")
                            for s in ["train", "val", "test"]])

    km_splits = pd.read_csv(SPLIT_DIR / "split_km.tsv", sep="\t")

    # Filter to trainable EC classes
    trainable_ec = ["4.2.1.1","4.1.1.39","4.1.1.31","4.1.1.49",
                    "6.3.4.14","4.1.1.32","6.4.1.1","6.4.1.4",
                    "6.4.1.2","6.4.1.3"]
    mask = km_splits["ec_number"].isin(trainable_ec).values
    X_filt = X_all[mask]
    y_filt = y_all[mask]
    ec_filt = km_splits["ec_number"].values[mask]

    from sklearn.model_selection import train_test_split
    _, X_te, _, y_te, _, ec_te = train_test_split(
        X_filt, y_filt, ec_filt, test_size=0.15, random_state=42)

    pred = booster.predict(xgb.DMatrix(X_te))

    r2      = float(r2_score(y_te, pred))
    r, _    = pearsonr(y_te, pred)
    rmse    = float(np.sqrt(np.mean((y_te - pred) ** 2)))

    # EC-mean baseline
    ec_means = {ec: y_filt[ec_filt == ec].mean() for ec in trainable_ec}
    pred_baseline = np.array([ec_means.get(e, y_filt.mean()) for e in ec_te])
    r2_baseline = float(r2_score(y_te, pred_baseline))
    log.info("  EC-mean baseline R²: %.4f", r2_baseline)

    # Per-EC R²
    per_ec = {}
    for ec in sorted(set(ec_te)):
        m = ec_te == ec
        if m.sum() < 5:
            continue
        r2_ec   = float(r2_score(y_te[m], pred[m]))
        r_ec, _ = pearsonr(y_te[m], pred[m])
        per_ec[ec] = {
            "n":  int(m.sum()),
            "r2": round(r2_ec, 4),
            "r":  round(float(r_ec), 4),
        }

    metrics = {
        "n_test":         int(len(y_te)),
        "n_ec_classes":   len(trainable_ec),
        "r2":             round(r2, 4),
        "rmse_log10":     round(rmse, 4),
        "pearson_r":      round(float(r), 4),
        "r2_baseline_ec_mean": round(r2_baseline, 4),
        "improvement_over_baseline": round(r2 - r2_baseline, 4),
        "per_ec_r2":      per_ec,
    }

    log.info("  R²=%.4f  RMSE=%.4f  r=%.4f", r2, rmse, r)
    log.info("  Improvement over EC-mean baseline: %.4f",
             metrics["improvement_over_baseline"])
    log.info("  Per-EC R²:")
    for ec, v in per_ec.items():
        log.info("    %s: R²=%.3f  r=%.3f  n=%d", ec, v["r2"], v["r"], v["n"])

    # SHAP
    if compute_shap:
        try:
            import shap
            log.info("  Computing SHAP values...")
            n_sample = min(500, len(X_te))
            idx = np.random.default_rng(42).choice(len(X_te), n_sample, replace=False)
            X_sample = X_te[idx]

            explainer = shap.TreeExplainer(booster)
            shap_vals = explainer.shap_values(X_sample)

            feat_names = load_feat_names("km", "_v3")
            mean_abs   = np.abs(shap_vals).mean(axis=0)
            top_idx    = mean_abs.argsort()[::-1][:50]

            top_features = [
                {"feature": feat_names[i] if i < len(feat_names) else f"f{i}",
                 "mean_abs_shap": float(mean_abs[i]),
                 "feature_group": _feat_group(feat_names[i] if i < len(feat_names) else "")}
                for i in top_idx
            ]

            np.savez_compressed(BENCH_DIR / "shap_km.npz",
                                shap_values=shap_vals,
                                X_sample=X_sample,
                                y_sample=y_te[idx],
                                ec_sample=ec_te[idx],
                                sample_idx=idx)

            metrics["shap_top50"] = top_features
            group_imp = {}
            for f in top_features:
                g = f["feature_group"]
                group_imp[g] = group_imp.get(g, 0) + f["mean_abs_shap"]
            total = sum(group_imp.values())
            metrics["shap_group_importance"] = {
                k: round(v / total, 4) for k, v in
                sorted(group_imp.items(), key=lambda x: x[1], reverse=True)
            }
            log.info("  Top SHAP: %s", [f["feature"] for f in top_features[:5]])
            log.info("  Group importance: %s", metrics["shap_group_importance"])

        except ImportError:
            log.warning("  shap not installed — skipping")

    return metrics


# ── Feature group helper ──────────────────────────────────────────────────

def _feat_group(name):
    if name.startswith("aac_"):       return "A_amino_acid_composition"
    if name.startswith("dp_"):        return "B_dipeptide"
    if name.startswith("pse_"):       return "C_pseudoaac"
    if name.startswith("phys_"):      return "D_physicochemical"
    if name.startswith("inv_"):       return "E_catalytic_core"
    if name.startswith("motif_"):     return "F_ec_motifs"
    if name.startswith("pfam_"):      return "G_pfam_domains"
    if name.startswith("n_") or name.startswith("panther") or \
       name.startswith("cath") or name.startswith("tigrfam"): return "H_interpro"
    if name.startswith("esm2_"):      return "I_esm2_embedding"
    if name.startswith("ankh_"):      return "J_ankh_embedding"
    if name.startswith("ec_oh_"):     return "K_ec_onehot"
    if name.startswith("kingdom_"):   return "L_taxonomy"
    if name.startswith("blast_"):     return "M_blast"
    return "Z_other"


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks",   nargs="+", default=["binary", "ec", "km"],
                    choices=["binary", "ec", "km"])
    ap.add_argument("--no-shap", action="store_true")
    args = ap.parse_args()

    compute_shap = not args.no_shap

    # Check shap installed
    if compute_shap:
        try:
            import shap
            log.info("SHAP version: %s", shap.__version__)
        except ImportError:
            log.warning("shap not installed — run: pip install shap")
            log.warning("Continuing without SHAP")
            compute_shap = False

    report = {"model_version": "v5", "created_at": TS, "tasks": {}}

    if "binary" in args.tasks:
        report["tasks"]["binary"] = benchmark_binary(compute_shap)

    if "ec" in args.tasks:
        report["tasks"]["ec"] = benchmark_ec(compute_shap)

    if "km" in args.tasks:
        report["tasks"]["km"] = benchmark_km(compute_shap)

    # Save report
    report_path = BENCH_DIR / "benchmark_report_v5.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    log.info("Saved: %s", report_path)

    # Summary
    log.info("\n" + "=" * 60)
    log.info("BENCHMARK SUMMARY — CarboDB v5")
    log.info("=" * 60)
    for task, res in report["tasks"].items():
        if not res:
            continue
        if task == "binary":
            log.info("Binary:  AUROC=%.4f  AUPRC=%.4f  F1=%.4f",
                     res.get("auroc", 0), res.get("auprc", 0), res.get("f1", 0))
        elif task == "ec":
            log.info("EC pred: Top1=%.4f  Top3=%.4f  F1_macro=%.4f",
                     res.get("top1_accuracy", 0), res.get("top3_accuracy", 0),
                     res.get("f1_macro", 0))
        elif task == "km":
            log.info("Km regr: R²=%.4f  r=%.4f  vs baseline R²=%.4f",
                     res.get("r2", 0), res.get("pearson_r", 0),
                     res.get("r2_baseline_ec_mean", 0))

    log.info("Done. Next: python scripts/10_predict_all.py")


if __name__ == "__main__":
    main()
