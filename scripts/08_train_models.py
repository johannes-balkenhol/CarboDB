#!/usr/bin/env python3
"""
08_train_models.py
==================
CarboDB — Step 08: Train XGBoost v5 models for all three tasks.

Tasks:
  1. Binary classification:  is this sequence a CO2-fixing carboxylase?
  2. EC class prediction:    which of the 27 EC classes?
  3. Km regression:          what is the log10(Km_CO2) in mM?

Input:  data/ml/  X_*.npz  y_*.npy  (from script 07)
Output: data/models/
  binary_v5.json          XGBoost binary classifier
  ec_v5.json              XGBoost EC multiclass classifier
  km_v5.json              XGBoost Km regressor
  training_report_v5.json metrics + feature importance summary

Usage:
  python scripts/08_train_models.py
  python scripts/08_train_models.py --tasks binary ec km
  python scripts/08_train_models.py --tasks km --n-estimators 2000
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import xgboost as xgb
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    f1_score,
    accuracy_score,
    mean_squared_error,
)
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("08_train_models")

ML_DIR    = PATHS.PRIMARY.parent / "ml"
MODEL_DIR = PATHS.PRIMARY.parent / "models"
MODEL_VERSION = "v5"


# ── Data loading ──────────────────────────────────────────────────────────

def load_split(task, split):
    suffix = "_fixed" if task == "ec" else ""
    X = np.load(ML_DIR / f"X_{task}_{split}{suffix}.npz")["X"]
    y = np.load(ML_DIR / f"y_{task}_{split}{suffix}.npy")
    return X, y


def load_ec_map():
    p = ML_DIR / "ec_label_map_fixed.json"
    if not p.exists():
        p = ML_DIR / "ec_label_map.json"
    with open(p) as f:
        return json.load(f)


def load_feature_names(task):
    with open(ML_DIR / f"feature_names_{task}.json") as f:
        return json.load(f)


# ── XGBoost params ────────────────────────────────────────────────────────

def binary_params(n_estimators, scale_pos_weight):
    return dict(
        n_estimators      = n_estimators,
        max_depth         = 8,
        learning_rate     = 0.05,
        subsample         = 0.8,
        colsample_bytree  = 0.8,
        min_child_weight  = 5,
        gamma             = 0.1,
        reg_alpha         = 0.1,
        reg_lambda        = 1.0,
        scale_pos_weight  = scale_pos_weight,
        objective         = "binary:logistic",
        eval_metric       = "aucpr",
        tree_method       = "hist",
        device            = "cuda",
        early_stopping_rounds = 50,
        verbosity         = 1,
        random_state      = 42,
    )


def ec_params(n_estimators, num_class):
    return dict(
        n_estimators      = n_estimators,
        max_depth         = 10,
        learning_rate     = 0.05,
        subsample         = 0.8,
        colsample_bytree  = 0.6,
        min_child_weight  = 3,
        gamma             = 0.1,
        reg_alpha         = 0.1,
        reg_lambda        = 1.0,
        objective         = "multi:softprob",
        num_class         = num_class,
        eval_metric       = "mlogloss",
        tree_method       = "hist",
        device            = "cuda",
        early_stopping_rounds = 50,
        verbosity         = 1,
        random_state      = 42,
    )


def km_params(n_estimators):
    return dict(
        n_estimators      = n_estimators,
        max_depth         = 8,
        learning_rate     = 0.02,
        subsample         = 0.8,
        colsample_bytree  = 0.8,
        min_child_weight  = 3,
        gamma             = 0.05,
        reg_alpha         = 0.1,
        reg_lambda        = 1.0,
        objective         = "reg:squarederror",
        eval_metric       = "rmse",
        tree_method       = "hist",
        device            = "cuda",
        early_stopping_rounds = 50,
        verbosity         = 1,
        random_state      = 42,
    )


# ── Evaluation ────────────────────────────────────────────────────────────

def eval_binary(model, X_test, y_test):
    probs = model.predict_proba(X_test)[:, 1]
    preds = (probs >= 0.5).astype(int)
    return {
        "auroc":     float(roc_auc_score(y_test, probs)),
        "auprc":     float(average_precision_score(y_test, probs)),
        "f1":        float(f1_score(y_test, preds)),
        "accuracy":  float(accuracy_score(y_test, preds)),
        "n_test":    int(len(y_test)),
        "n_pos":     int(y_test.sum()),
        "n_neg":     int((y_test == 0).sum()),
    }


def eval_ec(model, X_test, y_test, ec_map):
    probs = model.predict_proba(X_test)
    preds = probs.argmax(axis=1)
    top3  = np.argsort(probs, axis=1)[:, -3:]
    top3_correct = sum(y_test[i] in top3[i] for i in range(len(y_test)))

    # Per-class accuracy
    inv_map = {v: k for k, v in ec_map.items()}
    per_class = {}
    for cls_int in np.unique(y_test):
        mask = y_test == cls_int
        per_class[inv_map.get(cls_int, str(cls_int))] = {
            "n":        int(mask.sum()),
            "accuracy": float((preds[mask] == y_test[mask]).mean()),
        }

    return {
        "top1_accuracy": float(accuracy_score(y_test, preds)),
        "top3_accuracy": float(top3_correct / len(y_test)),
        "f1_macro":      float(f1_score(y_test, preds, average="macro")),
        "f1_weighted":   float(f1_score(y_test, preds, average="weighted")),
        "n_test":        int(len(y_test)),
        "n_classes":     int(len(ec_map)),
        "per_class":     per_class,
    }


def eval_km(model, X_test, y_test):
    preds = model.predict(X_test)
    mse   = mean_squared_error(y_test, preds)
    r, _  = pearsonr(y_test, preds)
    ss_res = np.sum((y_test - preds) ** 2)
    ss_tot = np.sum((y_test - y_test.mean()) ** 2)
    r2    = 1 - ss_res / ss_tot

    # Back-transform to mM for interpretability
    preds_mM = 10 ** preds
    true_mM  = 10 ** y_test
    fold_err = np.median(np.abs(preds_mM - true_mM) / true_mM)

    return {
        "r2":           float(r2),
        "rmse_log10":   float(np.sqrt(mse)),
        "pearson_r":    float(r),
        "median_fold_error_mM": float(fold_err),
        "n_test":       int(len(y_test)),
        "km_range_mM":  [float(10**y_test.min()), float(10**y_test.max())],
    }


def top_features(model, feature_names, n=20):
    scores = model.feature_importances_
    idx    = np.argsort(scores)[::-1][:n]
    return [{"feature": feature_names[i], "importance": float(scores[i])}
            for i in idx]


# ── Training ──────────────────────────────────────────────────────────────

def train_binary(n_estimators):
    log.info("══ Training BINARY classifier ══")
    t0 = time.time()

    X_train, y_train = load_split("binary", "train")
    X_val,   y_val   = load_split("binary", "val")
    X_test,  y_test  = load_split("binary", "test")

    log.info("  Train: %s  pos=%d  neg=%d",
             X_train.shape, int(y_train.sum()), int((y_train==0).sum()))
    log.info("  Val:   %s", X_val.shape)
    log.info("  Test:  %s", X_test.shape)

    # Class imbalance weight
    neg = (y_train == 0).sum()
    pos = (y_train == 1).sum()
    spw = neg / pos
    log.info("  scale_pos_weight: %.2f", spw)

    params = binary_params(n_estimators, spw)
    model  = xgb.XGBClassifier(**params)

    model.fit(
        X_train, y_train,
        eval_set=[(X_val, y_val)],
        verbose=100,
    )

    log.info("  Best iteration: %d", model.best_iteration)

    metrics = eval_binary(model, X_test, y_test)
    log.info("  AUROC:  %.4f", metrics["auroc"])
    log.info("  AUPRC:  %.4f", metrics["auprc"])
    log.info("  F1:     %.4f", metrics["f1"])
    log.info("  Acc:    %.4f", metrics["accuracy"])

    feat_names = load_feature_names("binary")
    top = top_features(model, feat_names)
    log.info("  Top features: %s", [f["feature"] for f in top[:5]])

    elapsed = time.time() - t0
    log.info("  Training time: %.1f min", elapsed / 60)

    return model, metrics, top


def train_ec(n_estimators):
    log.info("══ Training EC CLASS classifier ══")
    t0 = time.time()

    X_train, y_train = load_split("ec", "train")
    X_val,   y_val   = load_split("ec", "val")
    X_test,  y_test  = load_split("ec", "test")
    ec_map = load_ec_map()
    num_class = len(ec_map)

    log.info("  Train: %s  classes: %d", X_train.shape, num_class)
    log.info("  Class distribution (train): %s",
             dict(zip(*np.unique(y_train, return_counts=True))))

    params = ec_params(n_estimators, num_class)
    model  = xgb.XGBClassifier(**params)

    model.fit(
        X_train, y_train,
        eval_set=[(X_val, y_val)],
        verbose=100,
    )

    log.info("  Best iteration: %d", model.best_iteration)

    metrics = eval_ec(model, X_test, y_test, ec_map)
    log.info("  Top-1 accuracy: %.4f", metrics["top1_accuracy"])
    log.info("  Top-3 accuracy: %.4f", metrics["top3_accuracy"])
    log.info("  F1 macro:       %.4f", metrics["f1_macro"])

    feat_names = load_feature_names("ec")
    top = top_features(model, feat_names)
    log.info("  Top features: %s", [f["feature"] for f in top[:5]])

    elapsed = time.time() - t0
    log.info("  Training time: %.1f min", elapsed / 60)

    return model, metrics, top


def train_km(n_estimators):
    log.info("══ Training Km REGRESSOR ══")
    t0 = time.time()

    X_train, y_train = load_split("km", "train")
    X_val,   y_val   = load_split("km", "val")
    X_test,  y_test  = load_split("km", "test")

    log.info("  Train: %s  log10Km range: %.2f to %.2f",
             X_train.shape, y_train.min(), y_train.max())
    log.info("  Test:  %s", X_test.shape)

    params = km_params(n_estimators)
    model  = xgb.XGBRegressor(**params)

    model.fit(
        X_train, y_train,
        eval_set=[(X_val, y_val)],
        verbose=100,
    )

    log.info("  Best iteration: %d", model.best_iteration)

    metrics = eval_km(model, X_test, y_test)
    log.info("  R²:              %.4f", metrics["r2"])
    log.info("  RMSE (log10):    %.4f", metrics["rmse_log10"])
    log.info("  Pearson r:       %.4f", metrics["pearson_r"])
    log.info("  Median fold err: %.2f", metrics["median_fold_error_mM"])

    feat_names = load_feature_names("km")
    top = top_features(model, feat_names)
    log.info("  Top features: %s", [f["feature"] for f in top[:5]])

    elapsed = time.time() - t0
    log.info("  Training time: %.1f min", elapsed / 60)

    return model, metrics, top


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", nargs="+", default=["binary", "ec", "km"],
                    choices=["binary", "ec", "km"])
    ap.add_argument("--n-estimators", type=int, default=1000)
    args = ap.parse_args()

    MODEL_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Training CarboDB models %s", MODEL_VERSION)
    log.info("Tasks: %s  n_estimators: %d", args.tasks, args.n_estimators)
    log.info("GPU available: %s", xgb.XGBClassifier(
        device="cuda").get_params().get("device"))

    report = {
        "model_version": MODEL_VERSION,
        "created_at":    TS,
        "n_estimators":  args.n_estimators,
        "tasks":         {},
    }

    # ── Binary ──────────────────────────────────────────────────────────
    if "binary" in args.tasks:
        model, metrics, top = train_binary(args.n_estimators)
        model.get_booster().save_model(MODEL_DIR / f"binary_{MODEL_VERSION}.json")
        log.info("Saved: binary_%s.json", MODEL_VERSION)
        report["tasks"]["binary"] = {"metrics": metrics, "top_features": top}

    # ── EC class ────────────────────────────────────────────────────────
    if "ec" in args.tasks:
        model, metrics, top = train_ec(args.n_estimators)
        model.get_booster().save_model(MODEL_DIR / f"ec_{MODEL_VERSION}.json")
        log.info("Saved: ec_%s.json", MODEL_VERSION)
        report["tasks"]["ec"] = {"metrics": metrics, "top_features": top}

    # ── Km regression ───────────────────────────────────────────────────
    if "km" in args.tasks:
        model, metrics, top = train_km(args.n_estimators)
        model.get_booster().save_model(MODEL_DIR / f"km_{MODEL_VERSION}.json")
        log.info("Saved: km_%s.json", MODEL_VERSION)
        report["tasks"]["km"] = {"metrics": metrics, "top_features": top}

    # ── Save report ─────────────────────────────────────────────────────
    report_path = MODEL_DIR / f"training_report_{MODEL_VERSION}.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    log.info("Saved training report: %s", report_path)

    # ── Summary ─────────────────────────────────────────────────────────
    log.info("\n" + "=" * 60)
    log.info("TRAINING SUMMARY — CarboDB %s", MODEL_VERSION)
    log.info("=" * 60)
    for task, res in report["tasks"].items():
        m = res["metrics"]
        if task == "binary":
            log.info("Binary:  AUROC=%.4f  AUPRC=%.4f  F1=%.4f",
                     m["auroc"], m["auprc"], m["f1"])
        elif task == "ec":
            log.info("EC pred: Top1=%.4f  Top3=%.4f  F1_macro=%.4f",
                     m["top1_accuracy"], m["top3_accuracy"], m["f1_macro"])
        elif task == "km":
            log.info("Km regr: R²=%.4f  RMSE=%.4f  r=%.4f",
                     m["r2"], m["rmse_log10"], m["pearson_r"])

    log.info("Done. Next: python scripts/09_benchmark.py")


if __name__ == "__main__":
    main()
