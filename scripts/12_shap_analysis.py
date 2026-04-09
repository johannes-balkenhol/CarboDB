#!/usr/bin/env python3
"""
12_shap_analysis.py
===================
CarboDB — Step 12: Deep SHAP feature importance analysis.

Answers these specific questions:
  Q1. Which features drive CO2 carboxylase prediction (binary)?
      → Global SHAP on binary model, grouped by feature type
  Q2. Which features are important per EC class?
      → Per-class SHAP on EC model — top features for each of the 26 EC classes
  Q3. Which features are shared across all predicted carboxylases?
      → Intersection of top-N binary SHAP features
  Q4. Which features drive Km prediction (EC-independent)?
      → Global SHAP on Km model (already in script 09, extended here)
  Q5. Which features drive Km per EC class?
      → Subset SHAP: run separately on sequences of each EC class
  Q6. Feature overlap across EC classes for Km?
      → Pairwise Jaccard similarity of top-K Km SHAP features between EC classes
      → Heatmap + cluster groups

Outputs: data/shap/
  shap_binary_global.json          Q1 — top features + group importance
  shap_ec_per_class.json           Q2+Q3 — per-class top features
  shap_ec_shared.json              Q3 — features shared across all EC classes
  shap_km_global.json              Q4 — global Km feature importance
  shap_km_per_ec.json              Q5 — per-EC Km feature importance
  shap_km_overlap.json             Q6 — pairwise overlap matrix + cluster groups
  figures/                         TSV tables for plotting

Usage:
  python scripts/12_shap_analysis.py                    # all analyses
  python scripts/12_shap_analysis.py --tasks q1 q2      # specific questions
  python scripts/12_shap_analysis.py --top-k 30         # top-K features per class
  python scripts/12_shap_analysis.py --sample 1000      # SHAP sample size
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, ROOT, TS, setup_logging

log = setup_logging("12_shap_analysis")

ML_DIR    = ROOT / "data" / "ml"
MODEL_DIR = ROOT / "data" / "models"
BENCH_DIR = ROOT / "data" / "benchmark"
SHAP_DIR  = ROOT / "data" / "shap"
FIG_DIR   = SHAP_DIR / "figures"
SPLIT_DIR = ROOT / "data" / "splits"

SHAP_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Trainable EC classes for Km
KM_TRAINABLE_EC = [
    "4.2.1.1", "4.1.1.39", "4.1.1.31", "4.1.1.49",
    "6.3.4.14", "4.1.1.32", "6.4.1.1", "6.4.1.4",
    "6.4.1.2", "6.4.1.3",
]

EC_NAMES = {
    "4.1.1.39":  "RuBisCO",
    "4.2.1.1":   "Carbonic anhydrase",
    "6.3.4.16":  "ACC biotin carboxylase",
    "6.3.4.14":  "Pyruvate carboxylase",
    "6.3.5.5":   "Carbamoyl-P synthase",
    "6.3.4.18":  "3-MCC",
    "4.1.1.49":  "PEPC",
    "6.3.3.3":   "Dethiobiotin synthase",
    "4.1.1.31":  "PEPCK-CO2",
    "4.1.1.112": "2-OG carboxylase",
    "4.1.1.32":  "PEPCK-GTP",
    "6.4.1.1":   "Pyruvate carboxylase (6.4.1.1)",
    "6.4.1.2":   "ACC",
    "6.4.1.3":   "Propionyl-CoA carboxylase",
    "6.4.1.4":   "3-MCC (6.4.1.4)",
}


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_booster(fname: str):
    import xgboost as xgb
    path = MODEL_DIR / fname
    if not path.exists():
        log.error("Model not found: %s", path)
        return None
    b = xgb.Booster()
    b.load_model(str(path))
    log.info("Loaded: %s", fname)
    return b


def load_feat_names(task: str, suffix: str = "") -> list:
    p = ML_DIR / f"feature_names_{task}{suffix}.json"
    if p.exists():
        return json.load(open(p))
    log.warning("Feature names not found: %s", p)
    return []


def feat_group(name: str) -> str:
    """Map feature name → human-readable group label."""
    if name.startswith("aac_"):    return "Amino acid composition"
    if name.startswith("dp_"):     return "Dipeptide composition"
    if name.startswith("pse_"):    return "Pseudo-AAC"
    if name.startswith("phys_"):   return "Physicochemical"
    if name.startswith("inv_"):    return "Catalytic core motifs"
    if name.startswith("motif_"):  return "EC motifs"
    if name.startswith("pfam_"):   return "Pfam domains"
    if name.startswith("n_") or name.startswith("panther") or \
       name.startswith("cath") or name.startswith("tigrfam"): return "InterPro"
    if name.startswith("esm2_"):   return "ESM-2 embedding"
    if name.startswith("ankh_"):   return "Ankh embedding"
    if name.startswith("ec_oh_"):  return "EC one-hot"
    if name.startswith("kingdom_"): return "Kingdom"
    if name.startswith("blast_"):  return "BLAST"
    return "Other"


def top_features(shap_vals_2d: np.ndarray,
                 feat_names: list,
                 k: int = 30,
                 label: str = "global") -> list:
    """
    Compute top-K features by mean |SHAP| from a (n_samples, n_features) array.
    Returns list of dicts.
    """
    mean_abs = np.abs(shap_vals_2d).mean(axis=0)
    total    = mean_abs.sum()
    top_idx  = mean_abs.argsort()[::-1][:k]
    return [
        {
            "rank":           int(i + 1),
            "feature":        feat_names[idx] if idx < len(feat_names) else f"f{idx}",
            "group":          feat_group(feat_names[idx] if idx < len(feat_names) else ""),
            "mean_abs_shap":  round(float(mean_abs[idx]), 6),
            "pct_importance": round(float(mean_abs[idx] / total * 100), 2),
            "label":          label,
        }
        for i, idx in enumerate(top_idx)
    ]


def group_importance(feature_list: list) -> dict:
    """Aggregate feature-level importance into group-level percentages."""
    totals = {}
    grand  = 0.0
    for f in feature_list:
        g = f["group"]
        totals[g] = totals.get(g, 0.0) + f["mean_abs_shap"]
        grand += f["mean_abs_shap"]
    if grand == 0:
        return {}
    return {k: round(v / grand * 100, 2)
            for k, v in sorted(totals.items(), key=lambda x: x[1], reverse=True)}


def sample_X(X: np.ndarray, n: int, seed: int = 42) -> tuple:
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(X), min(n, len(X)), replace=False)
    return X[idx], idx


def save_tsv(rows: list, path: Path):
    """Save a list of dicts as a TSV file."""
    if not rows:
        return
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False)
    log.info("Saved TSV: %s", path)


# ══════════════════════════════════════════════════════════════════════════════
# Q1: Which features drive binary carboxylase prediction?
# ══════════════════════════════════════════════════════════════════════════════

def analysis_q1(n_sample: int, top_k: int):
    """
    Global SHAP on binary model.
    Also computes SHAP separately for true positives (confirmed carboxylases)
    vs false positives to show what distinguishes real carboxylases.
    """
    log.info("══ Q1: Binary — global feature importance ══")

    booster    = load_booster("binary_v5.json")
    feat_names = load_feat_names("binary")
    if booster is None:
        return {}

    import xgboost as xgb
    import shap

    X_te = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te = np.load(ML_DIR / "y_binary_test.npy")

    X_s, idx_s = sample_X(X_te, n_sample)
    y_s = y_te[idx_s]

    log.info("  Computing global SHAP (n=%d)...", len(X_s))
    explainer = shap.TreeExplainer(booster)
    sv        = explainer.shap_values(X_s)           # (n, n_features)

    # Save raw SHAP values for downstream analysis
    np.savez_compressed(SHAP_DIR / "shap_binary_raw.npz",
                        shap_values=sv, X_sample=X_s,
                        y_sample=y_s, sample_idx=idx_s)

    # Global top features
    top_global = top_features(sv, feat_names, k=top_k, label="global")

    # Separate: positives only (y=1) vs negatives (y=0)
    sv_pos = sv[y_s == 1]
    sv_neg = sv[y_s == 0]
    top_pos = top_features(sv_pos, feat_names, k=top_k, label="positives")
    top_neg = top_features(sv_neg, feat_names, k=top_k, label="negatives")

    result = {
        "task":             "binary",
        "n_sample":         int(len(X_s)),
        "n_features":       int(sv.shape[1]),
        "group_importance": group_importance(top_global),
        "top_global":       top_global,
        "top_carboxylases": top_pos,
        "top_non_carboxylases": top_neg,
    }

    log.info("  Group importance: %s", result["group_importance"])
    log.info("  Top-5 global: %s", [f["feature"] for f in top_global[:5]])

    # Save tables
    save_tsv(top_global, FIG_DIR / "q1_binary_global_top_features.tsv")
    save_tsv(top_pos,    FIG_DIR / "q1_binary_positives_top_features.tsv")

    # Save JSON
    out = SHAP_DIR / "shap_binary_global.json"
    json.dump(result, open(out, "w"), indent=2)
    log.info("  Saved: %s", out)

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Q2 + Q3: Per-EC-class feature importance + shared features
# ══════════════════════════════════════════════════════════════════════════════

def analysis_q2_q3(n_sample: int, top_k: int):
    """
    Q2: For each EC class, which features does the model use to identify it?
    Q3: Which features are shared across ALL EC classes (universal carboxylase signature)?
    """
    log.info("══ Q2+Q3: EC class — per-class + shared features ══")

    booster    = load_booster("ec_v5.json")
    feat_names = load_feat_names("ec")
    if booster is None:
        return {}

    import xgboost as xgb
    import shap

    X_te   = np.load(ML_DIR / "X_ec_test_fixed.npz")["X"]
    y_te   = np.load(ML_DIR / "y_ec_test_fixed.npy")
    ec_map = json.load(open(ML_DIR / "ec_label_map_fixed.json"))
    inv_map = {v: k for k, v in ec_map.items()}

    X_s, idx_s = sample_X(X_te, n_sample)
    y_s = y_te[idx_s]

    log.info("  Computing EC SHAP (n=%d, multiclass)...", len(X_s))
    explainer = shap.TreeExplainer(booster)
    sv_all    = explainer.shap_values(X_s)
    # sv_all: list of (n_samples, n_features) arrays, one per class
    # Shape: (n_classes, n_samples, n_features)
    # sv_all: list of arrays, one per class
    # Shape can vary by shap version — normalize to (n_classes, n_samples, n_features)
    # sv_all: list of arrays, one per class
    # Shape can vary by shap version — normalize to (n_classes, n_samples, n_features)
    sv_arr = np.array(sv_all)
    if sv_arr.ndim == 2:
        # shap returned (n_samples, n_features) for single output — wrap
        sv_arr = sv_arr[np.newaxis, :]
    log.info('  sv_arr shape: %s', sv_arr.shape)
    if sv_arr.ndim == 2:
        # shap returned (n_samples, n_features) for single output — wrap
        sv_arr = sv_arr[np.newaxis, :]
    log.info('  sv_arr shape: %s', sv_arr.shape)

    np.savez_compressed(SHAP_DIR / "shap_ec_raw.npz",
                        shap_values=sv_arr, X_sample=X_s,
                        y_sample=y_s, sample_idx=idx_s)

    per_class = {}
    all_top_feature_sets = []  # for Q3: intersection

    for cls_int, ec_str in sorted(inv_map.items(), key=lambda x: x[0]):
        cls_sv = sv_arr[:, :, cls_int]  # (n_samples, n_features) — shape is (n_samples, n_features, n_classes)

        # Option A: all samples (global view of what drives this class)
        tf_global = top_features(cls_sv, feat_names, k=top_k, label=f"{ec_str}_global")

        # Option B: only samples WHERE this class is the true label
        mask = y_s == cls_int
        if mask.sum() >= 5:
            cls_sv_true = cls_sv[mask]
            tf_true = top_features(cls_sv_true, feat_names, k=top_k,
                                   label=f"{ec_str}_true_positives")
        else:
            tf_true = []

        per_class[ec_str] = {
            "ec":              ec_str,
            "ec_name":         EC_NAMES.get(ec_str, ec_str),
            "n_test_samples":  int(mask.sum()),
            "group_importance": group_importance(tf_global),
            "top_features_global":    tf_global,
            "top_features_true_pos":  tf_true,
        }

        top_names = {f["feature"] for f in tf_global[:top_k]}
        all_top_feature_sets.append(top_names)
        log.info("  %s (%s): top=%s  groups=%s",
                 ec_str, EC_NAMES.get(ec_str, ""),
                 [f["feature"] for f in tf_global[:3]],
                 per_class[ec_str]["group_importance"])

    # Q3: Shared features — intersection of top-K across all EC classes
    if all_top_feature_sets:
        shared = set.intersection(*all_top_feature_sets)
        shared_50pct = set()  # features in top-K for ≥50% of classes
        n_classes = len(all_top_feature_sets)
        feature_counts = {}
        for fs in all_top_feature_sets:
            for f in fs:
                feature_counts[f] = feature_counts.get(f, 0) + 1
        shared_50pct = {f for f, c in feature_counts.items()
                        if c / n_classes >= 0.5}
        shared_75pct = {f for f, c in feature_counts.items()
                        if c / n_classes >= 0.75}

        # Get SHAP values for shared features (mean across all classes)
        sv_mean_all_classes = sv_arr.mean(axis=2)  # (n_samples, n_features) — mean over classes
        shared_features_ranked = top_features(
            sv_mean_all_classes, feat_names, k=top_k, label="shared_all_classes")
        shared_features_ranked = [f for f in shared_features_ranked
                                   if f["feature"] in shared_50pct]

    else:
        shared = set()
        shared_50pct = set()
        shared_features_ranked = []

    log.info("  Features in top-%d for ALL EC classes: %d", top_k, len(shared))
    log.info("  Features in top-%d for ≥50%% of classes: %d", top_k, len(shared_50pct))

    result_q2 = {
        "task":            "ec_per_class",
        "n_sample":        int(len(X_s)),
        "n_classes":       len(per_class),
        "top_k":           top_k,
        "per_class":       per_class,
    }

    result_q3 = {
        "task":            "ec_shared_features",
        "top_k_used":      top_k,
        "n_classes":       len(per_class),
        "n_shared_all":    len(shared),
        "n_shared_50pct":  len(shared_50pct),
        "n_shared_75pct":  len(shared_75pct),
        "shared_all_classes":    sorted(shared),
        "shared_50pct_classes":  sorted(shared_50pct),
        "shared_75pct_classes":  sorted(shared_75pct),
        "shared_features_ranked": shared_features_ranked,
        "feature_class_counts":  {f: c for f, c in
                                   sorted(feature_counts.items(),
                                          key=lambda x: x[1], reverse=True)[:50]},
    }

    # Save TSVs
    all_rows = []
    for ec_str, data in per_class.items():
        for f in data["top_features_global"]:
            f["ec"] = ec_str
            f["ec_name"] = data["ec_name"]
            all_rows.append(f)
    save_tsv(all_rows, FIG_DIR / "q2_ec_per_class_top_features.tsv")

    shared_rows = []
    for f, c in sorted(feature_counts.items(), key=lambda x: x[1], reverse=True)[:100]:
        shared_rows.append({
            "feature": f,
            "group": feat_group(f),
            "n_classes_in_top_k": c,
            "pct_classes": round(c / len(per_class) * 100, 1),
        })
    save_tsv(shared_rows, FIG_DIR / "q3_ec_shared_features.tsv")

    json.dump(result_q2, open(SHAP_DIR / "shap_ec_per_class.json", "w"), indent=2)
    json.dump(result_q3, open(SHAP_DIR / "shap_ec_shared.json", "w"), indent=2)
    log.info("  Saved: shap_ec_per_class.json, shap_ec_shared.json")

    return result_q2, result_q3


# ══════════════════════════════════════════════════════════════════════════════
# Q4: Global Km feature importance (EC-independent)
# ══════════════════════════════════════════════════════════════════════════════

def analysis_q4(n_sample: int, top_k: int):
    """
    Global SHAP on Km model — which features predict CO2 Km regardless of EC class.
    Extends the script-09 analysis with more detail.
    """
    log.info("══ Q4: Km — global EC-independent feature importance ══")

    booster    = load_booster("km_v5_weighted.json")
    feat_names = load_feat_names("km", "_v3")
    if booster is None:
        return {}

    import xgboost as xgb
    import shap
    from sklearn.model_selection import train_test_split

    X_all = np.vstack([np.load(ML_DIR / f"X_km_{s}_v3.npz")["X"]
                       for s in ["train", "val", "test"]])
    y_all = np.concatenate([np.load(ML_DIR / f"y_km_{s}_v3.npy")
                            for s in ["train", "val", "test"]])
    km_splits = pd.read_csv(SPLIT_DIR / "split_km.tsv", sep="\t")
    mask      = km_splits["ec_number"].isin(KM_TRAINABLE_EC).values
    X_filt    = X_all[mask]
    y_filt    = y_all[mask]
    ec_filt   = km_splits["ec_number"].values[mask]

    _, X_te, _, y_te, _, ec_te = train_test_split(
        X_filt, y_filt, ec_filt, test_size=0.15, random_state=42)

    X_s, idx_s = sample_X(X_te, n_sample)
    y_s  = y_te[idx_s]
    ec_s = ec_te[idx_s]

    log.info("  Computing Km SHAP (n=%d)...", len(X_s))
    explainer = shap.TreeExplainer(booster)
    sv        = explainer.shap_values(X_s)

    np.savez_compressed(SHAP_DIR / "shap_km_raw.npz",
                        shap_values=sv, X_sample=X_s,
                        y_sample=y_s, ec_sample=ec_s, sample_idx=idx_s)

    top_global = top_features(sv, feat_names, k=top_k, label="km_global")

    # Separate: low Km (top 25%) vs high Km (bottom 25%)
    low_mask  = y_s <= np.percentile(y_s, 25)
    high_mask = y_s >= np.percentile(y_s, 75)
    top_low_km  = top_features(sv[low_mask],  feat_names, k=top_k, label="low_km")
    top_high_km = top_features(sv[high_mask], feat_names, k=top_k, label="high_km")

    result = {
        "task":             "km_global",
        "n_sample":         int(len(X_s)),
        "n_features":       int(sv.shape[1]),
        "km_range_log10":   [round(float(y_s.min()), 3), round(float(y_s.max()), 3)],
        "group_importance": group_importance(top_global),
        "top_global":       top_global,
        "top_low_km":       top_low_km,    # features that push Km prediction low
        "top_high_km":      top_high_km,   # features that push Km prediction high
    }

    log.info("  Group importance: %s", result["group_importance"])
    log.info("  Top-5 Km: %s", [f["feature"] for f in top_global[:5]])

    save_tsv(top_global,   FIG_DIR / "q4_km_global_top_features.tsv")
    save_tsv(top_low_km,   FIG_DIR / "q4_km_low_km_features.tsv")
    save_tsv(top_high_km,  FIG_DIR / "q4_km_high_km_features.tsv")

    json.dump(result, open(SHAP_DIR / "shap_km_global.json", "w"), indent=2)
    log.info("  Saved: shap_km_global.json")

    # Return raw values for Q5/Q6
    return result, sv, feat_names, X_s, y_s, ec_s


# ══════════════════════════════════════════════════════════════════════════════
# Q5: Per-EC Km feature importance
# ══════════════════════════════════════════════════════════════════════════════

def analysis_q5(sv: np.ndarray, feat_names: list,
                X_s: np.ndarray, y_s: np.ndarray,
                ec_s: np.ndarray, top_k: int):
    """
    Q5: For each EC class, which features determine Km within that EC?
    Uses the pre-computed SHAP values from Q4, subsetted per EC.
    """
    log.info("══ Q5: Km — per-EC feature importance ══")

    per_ec = {}
    for ec in sorted(set(ec_s)):
        mask = ec_s == ec
        n    = mask.sum()
        if n < 10:
            log.info("  Skipping %s (n=%d < 10)", ec, n)
            continue

        sv_ec = sv[mask]
        y_ec  = y_s[mask]
        tf    = top_features(sv_ec, feat_names, k=top_k, label=f"km_{ec}")

        # Correlation of top SHAP feature with actual Km (directionality)
        enriched_features = []
        for f in tf:
            fname = f["feature"]
            if fname in feat_names:
                fidx  = feat_names.index(fname)
                if fidx < X_s.shape[1]:
                    x_col = X_s[mask, fidx]
                    # Pearson r of feature value with Km
                    if x_col.std() > 0:
                        r = float(np.corrcoef(x_col, y_ec)[0, 1])
                    else:
                        r = 0.0
                    f["km_pearson_r"] = round(r, 3)
                    f["direction"] = "high_km" if r > 0 else "low_km"
            enriched_features.append(f)

        per_ec[ec] = {
            "ec":              ec,
            "ec_name":         EC_NAMES.get(ec, ec),
            "n_samples":       int(n),
            "km_mean_log10":   round(float(y_ec.mean()), 3),
            "km_std_log10":    round(float(y_ec.std()), 3),
            "group_importance": group_importance(tf),
            "top_features":    enriched_features,
        }
        log.info("  %s (n=%d, mean_km=10^%.2f): top=%s  groups=%s",
                 ec, n, y_ec.mean(),
                 [f["feature"] for f in tf[:3]],
                 per_ec[ec]["group_importance"])

    result = {
        "task":     "km_per_ec",
        "top_k":    top_k,
        "per_ec":   per_ec,
    }

    # TSV: one row per (EC, feature)
    rows = []
    for ec, data in per_ec.items():
        for f in data["top_features"]:
            f["ec"] = ec
            f["ec_name"] = data["ec_name"]
            rows.append(f)
    save_tsv(rows, FIG_DIR / "q5_km_per_ec_top_features.tsv")

    json.dump(result, open(SHAP_DIR / "shap_km_per_ec.json", "w"), indent=2)
    log.info("  Saved: shap_km_per_ec.json")

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Q6: Feature overlap across EC classes for Km
# ══════════════════════════════════════════════════════════════════════════════

def analysis_q6(per_ec_result: dict, top_k: int):
    """
    Q6: How much do EC classes share Km-driving features?
    Computes pairwise Jaccard similarity of top-K Km features between EC classes.
    Identifies clusters of EC classes with similar Km mechanisms.
    """
    log.info("══ Q6: Km — feature overlap across EC classes ══")

    per_ec = per_ec_result["per_ec"]
    ec_list = sorted(per_ec.keys())

    if len(ec_list) < 2:
        log.warning("  Not enough EC classes for overlap analysis")
        return {}

    # Build feature sets per EC (top_k features by mean |SHAP|)
    feature_sets = {}
    for ec in ec_list:
        feature_sets[ec] = {f["feature"] for f in per_ec[ec]["top_features"][:top_k]}

    # Pairwise Jaccard similarity
    n = len(ec_list)
    jaccard_matrix = np.zeros((n, n))
    for i, ec_a in enumerate(ec_list):
        for j, ec_b in enumerate(ec_list):
            a, b = feature_sets[ec_a], feature_sets[ec_b]
            inter = len(a & b)
            union = len(a | b)
            jaccard_matrix[i, j] = inter / union if union > 0 else 0.0

    # Pairwise shared features (which specific features are shared)
    pairwise_shared = {}
    for i, ec_a in enumerate(ec_list):
        for j, ec_b in enumerate(ec_list):
            if i >= j:
                continue
            shared = sorted(feature_sets[ec_a] & feature_sets[ec_b])
            if shared:
                key = f"{ec_a}_vs_{ec_b}"
                pairwise_shared[key] = {
                    "ec_a": ec_a,
                    "ec_b": ec_b,
                    "jaccard": round(float(jaccard_matrix[i, j]), 3),
                    "n_shared": len(shared),
                    "shared_features": shared,
                    "shared_groups": list(set(feat_group(f) for f in shared)),
                }

    # Cluster EC classes by Jaccard similarity (simple linkage)
    try:
        from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
        from scipy.spatial.distance import squareform

        dist_matrix = 1.0 - jaccard_matrix
        np.fill_diagonal(dist_matrix, 0.0)
        condensed = squareform(dist_matrix)
        Z = linkage(condensed, method="average")

        # Cut tree at distance 0.5 (Jaccard >= 0.5 = same cluster)
        labels = fcluster(Z, t=0.5, criterion="distance")
        clusters = {}
        for ec, lbl in zip(ec_list, labels):
            clusters.setdefault(int(lbl), []).append(ec)
        log.info("  Found %d Km mechanism clusters:", len(clusters))
        for cid, members in sorted(clusters.items()):
            log.info("    Cluster %d: %s", cid,
                     [f"{ec} ({EC_NAMES.get(ec,ec)})" for ec in members])
    except ImportError:
        log.warning("  scipy not available for clustering")
        clusters = {i+1: [ec] for i, ec in enumerate(ec_list)}

    # Global feature frequency (how many EC classes does each feature appear in?)
    feature_freq = {}
    for ec in ec_list:
        for f in feature_sets[ec]:
            feature_freq[f] = feature_freq.get(f, 0) + 1

    universal_km_features = sorted(
        [(f, c) for f, c in feature_freq.items()],
        key=lambda x: x[1], reverse=True)[:top_k]

    result = {
        "task":                    "km_overlap",
        "ec_classes":              ec_list,
        "top_k":                   top_k,
        "jaccard_matrix":          jaccard_matrix.tolist(),
        "clusters":                {str(k): v for k, v in clusters.items()},
        "pairwise_shared":         pairwise_shared,
        "universal_km_features":   [
            {"feature": f, "group": feat_group(f),
             "n_ec_classes": c,
             "pct_ec_classes": round(c / len(ec_list) * 100, 1)}
            for f, c in universal_km_features
        ],
    }

    # Summary stats
    jac_vals = jaccard_matrix[np.triu_indices(n, k=1)]
    result["jaccard_stats"] = {
        "mean":   round(float(jac_vals.mean()), 3),
        "median": round(float(np.median(jac_vals)), 3),
        "min":    round(float(jac_vals.min()), 3),
        "max":    round(float(jac_vals.max()), 3),
    }
    log.info("  Jaccard stats: mean=%.3f median=%.3f min=%.3f max=%.3f",
             result["jaccard_stats"]["mean"], result["jaccard_stats"]["median"],
             result["jaccard_stats"]["min"],  result["jaccard_stats"]["max"])
    log.info("  Universal Km features (in all EC classes): %s",
             [f for f, c in universal_km_features if c == len(ec_list)])

    # TSVs
    # Jaccard matrix
    jac_df = pd.DataFrame(jaccard_matrix, index=ec_list, columns=ec_list)
    jac_df.to_csv(FIG_DIR / "q6_km_jaccard_matrix.tsv", sep="\t")

    # Pairwise shared
    pw_rows = [v for v in pairwise_shared.values()]
    for r in pw_rows:
        r["shared_features"] = "|".join(r["shared_features"])
        r["shared_groups"]   = "|".join(r["shared_groups"])
    save_tsv(pw_rows, FIG_DIR / "q6_km_pairwise_shared.tsv")

    # Universal features
    save_tsv(result["universal_km_features"],
             FIG_DIR / "q6_km_universal_features.tsv")

    json.dump(result, open(SHAP_DIR / "shap_km_overlap.json", "w"), indent=2)
    log.info("  Saved: shap_km_overlap.json")

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="Deep SHAP analysis for CarboDB v5 models.")
    ap.add_argument("--tasks",   nargs="+",
                    default=["q1", "q2", "q3", "q4", "q5", "q6"],
                    choices=["q1", "q2", "q3", "q4", "q5", "q6"],
                    help="Which analyses to run (default: all)")
    ap.add_argument("--sample",  type=int, default=2000,
                    help="SHAP sample size (default: 2000)")
    ap.add_argument("--top-k",   type=int, default=30,
                    help="Top-K features per class (default: 30)")
    args = ap.parse_args()

    try:
        import shap
        log.info("SHAP version: %s", shap.__version__)
    except ImportError:
        log.error("shap not installed — run: pip install shap")
        sys.exit(1)

    tasks   = set(args.tasks)
    summary = {}

    # Q1
    if "q1" in tasks:
        summary["q1"] = analysis_q1(args.sample, args.top_k)

    # Q2 + Q3 (always together)
    if "q2" in tasks or "q3" in tasks:
        res23 = analysis_q2_q3(args.sample, args.top_k)
        if res23:
            r2, r3 = res23
            summary["q2"] = r2
            summary["q3"] = r3

    # Q4 + Q5 + Q6 (share the same SHAP computation)
    sv = feat_names_km = X_s = y_s = ec_s = None

    if any(q in tasks for q in ["q4", "q5", "q6"]):
        q4_result, sv, feat_names_km, X_s, y_s, ec_s = analysis_q4(
            args.sample, args.top_k)
        summary["q4"] = q4_result

    if "q5" in tasks and sv is not None:
        q5_result = analysis_q5(sv, feat_names_km, X_s, y_s, ec_s, args.top_k)
        summary["q5"] = q5_result

        if "q6" in tasks:
            summary["q6"] = analysis_q6(q5_result, args.top_k)

    elif "q6" in tasks:
        # Load q5 result from disk if it exists
        p = SHAP_DIR / "shap_km_per_ec.json"
        if p.exists():
            q5_result = json.load(open(p))
            summary["q6"] = analysis_q6(q5_result, args.top_k)
        else:
            log.error("q5 results not found — run with --tasks q5 q6")

    # ── Final summary ─────────────────────────────────────────────────────────
    log.info("\n" + "=" * 60)
    log.info("SHAP ANALYSIS SUMMARY — CarboDB v5")
    log.info("=" * 60)

    if "q1" in summary and summary["q1"]:
        log.info("Q1 Binary group importance: %s",
                 summary["q1"].get("group_importance", {}))

    if "q3" in summary and summary["q3"]:
        log.info("Q3 Shared features across ALL EC classes: %d",
                 summary["q3"].get("n_shared_all", 0))
        log.info("Q3 Shared features in ≥50%% of EC classes: %d",
                 summary["q3"].get("n_shared_50pct", 0))

    if "q4" in summary and summary["q4"]:
        log.info("Q4 Km group importance: %s",
                 summary["q4"].get("group_importance", {}))

    if "q6" in summary and summary["q6"]:
        js = summary["q6"].get("jaccard_stats", {})
        log.info("Q6 Km Jaccard: mean=%.3f median=%.3f",
                 js.get("mean", 0), js.get("median", 0))
        log.info("Q6 Clusters: %s",
                 {k: v for k, v in summary["q6"].get("clusters", {}).items()})

    log.info("Outputs saved to: %s", SHAP_DIR)
    log.info("TSV tables in:    %s", FIG_DIR)
    log.info("Done.")


if __name__ == "__main__":
    main()
