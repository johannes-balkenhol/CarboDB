#!/usr/bin/env python3
"""
Run the v5 cascade (binary -> EC -> Km) on the Pfam-filtered candidate sequences.

Adapted from scripts/22_metagenome_scan.py with the TSV-writing bug fixed
(uses parquet output instead of TSV with broken column separation).

Reads:  data/gmgc/interim/candidates.faa
        data/gmgc/interim/candidates_pfam.tsv
        data/models/binary_v5.json
        data/models/ec_v5.json
        data/models/km_v5_weighted.json
        data/ml/feature_names_binary.json
        data/ml/feature_names_km.json
        data/ml/ec_label_map.json
Writes: data/gmgc/results/predictions_pilot.parquet
"""
import json
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

ROOT     = Path.home() / "Projects_shared/CarboDB_v3"
ML_DIR   = ROOT / "data/ml"
MOD_DIR  = ROOT / "data/models"
FASTA    = ROOT / "data/gmgc/interim/candidates.faa"
PFAM_TSV = ROOT / "data/gmgc/interim/candidates_pfam.tsv"
OUT_DIR  = ROOT / "data/gmgc/results"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FILE = OUT_DIR / "predictions_pilot.parquet"

BATCH      = 2000
BIN_THRESH = 0.5

# ── 1. LOAD MODELS ────────────────────────────────────────────────────
print("Loading v5 models...", file=sys.stderr)
bin_model = xgb.Booster(); bin_model.load_model(str(MOD_DIR / "binary_v5.json"))
ec_model  = xgb.Booster(); ec_model.load_model(str(MOD_DIR / "ec_v5.json"))
km_model  = xgb.Booster(); km_model.load_model(str(MOD_DIR / "km_v5_weighted.json"))

feat_names_bin = json.load(open(ML_DIR / "feature_names_binary.json"))
feat_names_km  = json.load(open(ML_DIR / "feature_names_km.json"))

# The Tara script used ec_label_map_fixed.json; fall back to ec_label_map.json
ec_map_path = ML_DIR / "ec_label_map_fixed.json"
if not ec_map_path.exists():
    ec_map_path = ML_DIR / "ec_label_map.json"
ec_label_map = json.load(open(ec_map_path))
ec_classes = [ec for ec, _ in sorted(ec_label_map.items(), key=lambda x: x[1])]

print(f"  Binary features: {len(feat_names_bin)}", file=sys.stderr)
print(f"  Km features:     {len(feat_names_km)}", file=sys.stderr)
print(f"  EC classes:      {len(ec_classes)}", file=sys.stderr)

# ── 2. FEATURE FUNCTIONS ─────────────────────────────────────────────
AAs = list("ACDEFGHIKLMNPQRSTVWY")

MW_MAP = {"A":89,"C":121,"D":133,"E":147,"F":165,"G":75,"H":155,"I":131,
          "K":146,"L":131,"M":149,"N":132,"P":115,"Q":146,"R":174,"S":105,
          "T":119,"V":117,"W":204,"Y":181}
CHARGE_MAP = {"D":-1,"E":-1,"K":1,"R":1,"H":0.1}
GRAVY_MAP = {"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,
             "I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,
             "R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}

def composition_features(seq):
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())
    n = len(seq) or 1
    aac = {f"aac_{a}": seq.count(a) / n for a in AAs}
    dp = {}
    for i in range(len(seq) - 1):
        k = f"dp_{seq[i]}{seq[i+1]}"
        dp[k] = dp.get(k, 0) + 1
    dp = {k: v / max(n - 1, 1) for k, v in dp.items()}
    phys = {
        "phys_length": n,
        "phys_mw":     sum(MW_MAP.get(a, 110) for a in seq) / n,
        "phys_charge": sum(CHARGE_MAP.get(a, 0) for a in seq),
        "phys_gravy":  sum(GRAVY_MAP.get(a, 0) for a in seq) / n,
    }
    return {**aac, **dp, **phys}


def build_pfam_features(pfam_hits_per_seq, feat_names):
    """pfam_hits_per_seq[seq_id] -> {pfam_acc: best_bitscore}"""
    pfam_cols = [f for f in feat_names if f.startswith("pfam_")]
    pfam_domains = [f.replace("pfam_", "") for f in pfam_cols
                    if f.replace("pfam_", "") != "n_hits"]
    out = {}
    for sid, hits in pfam_hits_per_seq.items():
        vec = {f"pfam_{d}": hits.get(d, 0.0) for d in pfam_domains}
        vec["pfam_n_hits"] = len(hits)
        out[sid] = vec
    return out


def build_feature_vector(seq, pfam_feats, feat_names):
    comp = composition_features(seq)
    feat = {}
    for f in feat_names:
        if f in comp:
            feat[f] = comp[f]
        elif f in pfam_feats:
            feat[f] = pfam_feats[f]
        else:
            feat[f] = 0.0   # ESM-2 / EC one-hot / kingdom features zero-filled
    return np.array([feat.get(f, 0.0) for f in feat_names], dtype=np.float32)


# ── 3. LOAD PFAM HITS PER SEQUENCE ───────────────────────────────────
print(f"Loading Pfam hits from {PFAM_TSV}", file=sys.stderr)
pfam_hits = defaultdict(dict)
with open(PFAM_TSV) as f:
    next(f)  # skip header
    for line in f:
        sid, acc, ev, bs = line.strip().split("\t")
        # Take best bitscore per (seq, Pfam) combination
        if acc not in pfam_hits[sid] or float(bs) > pfam_hits[sid][acc]:
            pfam_hits[sid][acc] = float(bs)
print(f"  {len(pfam_hits):,} sequences with Pfam hits", file=sys.stderr)

# ── 4. PARSE FASTA ────────────────────────────────────────────────────
print(f"Loading sequences from {FASTA}", file=sys.stderr)
seqs = {}
with open(FASTA) as f:
    sid = None
    buf = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if sid:
                seqs[sid] = "".join(buf)
            sid = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if sid:
        seqs[sid] = "".join(buf)
print(f"  {len(seqs):,} candidate sequences loaded", file=sys.stderr)

# ── 5. PREDICT IN BATCHES ─────────────────────────────────────────────
pfam_feats_per_seq = build_pfam_features(pfam_hits, feat_names_bin)

seq_ids = list(seqs.keys())
n_total = len(seq_ids)
print(f"\nRunning v5 cascade on {n_total:,} candidates in batches of {BATCH}...", file=sys.stderr)

results = []
t0 = time.time()
for i in range(0, n_total, BATCH):
    batch_ids  = seq_ids[i:i+BATCH]
    batch_seqs = [seqs[s] for s in batch_ids]
    batch_pfam = [pfam_feats_per_seq.get(s, {}) for s in batch_ids]

    # Binary prediction
    X_bin = np.stack([
        build_feature_vector(s, p, feat_names_bin)
        for s, p in zip(batch_seqs, batch_pfam)
    ])
    dbin = xgb.DMatrix(X_bin, feature_names=feat_names_bin)
    bin_probs = bin_model.predict(dbin)

    # EC + Km on all (regardless of binary threshold — we keep all rows in
    # the parquet output so downstream filtering is flexible)
    X_km = np.stack([
        build_feature_vector(s, p, feat_names_km)
        for s, p in zip(batch_seqs, batch_pfam)
    ])
    dec = xgb.DMatrix(X_km, feature_names=feat_names_km)
    ec_probs_mat = ec_model.predict(dec)
    km_preds_log = km_model.predict(dec)

    for j in range(len(batch_ids)):
        ec_argmax = int(np.argmax(ec_probs_mat[j]))
        results.append({
            "seq_id":         batch_ids[j],
            "seq_len":        len(batch_seqs[j]),
            "bin_prob":       float(bin_probs[j]),
            "is_carboxylase": bool(bin_probs[j] >= BIN_THRESH),
            "ec_pred":        ec_classes[ec_argmax],
            "ec_conf":        float(ec_probs_mat[j, ec_argmax]),
            "km_pred_log10":  float(km_preds_log[j]),
            "km_pred_mM":     float(10 ** km_preds_log[j]),
            "n_carboxy_pfam": len(pfam_hits.get(batch_ids[j], {})),
        })

    if (i // BATCH) % 10 == 0:
        elapsed = time.time() - t0
        print(f"  Batch {i // BATCH:4d} / {n_total // BATCH}  "
              f"({i+len(batch_ids):,}/{n_total:,})  "
              f"elapsed {elapsed:.1f}s", file=sys.stderr)

# ── 6. WRITE PARQUET ─────────────────────────────────────────────────
df = pd.DataFrame(results)
print(f"\nWriting {len(df):,} predictions to {OUT_FILE}", file=sys.stderr)
df.to_parquet(OUT_FILE, index=False)

print(f"\n=== Summary ===", file=sys.stderr)
print(f"Total candidates scored:   {len(df):,}", file=sys.stderr)
n_carboxy = df["is_carboxylase"].sum()
print(f"Predicted carboxylases:    {n_carboxy:,} ({100*n_carboxy/len(df):.1f}%)", file=sys.stderr)
print(f"\nTop EC classes (among predicted carboxylases):", file=sys.stderr)
print(df[df["is_carboxylase"]]["ec_pred"].value_counts().head(10), file=sys.stderr)
print(f"\nKm distribution (mM, predicted carboxylases):", file=sys.stderr)
print(df[df["is_carboxylase"]]["km_pred_mM"].describe(), file=sys.stderr)

elapsed = time.time() - t0
print(f"\nTotal runtime: {elapsed:.1f}s", file=sys.stderr)
