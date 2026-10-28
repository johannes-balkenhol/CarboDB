#!/usr/bin/env python3
"""
Novelty check — run v5 binary classifier on Pfam-NEGATIVE marine sequences.

Tests whether v5 can find carboxylase candidates that don't hit any of the
18 known carboxylase Pfams. If v5 calls any of these as carboxylases with
high confidence (bin_prob >= 0.5), they are candidate novel carboxylases.

Reads:  data/gmgc/interim/pilot_1M.faa  (the 1M random marine sequences)
        data/gmgc/interim/candidates_pfam.tsv  (gives us the 5,182 Pfam-positive IDs)
Writes: data/gmgc/results/novelty_check.parquet
        data/gmgc/results/novelty_candidates.faa  (any positives, for follow-up)
"""
import json
import re
import sys
import time
from pathlib import Path
import random

import numpy as np
import pandas as pd
import xgboost as xgb

ROOT     = Path.home() / "Projects_shared/CarboDB_v3"
ML_DIR   = ROOT / "data/ml"
MOD_DIR  = ROOT / "data/models"
PILOT_FA = ROOT / "data/gmgc/interim/pilot_1M.faa"
PFAM_TSV = ROOT / "data/gmgc/interim/candidates_pfam.tsv"
OUT_DIR  = ROOT / "data/gmgc/results"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_FILE = OUT_DIR / "novelty_check.parquet"
NOVEL_FA = OUT_DIR / "novelty_candidates.faa"

N_SAMPLES = 100_000   # subsample size from Pfam-negative pool
BATCH     = 2000
BIN_THRESH = 0.5
SEED      = 42

random.seed(SEED)
np.random.seed(SEED)

# ── 1. LOAD MODELS ───────────────────────────────────────────────────
print("Loading v5 binary model...", file=sys.stderr)
bin_model = xgb.Booster()
bin_model.load_model(str(MOD_DIR / "binary_v5.json"))
feat_names = json.load(open(ML_DIR / "feature_names_binary.json"))
print(f"  Binary features: {len(feat_names)}", file=sys.stderr)

# ── 2. IDENTIFY Pfam-NEGATIVE SEQUENCES ──────────────────────────────
print("Loading list of Pfam-positive sequences...", file=sys.stderr)
pfam_positive = set(
    pd.read_csv(PFAM_TSV, sep="\t", usecols=["seq_id"])["seq_id"].unique()
)
print(f"  Pfam-positive (carboxylase): {len(pfam_positive):,}", file=sys.stderr)

# ── 3. STREAM-COLLECT Pfam-NEGATIVE IDS AND THEIR SEQUENCES ──────────
print("Streaming pilot FASTA, collecting Pfam-negative sequences...", file=sys.stderr)
neg_ids = []
neg_seqs = {}
with open(PILOT_FA) as f:
    sid = None
    buf = []
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if sid and sid not in pfam_positive:
                neg_seqs[sid] = "".join(buf)
                neg_ids.append(sid)
            sid = line[1:].split()[0]
            buf = []
        else:
            buf.append(line)
    if sid and sid not in pfam_positive:
        neg_seqs[sid] = "".join(buf)
        neg_ids.append(sid)

print(f"  Pfam-negative total: {len(neg_ids):,}", file=sys.stderr)

# ── 4. RANDOM SUBSAMPLE ───────────────────────────────────────────────
n_to_sample = min(N_SAMPLES, len(neg_ids))
sampled = random.sample(neg_ids, n_to_sample)
print(f"  Subsampled: {n_to_sample:,}", file=sys.stderr)

# ── 5. FEATURE FUNCTIONS (identical to script 05) ────────────────────
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

def build_feature_vector(seq, feat_names):
    comp = composition_features(seq)
    feat = {}
    for f in feat_names:
        if f in comp:
            feat[f] = comp[f]
        else:
            # Pfam features are all 0 (we know these are Pfam-negative)
            # ESM-2, EC one-hot, kingdom features also zero
            feat[f] = 0.0
    return np.array([feat.get(f, 0.0) for f in feat_names], dtype=np.float32)

# ── 6. PREDICT IN BATCHES ─────────────────────────────────────────────
print(f"\nRunning v5 binary classifier on {n_to_sample:,} Pfam-negative sequences...",
      file=sys.stderr)
results = []
t0 = time.time()
for i in range(0, n_to_sample, BATCH):
    batch_ids = sampled[i:i+BATCH]
    batch_seqs = [neg_seqs[s] for s in batch_ids]

    X = np.stack([build_feature_vector(s, feat_names) for s in batch_seqs])
    dmat = xgb.DMatrix(X, feature_names=feat_names)
    bin_probs = bin_model.predict(dmat)

    for j in range(len(batch_ids)):
        results.append({
            "seq_id":   batch_ids[j],
            "seq_len":  len(batch_seqs[j]),
            "bin_prob": float(bin_probs[j]),
            "is_carboxylase_called": bool(bin_probs[j] >= BIN_THRESH),
        })

    if (i // BATCH) % 10 == 0:
        elapsed = time.time() - t0
        print(f"  Batch {i//BATCH:4d}/{n_to_sample//BATCH}  "
              f"({i+len(batch_ids):,}/{n_to_sample:,})  "
              f"elapsed {elapsed:.1f}s", file=sys.stderr)

df = pd.DataFrame(results)

# ── 7. WRITE RESULTS ──────────────────────────────────────────────────
df.to_parquet(OUT_FILE, index=False)
print(f"\nWrote {OUT_FILE}", file=sys.stderr)

# ── 8. SUMMARY ────────────────────────────────────────────────────────
n_pos = df["is_carboxylase_called"].sum()
print(f"\n=== Novelty check summary ===", file=sys.stderr)
print(f"Pfam-negative sequences scored:       {len(df):,}", file=sys.stderr)
print(f"v5 binary positive (bin_prob >= 0.5): {n_pos:,} ({100*n_pos/len(df):.3f}%)", file=sys.stderr)
print(f"  ... at >= 0.9: {(df['bin_prob']>=0.9).sum():,}", file=sys.stderr)
print(f"  ... at >= 0.99: {(df['bin_prob']>=0.99).sum():,}", file=sys.stderr)

if n_pos > 0:
    # Write the candidate novel carboxylases as a FASTA for inspection
    print(f"\nWriting {n_pos} novelty candidates to {NOVEL_FA}", file=sys.stderr)
    pos_df = df[df["is_carboxylase_called"]].sort_values("bin_prob", ascending=False)
    with open(NOVEL_FA, "w") as f:
        for _, row in pos_df.iterrows():
            sid = row["seq_id"]
            f.write(f">{sid} bin_prob={row['bin_prob']:.4f}\n")
            seq = neg_seqs[sid]
            for k in range(0, len(seq), 80):
                f.write(seq[k:k+80] + "\n")
    print(f"Top 20 candidates by bin_prob:", file=sys.stderr)
    print(pos_df[["seq_id", "bin_prob", "seq_len"]].head(20).to_string(index=False),
          file=sys.stderr)
else:
    print("No Pfam-negative sequences scored as carboxylase. Pfam prefilter is safe.",
          file=sys.stderr)

elapsed = time.time() - t0
print(f"\nTotal runtime: {elapsed:.1f}s", file=sys.stderr)
