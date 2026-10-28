#!/usr/bin/env python3
"""
Random Forest predicting Km from per-Pfam BITSCORES (not just EC class).

Tests the hypothesis: Pfam bitscores are continuous predictors of
EC-specific Km, with stronger Pfam matches → more canonical kinetics.

Replaces the previous RF (which used ec_pred one-hot + n_carboxy_pfam)
with the actual feature representation v5 uses internally: 18 per-Pfam
bitscore features.

Reads:  data/gmgc/results/predictions_joined.parquet
        data/gmgc/interim/candidates_pfam.tsv
Writes: data/gmgc/figures/12_bitscore_rf_importance.png
        data/gmgc/results/bitscore_rf_importance.tsv
"""
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score

ROOT  = Path.home() / "Projects_shared/CarboDB_v3"
PRED  = ROOT / "data/gmgc/results/predictions_joined.parquet"
PFAM  = ROOT / "data/gmgc/interim/candidates_pfam.tsv"
OUT_T = ROOT / "data/gmgc/results/bitscore_rf_importance.tsv"
OUT_F = ROOT / "data/gmgc/figures/12_bitscore_rf_importance.png"

CARBOXY_PFAMS = [
    "PF00016","PF00101","PF00194","PF00289","PF00311","PF00364",
    "PF00384","PF00682","PF00821","PF01039","PF01309","PF02785",
    "PF02786","PF02787","PF02788","PF03119","PF03590","PF03599",
]
PFAM_NAMES = {
    "PF00016":"RuBisCO_large","PF00101":"RuBisCO_small","PF00194":"Carb_anhydrase",
    "PF00289":"Biotin_carb_N","PF00311":"PEPcase","PF00364":"Biotin_lipoyl",
    "PF00384":"Molybdopterin","PF00682":"HMGL-like","PF00821":"PEPCK_GTP",
    "PF01039":"Carboxyl_trans","PF01309":"EAV_GS","PF02785":"Biotin_carb_C",
    "PF02786":"CPSase_L_D2","PF02787":"CPSase_L_D3","PF02788":"RuBisCO_large_N",
    "PF03119":"DNA_ligase_ZBD","PF03590":"AsnA","PF03599":"CdhD",
}

# ── 1. LOAD PREDICTIONS + PFAM BITSCORES ───────────────────────────
print("Loading predictions...", file=sys.stderr)
df = pd.read_parquet(PRED)
df = df[df['is_carboxylase']].copy()
print(f"  Carboxylase predictions: {len(df):,}", file=sys.stderr)

print("Loading per-sequence Pfam bitscores...", file=sys.stderr)
pfam_bs = defaultdict(dict)
with open(PFAM) as f:
    next(f)  # header
    for line in f:
        seq_id, acc, ev, bs = line.strip().split('\t')
        bs = float(bs)
        # Keep best bitscore per (seq, Pfam) pair
        if acc not in pfam_bs[seq_id] or bs > pfam_bs[seq_id][acc]:
            pfam_bs[seq_id][acc] = bs

print(f"  Sequences with Pfam hits: {len(pfam_bs):,}", file=sys.stderr)

# ── 2. BUILD FEATURE MATRIX ───────────────────────────────────────
print("\nBuilding bitscore matrix...", file=sys.stderr)
features = []
labels = []
seq_ids = []
for _, row in df.iterrows():
    sid = row['seq_id']
    hits = pfam_bs.get(sid, {})
    bs_vec = [hits.get(pf, 0.0) for pf in CARBOXY_PFAMS]
    # Also include seq_len + n_carboxy_pfam for orthogonal comparison
    bs_vec.append(row['seq_len'])
    bs_vec.append(row['n_carboxy_pfam'])
    bs_vec.append(row['latitude'] if pd.notna(row['latitude']) else 0)
    bs_vec.append(row['longitude'] if pd.notna(row['longitude']) else 0)
    features.append(bs_vec)
    labels.append(np.log10(row['km_pred_mM']))
    seq_ids.append(sid)

X = np.array(features)
y = np.array(labels)
feat_names = [f"bs_{pf}_{PFAM_NAMES[pf]}" for pf in CARBOXY_PFAMS] + \
             ["seq_len", "n_carboxy_pfam", "latitude", "longitude"]
print(f"  Feature matrix: {X.shape}", file=sys.stderr)

# ── 3. CROSS-VALIDATED RF ─────────────────────────────────────────
rf = RandomForestRegressor(n_estimators=400, max_depth=12, min_samples_leaf=10,
                            random_state=42, n_jobs=-1)
print("\nRunning 5-fold CV...", file=sys.stderr)
cv_r2 = cross_val_score(rf, X, y, cv=5, scoring='r2', n_jobs=-1)
print(f"  Mean R²: {cv_r2.mean():.4f} +/- {cv_r2.std():.4f}", file=sys.stderr)

rf.fit(X, y)
print(f"  Training R²: {rf.score(X, y):.4f}", file=sys.stderr)

# ── 4. FEATURE IMPORTANCE ─────────────────────────────────────────
imp = pd.Series(rf.feature_importances_, index=feat_names).sort_values(ascending=False)
imp.to_csv(OUT_T, sep='\t', header=['importance'])
print(f"\nAll feature importances (sorted):", file=sys.stderr)
for f, v in imp.items():
    print(f"  {f:45s} {v:.4f}", file=sys.stderr)

# Aggregate by category
bitscore_total = imp[imp.index.str.startswith('bs_')].sum()
print(f"\nAggregate:", file=sys.stderr)
print(f"  Pfam bitscores (18 total):   {bitscore_total:.4f}", file=sys.stderr)
print(f"  seq_len:                     {imp.get('seq_len', 0):.4f}", file=sys.stderr)
print(f"  n_carboxy_pfam:              {imp.get('n_carboxy_pfam', 0):.4f}", file=sys.stderr)
print(f"  latitude:                    {imp.get('latitude', 0):.4f}", file=sys.stderr)
print(f"  longitude:                   {imp.get('longitude', 0):.4f}", file=sys.stderr)

# ── 5. PLOT ────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 8))
top = imp.head(22)[::-1]
colors = ['steelblue' if f.startswith('bs_') else 'salmon' for f in top.index]
ax.barh(range(len(top)), top.values, color=colors)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top.index, fontsize=9)
ax.set_xlabel("Random Forest feature importance")
ax.set_title(f"Km predictors — Pfam bitscores + metadata (CV R² = {cv_r2.mean():.3f})")

from matplotlib.patches import Patch
ax.legend(handles=[
    Patch(facecolor='steelblue', label='Per-Pfam bitscore'),
    Patch(facecolor='salmon', label='Sequence / sample metadata'),
], loc='lower right')
plt.tight_layout()
plt.savefig(OUT_F, dpi=120)
plt.close()
print(f"\nSaved {OUT_F}", file=sys.stderr)
