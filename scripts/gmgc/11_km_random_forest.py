#!/usr/bin/env python3
"""
Random Forest regression to identify what drives predicted Km.

Trains a RandomForestRegressor on the pilot predictions using features:
  - ec_pred (one-hot)
  - tax_name (top 30 genera one-hot)
  - latitude, longitude
  - GCfrac, seq_len, n_carboxy_pfam, is_complete

Reports feature importance to identify what actually predicts Km within
our dataset. If EC and genus dominate, no other metadata adds signal.
If latitude/GCfrac add signal, that's a finding.

Reads:  data/gmgc/results/predictions_joined.parquet
Writes: data/gmgc/figures/11_rf_feature_importance.png
        data/gmgc/results/rf_feature_importance.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import OneHotEncoder

ROOT  = Path.home() / "Projects_shared/CarboDB_v3"
PRED  = ROOT / "data/gmgc/results/predictions_joined.parquet"
OUT_T = ROOT / "data/gmgc/results/rf_feature_importance.tsv"
OUT_F = ROOT / "data/gmgc/figures/11_rf_feature_importance.png"

df = pd.read_parquet(PRED)
df = df[df['is_carboxylase']].copy()
print(f"Carboxylase predictions: {len(df):,}", file=sys.stderr)

# Target: log10(Km)
y = np.log10(df['km_pred_mM'].values)

# Build features
# 1. EC one-hot
ec_dummies = pd.get_dummies(df['ec_pred'], prefix='ec')

# 2. Top 30 genera one-hot (others -> "Other")
top_genera = df['tax_name'].value_counts().head(30).index.tolist()
genus_col = df['tax_name'].where(df['tax_name'].isin(top_genera), 'Other')
gen_dummies = pd.get_dummies(genus_col, prefix='gen')

# 3. Numeric features (impute missing with median)
num_features = ['latitude', 'longitude', 'GCfrac', 'seq_len', 'n_carboxy_pfam']
nums = df[num_features].copy()
for col in num_features:
    nums[col] = nums[col].fillna(nums[col].median())
nums['is_complete'] = df['is_complete'].fillna(False).astype(int)

X = pd.concat([ec_dummies, gen_dummies, nums], axis=1)
print(f"Feature matrix: {X.shape}", file=sys.stderr)

# ── RandomForest with cross-validation ──────────────────────────────
rf = RandomForestRegressor(n_estimators=300, max_depth=10, min_samples_leaf=10,
                            random_state=42, n_jobs=-1)
print("Running 5-fold CV...", file=sys.stderr)
cv_r2 = cross_val_score(rf, X, y, cv=5, scoring='r2', n_jobs=-1)
print(f"  Mean R²: {cv_r2.mean():.4f} +/- {cv_r2.std():.4f}", file=sys.stderr)

# Fit on all data for feature importance
rf.fit(X, y)
print(f"  Training R²: {rf.score(X, y):.4f}", file=sys.stderr)

# Feature importance
imp = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False)
imp.to_csv(OUT_T, sep='\t', header=['importance'])
print(f"\nTop 25 features:", file=sys.stderr)
print(imp.head(25).to_string(), file=sys.stderr)

# Plot top 25
top_n = 25
fig, ax = plt.subplots(figsize=(10, 8))
top = imp.head(top_n)[::-1]
colors = ['seagreen' if f.startswith('ec_') else
          'steelblue' if f.startswith('gen_') else 'salmon' for f in top.index]
ax.barh(range(len(top)), top.values, color=colors)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top.index, fontsize=9)
ax.set_xlabel("Random Forest feature importance")
ax.set_title(f"Top {top_n} Km drivers (5-fold CV R² = {cv_r2.mean():.3f})")
# Custom legend
from matplotlib.patches import Patch
legend_elems = [
    Patch(facecolor='seagreen', label='EC class'),
    Patch(facecolor='steelblue', label='Genus'),
    Patch(facecolor='salmon', label='Sample / sequence metadata'),
]
ax.legend(handles=legend_elems, loc='lower right')
plt.tight_layout()
plt.savefig(OUT_F, dpi=120)
plt.close()
print(f"\nSaved {OUT_F}", file=sys.stderr)

# Summary by feature category
print(f"\nAggregate importance by feature category:", file=sys.stderr)
cat_sums = pd.Series({
    'EC class (one-hot)':     imp[imp.index.str.startswith('ec_')].sum(),
    'Genus (one-hot, top 30)': imp[imp.index.str.startswith('gen_')].sum(),
    'Latitude':                imp.get('latitude', 0),
    'Longitude':               imp.get('longitude', 0),
    'GCfrac':                  imp.get('GCfrac', 0),
    'seq_len':                 imp.get('seq_len', 0),
    'n_carboxy_pfam':          imp.get('n_carboxy_pfam', 0),
    'is_complete':             imp.get('is_complete', 0),
}).sort_values(ascending=False)
print(cat_sums, file=sys.stderr)

print("\nDone.", file=sys.stderr)
