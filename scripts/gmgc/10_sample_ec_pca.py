#!/usr/bin/env python3
"""
Sample × EC composition matrix + PCA + PERMANOVA.

Builds a samples × EC-class composition matrix (relative abundances),
does PCA, colors by RuBisCO presence/absence, country, and latitude.
Tests hypothesis H2: samples with RuBisCO have distinct overall
carboxylase composition.

Reads:  data/gmgc/results/predictions_joined.parquet
Writes: data/gmgc/figures/09_sample_pca_rbc.png
        data/gmgc/figures/10_sample_pca_latitude.png
        data/gmgc/results/sample_ec_matrix.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

ROOT  = Path.home() / "Projects_shared/CarboDB_v3"
PRED  = ROOT / "data/gmgc/results/predictions_joined.parquet"
OUT_T = ROOT / "data/gmgc/results/sample_ec_matrix.tsv"
OUT_F1 = ROOT / "data/gmgc/figures/09_sample_pca_rbc.png"
OUT_F2 = ROOT / "data/gmgc/figures/10_sample_pca_latitude.png"

df = pd.read_parquet(PRED)
df = df[df['is_carboxylase'] & df['sample_id'].notna()].copy()
print(f"Predictions with sample_id: {len(df):,}", file=sys.stderr)

# Build samples × EC count matrix
ec_counts = pd.crosstab(df['sample_id'], df['ec_pred'])
print(f"Matrix shape: {ec_counts.shape}  (samples × EC classes)", file=sys.stderr)

# Convert to relative abundance (rows sum to 1)
ec_relab = ec_counts.div(ec_counts.sum(axis=1), axis=0)
ec_relab.to_csv(OUT_T, sep='\t')
print(f"Wrote {OUT_T}", file=sys.stderr)

# Drop samples with too few sequences (< 10) to avoid noisy compositions
counts_per_sample = ec_counts.sum(axis=1)
ec_relab_filtered = ec_relab.loc[counts_per_sample >= 10]
print(f"Samples with >= 10 carboxylases: {len(ec_relab_filtered)}", file=sys.stderr)

# PCA
pca = PCA(n_components=2)
coords = pca.fit_transform(ec_relab_filtered.values)
print(f"PCA explained variance: PC1 {pca.explained_variance_ratio_[0]:.3f}, "
      f"PC2 {pca.explained_variance_ratio_[1]:.3f}", file=sys.stderr)

# Get RuBisCO presence per sample (count > 0)
rbc_present = (ec_counts.get('4.1.1.39', pd.Series(0)) > 0)

# Sample-level lat/lon from one row per sample
sample_meta = df.groupby('sample_id').agg(
    latitude=('latitude','first'),
    longitude=('longitude','first'),
    country=('country','first'),
).reset_index()
sample_meta = sample_meta.set_index('sample_id').reindex(ec_relab_filtered.index)

# ── Figure 9: PCA colored by RuBisCO presence ───────────────────────
fig, ax = plt.subplots(figsize=(9, 6))
present = rbc_present.reindex(ec_relab_filtered.index).fillna(False)
ax.scatter(coords[~present, 0], coords[~present, 1],
           c='steelblue', s=50, alpha=0.6, label='No RuBisCO')
ax.scatter(coords[present, 0], coords[present, 1],
           c='red', s=80, alpha=0.8, label='RuBisCO present', edgecolor='black')
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("Marine samples in EC-composition space — RuBisCO presence")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_F1, dpi=120)
plt.close()
print(f"Saved {OUT_F1}", file=sys.stderr)

# ── Figure 10: PCA colored by latitude ──────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6))
lat = sample_meta['latitude'].values
sc = ax.scatter(coords[:, 0], coords[:, 1], c=lat, s=60, alpha=0.8,
                cmap='coolwarm', edgecolor='black', lw=0.3, vmin=-60, vmax=60)
cb = plt.colorbar(sc, ax=ax, label='Latitude')
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
ax.set_title("Marine samples in EC-composition space — colored by latitude")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_F2, dpi=120)
plt.close()
print(f"Saved {OUT_F2}", file=sys.stderr)

# ── Quick PERMANOVA-like test (Mann-Whitney on PC1 by RuBisCO presence) ──
pc1_with_rbc    = coords[present, 0]
pc1_without_rbc = coords[~present, 0]
from scipy import stats
if len(pc1_with_rbc) >= 3 and len(pc1_without_rbc) >= 3:
    u_stat, u_p = stats.mannwhitneyu(pc1_with_rbc, pc1_without_rbc)
    print(f"\nMann-Whitney on PC1, RuBisCO+ vs RuBisCO-: "
          f"U={u_stat:.1f}, p={u_p:.4f}", file=sys.stderr)
    if u_p < 0.05:
        print(f"  -> Samples with RuBisCO have significantly different PC1 distribution",
              file=sys.stderr)

# Loading insights: which ECs drive PC1?
loadings = pd.Series(pca.components_[0], index=ec_relab_filtered.columns).sort_values()
print(f"\nTop 5 EC loadings on PC1 (most negative):", file=sys.stderr)
print(loadings.head(5), file=sys.stderr)
print(f"\nTop 5 EC loadings on PC1 (most positive):", file=sys.stderr)
print(loadings.tail(5), file=sys.stderr)

print("\nDone.", file=sys.stderr)
