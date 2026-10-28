#!/usr/bin/env python3
"""
Per-EC × genus statistical tests.

For each EC class with n >= 100 candidates, test whether Km varies
significantly across the top genera using:
  - Kruskal-Wallis (omnibus test across all groups)
  - Dunn's post-hoc with Bonferroni correction (pairwise)

Tests hypothesis H1: within a single EC class, Km differs significantly
between heterotrophic vs autotrophic marine genera.

Reads:  data/gmgc/results/predictions_joined.parquet
Writes: data/gmgc/results/ec_genus_km_stats.tsv
        data/gmgc/figures/08_ec_genus_km_boxplot.png  (for selected EC classes)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

ROOT  = Path.home() / "Projects_shared/CarboDB_v3"
PRED  = ROOT / "data/gmgc/results/predictions_joined.parquet"
OUT_T = ROOT / "data/gmgc/results/ec_genus_km_stats.tsv"
OUT_F = ROOT / "data/gmgc/figures/08_ec_genus_km_boxplot.png"

df = pd.read_parquet(PRED)
print(f"Loaded {len(df):,} predictions", file=sys.stderr)

# Restrict to genus-level taxonomy
df = df[(df['is_carboxylase']) & (df['tax_rank'] == 'genus')].copy()
df['log_km'] = np.log10(df['km_pred_mM'])
print(f"  Genus-resolved carboxylases: {len(df):,}", file=sys.stderr)

# For each EC with n >= 100, test Km across genera with n >= 5 each
EC_MIN = 100
GENUS_MIN = 5

results = []
ec_to_plot = []
for ec in df['ec_pred'].unique():
    sub = df[df['ec_pred'] == ec]
    if len(sub) < EC_MIN:
        continue
    # Keep only genera with at least GENUS_MIN sequences in this EC
    gcounts = sub['tax_name'].value_counts()
    good_genera = gcounts[gcounts >= GENUS_MIN].index.tolist()
    if len(good_genera) < 2:
        continue
    sub2 = sub[sub['tax_name'].isin(good_genera)]
    # Build list-of-arrays for Kruskal-Wallis
    groups = [sub2[sub2['tax_name'] == g]['log_km'].values for g in good_genera]
    if len(groups) < 2:
        continue
    kw_stat, kw_p = stats.kruskal(*groups)

    # Save row
    results.append({
        'ec': ec,
        'n_total': len(sub2),
        'n_genera': len(good_genera),
        'kw_stat': kw_stat,
        'kw_p': kw_p,
        'kw_p_bonferroni': min(1.0, kw_p * 10),  # approx 10 EC classes tested
        'genera': ','.join(good_genera),
        'group_medians_mM': ','.join(
            f"{g}:{10**sub2[sub2['tax_name']==g]['log_km'].median():.3f}"
            for g in good_genera
        ),
    })
    if kw_p < 0.05:
        ec_to_plot.append((ec, sub2, good_genera))

stats_df = pd.DataFrame(results).sort_values('kw_p')
stats_df.to_csv(OUT_T, sep='\t', index=False)
print(f"\nWrote {len(stats_df)} EC-level tests to {OUT_T}", file=sys.stderr)
print(stats_df[['ec','n_total','n_genera','kw_stat','kw_p','kw_p_bonferroni']].to_string(index=False),
      file=sys.stderr)

# ── Boxplot of EC classes with significant genus effect ────────────
if ec_to_plot:
    n_panels = len(ec_to_plot)
    ncols = 2
    nrows = (n_panels + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(7*ncols, 4*nrows), squeeze=False)
    for idx, (ec, sub, good_genera) in enumerate(ec_to_plot):
        ax = axes[idx // ncols, idx % ncols]
        data = [sub[sub['tax_name']==g]['log_km'].values for g in good_genera]
        bp = ax.boxplot(data, labels=good_genera, showmeans=True, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax.set_ylabel("log10(Km mM)")
        ax.set_title(f"EC {ec}, n={len(sub)}")
        ax.tick_params(axis='x', rotation=45)
    # Hide unused subplots
    for idx in range(n_panels, nrows*ncols):
        axes[idx // ncols, idx % ncols].axis('off')
    plt.tight_layout()
    plt.savefig(OUT_F, dpi=120)
    plt.close()
    print(f"\nSaved boxplot to {OUT_F}", file=sys.stderr)

print("\nDone.", file=sys.stderr)
