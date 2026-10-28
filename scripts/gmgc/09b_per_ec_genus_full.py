#!/usr/bin/env python3
"""Per-EC × genus Kruskal-Wallis at full marine scale."""
import sys
from pathlib import Path
import pandas as pd
from scipy.stats import kruskal

ROOT = Path.home() / "Projects_shared/CarboDB_v3"
df = pd.read_parquet(ROOT / "data/gmgc/results/predictions_joined_full.parquet")
df = df[df['is_carboxylase'] & (df['habitat']=='marine')].copy()

# Restrict to high-confidence ec_pred >= 0.7 for cleaner biology
df = df[df['ec_conf'] >= 0.7]
print(f"After ec_conf>=0.7 filter: {len(df):,}", file=sys.stderr)

results = []
for ec, sub in df.groupby('ec_pred'):
    if len(sub) < 100: continue
    # Top genera with n>=30 in this EC
    g = sub[sub['tax_rank']=='genus'].groupby('tax_name').filter(lambda x: len(x)>=30)
    if g['tax_name'].nunique() < 2: continue
    groups = [grp['km_pred_mM'].values for _, grp in g.groupby('tax_name')]
    if min(len(x) for x in groups) < 30: continue
    try:
        stat, p = kruskal(*groups)
    except Exception:
        continue
    medians = g.groupby('tax_name')['km_pred_mM'].median().sort_values()
    results.append({
        'ec': ec, 'n_total': len(sub), 'n_genera_tested': g['tax_name'].nunique(),
        'kw_stat': stat, 'kw_p': p,
        'group_medians': '; '.join(f'{k}:{v:.3f}' for k,v in medians.items())
    })

res = pd.DataFrame(results).sort_values('kw_p')
# Bonferroni correction
res['kw_p_bonferroni'] = res['kw_p'] * len(res)
res['kw_p_bonferroni'] = res['kw_p_bonferroni'].clip(upper=1.0)
res['significant'] = res['kw_p_bonferroni'] < 0.05

print(res[['ec','n_total','n_genera_tested','kw_stat','kw_p','kw_p_bonferroni','significant']]
      .to_string(index=False), file=sys.stderr)
res.to_csv(ROOT / "data/gmgc/results/ec_genus_km_stats_full.tsv", sep='\t', index=False)
print(f"\nWrote {ROOT}/data/gmgc/results/ec_genus_km_stats_full.tsv", file=sys.stderr)
