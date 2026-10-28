#!/usr/bin/env python3
"""Ocean basin analysis with CORRECTED basin classifier.

Fix: Original script 15 tested 'if lon < -70: return Atlantic' before
the Pacific condition (lon < -100), causing eastern Pacific samples
(longitude -100 to -159) to be lumped as Atlantic. Pacific tested FIRST.
"""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path.home() / "Projects_shared/CarboDB_v3"

def basin(lat, lon):
    """Corrected basin classifier."""
    if pd.isna(lat) or pd.isna(lon):
        return 'unknown'
    if lat < -60:                                  return 'Southern Ocean'
    if lat > 60 and (lon < -100 or lon > 100):    return 'Arctic'
    if 30 <= lat <= 47 and -10 <= lon <= 40:      return 'Mediterranean'
    if -60 < lat < 30 and 40 < lon < 110:         return 'Indian Ocean'
    # PACIFIC tested BEFORE Atlantic (the fix)
    if lon > 110 or lon < -100:                   return 'Pacific'
    if -100 <= lon < 20 or (lon > 20 and lon < 25):
                                                  return 'Atlantic'
    return 'unknown'

df = pd.read_parquet(ROOT / "data/gmgc/results/predictions_joined_full.parquet")
df = df[df['is_carboxylase'] & (df['habitat']=='marine')].copy()
df['ocean_basin'] = df.apply(lambda r: basin(r['latitude'], r['longitude']), axis=1)

print("═══ Sample distribution by basin (CORRECTED) ═══")
basin_samples = df.groupby('ocean_basin')['sample_id'].nunique().sort_values(ascending=False)
print(basin_samples)
print(f"\nTotal samples: {basin_samples.sum()}")
print()

print("═══ Carboxylase sequence count by basin (CORRECTED) ═══")
basin_counts = df.groupby('ocean_basin').size().sort_values(ascending=False)
print(basin_counts)
print(f"\nTotal carboxylases: {basin_counts.sum():,}")
print()

print("═══ EC composition by basin (CORRECTED, top ECs, % of basin total) ═══")
top_ecs = df['ec_pred'].value_counts().head(8).index.tolist()
basin_ec = df[df['ec_pred'].isin(top_ecs)].groupby(['ocean_basin','ec_pred']).size().unstack(fill_value=0)
basin_ec_frac = basin_ec.div(basin_ec.sum(axis=1), axis=0)
print((basin_ec_frac * 100).round(1))
print()

print("═══ RuBisCO Form distribution by basin (CORRECTED) ═══")
rbc = pd.read_csv(ROOT / "data/gmgc/results/rubisco_forms_v2.tsv", sep='\t')
gps = df[['seq_id','latitude','longitude']].drop_duplicates('seq_id')
rbc = rbc.merge(gps, on='seq_id', how='left').dropna(subset=['latitude','longitude'])
rbc['ocean_basin'] = rbc.apply(lambda r: basin(r['latitude'], r['longitude']), axis=1)
rbc_basin = rbc.groupby(['ocean_basin','rubisco_form']).size().unstack(fill_value=0)
print(rbc_basin)
print()

print("═══ Pelagibacter CA Km by basin (CORRECTED) ═══")
pel_ca = df[(df['tax_name']=='Candidatus Pelagibacter') & (df['ec_pred']=='4.2.1.1')]
print(pel_ca.groupby('ocean_basin')['km_pred_mM'].agg(['count','median','mean']).round(3))
print()

# Save corrected results
df[['seq_id','sample_id','tax_name','ec_pred','km_pred_mM','latitude','longitude','ocean_basin']]\
  .to_parquet(ROOT / "data/gmgc/results/predictions_with_basin_fixed.parquet")

# Stacked bar plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
(basin_ec_frac * 100).plot(kind='bar', stacked=True, ax=axes[0], colormap='Set2', width=0.8)
axes[0].set_ylabel('% of basin carboxylases')
axes[0].set_xlabel('Ocean basin (corrected)')
axes[0].set_title('EC composition by basin (Pacific reassigned)')
axes[0].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8)
axes[0].tick_params(axis='x', rotation=30)

rbc_pct = rbc_basin.div(rbc_basin.sum(axis=1), axis=0) * 100
rbc_pct.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab10', width=0.8)
axes[1].set_ylabel('% of basin RuBisCOs')
axes[1].set_xlabel('Ocean basin (corrected)')
axes[1].set_title('RuBisCO Form composition by basin (Pacific reassigned)')
axes[1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8)
axes[1].tick_params(axis='x', rotation=30)

plt.tight_layout()
plt.savefig(ROOT / "data/gmgc/figures/15b_ocean_basin_composition_fixed.png",
             dpi=120, bbox_inches='tight')
print("Saved figure: 15b_ocean_basin_composition_fixed.png", file=sys.stderr)
