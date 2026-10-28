#!/usr/bin/env python3
"""Carboxylase composition by ocean basin (Tara-style biogeography)."""
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path.home() / "Projects_shared/CarboDB_v3"
df = pd.read_parquet(ROOT / "data/gmgc/results/predictions_joined_full.parquet")
df = df[df['is_carboxylase'] & (df['habitat']=='marine')].copy()

# Simple ocean basin classifier from lat/lon
def basin(lat, lon):
    if pd.isna(lat) or pd.isna(lon):
        return 'unknown'
    if lat < -60:                                  return 'Southern Ocean'
    if lat > 60 and (lon < -100 or lon > 100):    return 'Arctic'
    if 30 <= lat <= 47 and -10 <= lon <= 40:      return 'Mediterranean'
    if -60 < lat < 30 and 40 < lon < 110:         return 'Indian Ocean'
    if lon < -70 or (lon > 20 and lon < 25):      return 'Atlantic'
    if lon > 110 or lon < -100:                    return 'Pacific'
    return 'Atlantic'

df['ocean_basin'] = df.apply(lambda r: basin(r['latitude'], r['longitude']), axis=1)

print("═══ Sample distribution by basin ═══")
basin_samples = df.groupby('ocean_basin')['sample_id'].nunique().sort_values(ascending=False)
print(basin_samples)
print()

print("═══ Carboxylase sequence count by basin ═══")
basin_counts = df.groupby('ocean_basin').size().sort_values(ascending=False)
print(basin_counts)
print()

print("═══ EC composition by basin (top ECs, fraction of basin total) ═══")
top_ecs = df['ec_pred'].value_counts().head(8).index.tolist()
basin_ec = df[df['ec_pred'].isin(top_ecs)].groupby(['ocean_basin','ec_pred']).size().unstack(fill_value=0)
basin_ec_frac = basin_ec.div(basin_ec.sum(axis=1), axis=0)
print((basin_ec_frac * 100).round(1))
print()

print("═══ RuBisCO Form distribution by basin ═══")
rbc = pd.read_csv(ROOT / "data/gmgc/results/rubisco_forms_v2.tsv", sep='\t')
gps = df[['seq_id','latitude','longitude']].drop_duplicates('seq_id')
rbc = rbc.merge(gps, on='seq_id', how='left').dropna(subset=['latitude'])
rbc['ocean_basin'] = rbc.apply(lambda r: basin(r['latitude'], r['longitude']), axis=1)
rbc_basin = rbc.groupby(['ocean_basin','rubisco_form']).size().unstack(fill_value=0)
print(rbc_basin)
print()

print("═══ Pelagibacter CA Km by basin ═══")
pel_ca = df[(df['tax_name']=='Candidatus Pelagibacter') & (df['ec_pred']=='4.2.1.1')]
print(pel_ca.groupby('ocean_basin')['km_pred_mM'].agg(['count','median','mean']))
print()

# Stacked bar — EC composition by basin
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
ec_pct = (basin_ec_frac * 100)
ec_pct.plot(kind='bar', stacked=True, ax=axes[0], colormap='Set2', width=0.8)
axes[0].set_ylabel('% of basin carboxylases')
axes[0].set_xlabel('Ocean basin')
axes[0].set_title(f'EC composition by ocean basin (top 8 ECs)')
axes[0].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8)
axes[0].tick_params(axis='x', rotation=30)

# RuBisCO forms by basin
rbc_pct = rbc_basin.div(rbc_basin.sum(axis=1), axis=0) * 100
rbc_pct.plot(kind='bar', stacked=True, ax=axes[1], colormap='tab10', width=0.8)
axes[1].set_ylabel('% of basin RuBisCOs')
axes[1].set_xlabel('Ocean basin')
axes[1].set_title('RuBisCO Form composition by ocean basin')
axes[1].legend(loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=8)
axes[1].tick_params(axis='x', rotation=30)

plt.tight_layout()
plt.savefig(ROOT / "data/gmgc/figures/15_ocean_basin_composition.png",
             dpi=120, bbox_inches='tight')
print("Saved figure", file=sys.stderr)
