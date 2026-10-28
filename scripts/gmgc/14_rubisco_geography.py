#!/usr/bin/env python3
"""Plot RuBisCO Form distribution on world map (uses rubisco_forms_v2.tsv)."""
import sys
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path.home() / "Projects_shared/CarboDB_v3"
rbc = pd.read_csv(ROOT / "data/gmgc/results/rubisco_forms_v2.tsv", sep='\t')
df = pd.read_parquet(ROOT / "data/gmgc/results/predictions_joined_full.parquet")

# Merge GPS in
gps = df[['seq_id','latitude','longitude']].drop_duplicates('seq_id')
rbc = rbc.merge(gps, on='seq_id', how='left')
rbc = rbc.dropna(subset=['latitude','longitude'])
print(f"RuBisCO predictions with GPS: {len(rbc):,}", file=sys.stderr)

# Plot
forms = sorted(rbc['rubisco_form'].unique())
fig, axes = plt.subplots(2, 2, figsize=(15, 9))
axes = axes.flatten()

# Major form categories
major_forms = ['I (cyanobacterial/canonical)', 'II (alpha-proteobacterial)',
                'III (archaeal)', 'IV (RuBisCO-like protein candidate)']

for ax, form in zip(axes, major_forms):
    sub = rbc[rbc['rubisco_form'] == form]
    ax.scatter(rbc['longitude'], rbc['latitude'], c='lightgray', s=2, alpha=0.4,
                label='All RuBisCO')
    if len(sub) > 0:
        ax.scatter(sub['longitude'], sub['latitude'], c='C3', s=8, alpha=0.7,
                    label=f'{form} (n={len(sub)})')
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude'); ax.set_ylabel('Latitude')
    ax.set_title(form.split('(')[0].strip())
    ax.legend(loc='lower left', fontsize=9)
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(ROOT / "data/gmgc/figures/14_rubisco_forms_world.png", dpi=120,
             bbox_inches='tight')
print(f"Saved figure", file=sys.stderr)

# Latitudinal density per form
fig, ax = plt.subplots(figsize=(10, 5))
for form in major_forms:
    sub = rbc[rbc['rubisco_form'] == form]
    if len(sub) >= 5:
        ax.hist(sub['latitude'], bins=18, alpha=0.5, label=f'{form.split("(")[0].strip()} (n={len(sub)})')
ax.set_xlabel('Latitude'); ax.set_ylabel('Sequence count')
ax.set_title('RuBisCO Form by latitude')
ax.legend(); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(ROOT / "data/gmgc/figures/14_rubisco_forms_latitude.png", dpi=120)
print(f"Saved latitude histogram", file=sys.stderr)
