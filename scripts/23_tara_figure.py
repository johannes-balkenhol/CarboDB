"""
Script 23: TARA Oceans Metagenome Results Figure
Run: python scripts/23_tara_figure.py
Output: figures/tara_metagenome_figure.pdf/.png
"""
import json, glob, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from collections import Counter
import warnings; warnings.filterwarnings('ignore')

# ── Load results ──────────────────────────────────────────────────────────
results = []
for f in sorted(glob.glob('data/metagenome/results_random/chunk_*.json')):
    d = json.load(open(f))
    results.extend(d if isinstance(d, list) else [d])

meta = {
    'MGYA00679207': {'temp': None,  'depth': 40,  'ocean': 'N Pacific\n(40m)'},
    'MGYA00679210': {'temp': 18.17, 'depth': 250, 'ocean': 'N Atlantic\n(250m deep)'},
    'MGYA00679214': {'temp': 7.35,  'depth': 5,   'ocean': 'Southern Ocean\n(7.4°C)'},
    'MGYA00679222': {'temp': 21.47, 'depth': 5,   'ocean': 'Mediterranean\n(21.5°C)'},
    'MGYA00679225': {'temp': 25.15, 'depth': 5,   'ocean': 'S Pacific\n(25.2°C)'},
}
sample_colors = {
    'MGYA00679207': '#1565C0',
    'MGYA00679210': '#0097A7',
    'MGYA00679214': '#2E7D32',
    'MGYA00679222': '#E65100',
    'MGYA00679225': '#B71C1C',
}

plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 7,
    'axes.linewidth': 0.6, 'axes.spines.right': False,
    'axes.spines.top': False, 'figure.dpi': 300,
    'axes.titlesize': 8, 'axes.titleweight': 'bold',
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
})
BG = '#FAFAFA'

fig = plt.figure(figsize=(11.69, 8.27))
fig.suptitle('CarboDB v5 — TARA Oceans metagenome scan: carboxylase Km landscape from 5 ocean samples',
             fontsize=8.5, fontweight='bold', y=0.98)
gs = gridspec.GridSpec(2, 3, figure=fig, left=0.08, right=0.97,
                       top=0.91, bottom=0.10, wspace=0.42, hspace=0.55)
axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]

def style(ax):
    ax.set_facecolor(BG)
    ax.grid(axis='both', color='#e8e8e8', linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)

def lab(ax, l):
    ax.text(-0.16, 1.10, l, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')

# ── A: Overall EC distribution ────────────────────────────────────────────
ax = axes[0]; lab(ax, 'A'); style(ax)
ax.set_title('EC class distribution (CO₂-interacting)')
co2 = [r for r in results if r.get('is_carboxylase')]
ec_counts = Counter(r['ec_predicted'] for r in co2)
EC_NAMES = {
    '4.2.1.1': 'Carbonic\nanhydrase',
    '4.1.1.39': 'RuBisCO',
    '6.3.5.5': 'CPS',
    '6.3.3.3': 'PFAS',
    '4.1.1.112': 'β-CA',
    '6.3.4.14': 'BC',
    '4.1.1.21': 'PGD',
    '6.3.4.16': 'CbaA',
}
EC_COLORS = {
    '4.2.1.1': '#1565C0', '4.1.1.39': '#2E7D32', '6.3.5.5': '#E65100',
    '6.3.3.3': '#7B1FA2', '4.1.1.112': '#0097A7', '6.3.4.14': '#F57F17',
    '4.1.1.21': '#37474F', '6.3.4.16': '#AD1457',
}
top_ecs = ec_counts.most_common(8)
labels = [EC_NAMES.get(ec, ec) for ec, _ in top_ecs]
vals = [n for _, n in top_ecs]
colors = [EC_COLORS.get(ec, '#999999') for ec, _ in top_ecs]
bars = ax.bar(range(len(vals)), vals, color=colors, alpha=0.85, linewidth=0)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=5.5, rotation=30, ha='right')
ax.set_ylabel('Count')
for i, v in enumerate(vals):
    ax.text(i, v + 0.3, str(v), ha='center', fontsize=5.5)
ax.text(0.97, 0.97, f'Total CO₂-interacting:\n{len(co2)}/10,000 (4.5%)',
        transform=ax.transAxes, ha='right', va='top', fontsize=5.5,
        style='italic', color='#333333')

# ── B: RuBisCO Km distribution ────────────────────────────────────────────
ax = axes[1]; lab(ax, 'B'); style(ax)
ax.set_title('RuBisCO (4.1.1.39) predicted Km distribution')
rubisco = [r for r in results if r.get('km_predicted_mM') and r['ec_predicted'] == '4.1.1.39']
kms = [r['km_predicted_mM'] * 1000 for r in rubisco]
log_kms = np.log10(kms)
ax.hist(log_kms, bins=15, color='#2E7D32', alpha=0.8, edgecolor='white', linewidth=0.4)
ax.set_xlabel('Predicted Km (log₁₀ µM)')
ax.set_ylabel('Count')
ax.axvline(np.log10(10.5), color='#E53935', linewidth=1.5, linestyle='--',
           label='Top hit 10.5 µM\n(Holococcolithophora)')
ax.axvline(np.log10(40), color='#FF8F00', linewidth=1.0, linestyle=':',
           label='Known coccolithophore\n~40 µM (Emiliania)')
ax.legend(fontsize=5, loc='upper right')
ax.text(0.02, 0.97, f'n={len(rubisco)}\nmean={np.mean(kms):.0f} µM\nrange={min(kms):.0f}–{max(kms):.0f} µM',
        transform=ax.transAxes, va='top', fontsize=5.5, style='italic')
ticks = [10, 50, 100, 500, 1000, 2000]
ax.set_xticks([np.log10(t) for t in ticks])
ax.set_xticklabels([str(t) for t in ticks])

# ── C: Per-sample RuBisCO Km ─────────────────────────────────────────────
ax = axes[2]; lab(ax, 'C'); style(ax)
ax.set_title('RuBisCO Km per ocean sample')
sample_order = ['MGYA00679214', 'MGYA00679210', 'MGYA00679207', 'MGYA00679222', 'MGYA00679225']
positions = []
all_kms = []
sample_labels = []
for i, s in enumerate(sample_order):
    sample_rubisco = [r for r in rubisco if r['cdb_query_id'].startswith(s)]
    if sample_rubisco:
        skms = [r['km_predicted_mM'] * 1000 for r in sample_rubisco]
        all_kms.append(skms)
        positions.append(i)
        sample_labels.append(meta[s]['ocean'])

bp = ax.boxplot(all_kms, positions=positions, widths=0.5,
                patch_artist=True, showfliers=True,
                medianprops={'color': 'white', 'linewidth': 1.5},
                flierprops={'marker': 'o', 'markersize': 3, 'alpha': 0.5})
for i, (patch, s) in enumerate(zip(bp['boxes'], sample_order)):
    patch.set_facecolor(sample_colors[s])
    patch.set_alpha(0.8)
ax.set_xticks(positions)
ax.set_xticklabels(sample_labels, fontsize=5.5)
ax.set_ylabel('Predicted Km (µM)')
# Highlight top hit
ax.axhline(10.5, color='#E53935', linewidth=1.0, linestyle='--', alpha=0.7)
ax.text(0.97, 0.06, '10.5 µM\nHolococcolithophora',
        transform=ax.transAxes, ha='right', fontsize=5, color='#E53935', style='italic')
# Add temperature annotations
temps = [7.35, 18.17, None, 21.47, 25.15]
for i, (pos, t) in enumerate(zip(positions, temps)):
    if t:
        ax.text(pos, ax.get_ylim()[0] - ax.get_ylim()[1] * 0.05,
                f'{t}°C', ha='center', fontsize=4.8, color='#555555')

# ── D: CA Km per sample ───────────────────────────────────────────────────
ax = axes[3]; lab(ax, 'D'); style(ax)
ax.set_title('Carbonic anhydrase Km per ocean sample')
ca_seqs = [r for r in results if r.get('km_predicted_mM') and r['ec_predicted'] == '4.2.1.1']
ca_kms_per_sample = []
ca_positions = []
for i, s in enumerate(sample_order):
    sca = [r for r in ca_seqs if r['cdb_query_id'].startswith(s)]
    if sca:
        ca_kms_per_sample.append([r['km_predicted_mM'] * 1000 for r in sca])
        ca_positions.append(i)

bp2 = ax.boxplot(ca_kms_per_sample, positions=ca_positions, widths=0.5,
                 patch_artist=True, showfliers=False,
                 medianprops={'color': 'white', 'linewidth': 1.5})
for patch, s in zip(bp2['boxes'], sample_order):
    patch.set_facecolor(sample_colors[s])
    patch.set_alpha(0.8)
ax.set_xticks(ca_positions)
ax.set_xticklabels([meta[s]['ocean'] for s in sample_order], fontsize=5.5)
ax.set_ylabel('Predicted Km (µM)')

# ── E: Temperature vs RuBisCO Km ──────────────────────────────────────────
ax = axes[4]; lab(ax, 'E'); style(ax)
ax.set_title('Temperature vs mean RuBisCO Km (n=4 samples)')
for s, m in meta.items():
    if m['temp']:
        srb = [r for r in rubisco if r['cdb_query_id'].startswith(s)]
        if srb:
            mean_km = np.mean([r['km_predicted_mM'] * 1000 for r in srb])
            ax.scatter(m['temp'], mean_km, color=sample_colors[s], s=80,
                       zorder=5, label=m['ocean'].replace('\n', ' '))
            ax.annotate(m['ocean'].replace('\n', ' '),
                        (m['temp'], mean_km), fontsize=4.8,
                        xytext=(4, 2), textcoords='offset points')

ax.set_xlabel('Ocean temperature (°C)')
ax.set_ylabel('Mean RuBisCO Km (µM)')
ax.text(0.05, 0.97, 'r=0.039, p=0.96 (n=4, underpowered)\nNeed 50+ samples for significance',
        transform=ax.transAxes, va='top', fontsize=5, style='italic', color='#E53935')

# ── F: Top hit spotlight ──────────────────────────────────────────────────
ax = axes[5]; lab(ax, 'F'); style(ax)
ax.set_title('Top hit: Holococcolithophora RuBisCO (Mediterranean)')
ax.axis('off')

info = [
    ('Sequence ID', 'ERZ17499738.381682-NODE-1124297'),
    ('Sample', 'MGYA00679222 — Mediterranean'),
    ('Temperature', '21.47°C (warm surface water)'),
    ('Depth', '5 m'),
    ('Predicted Km', '10.5 µM ← lowest RuBisCO in scan'),
    ('Binary probability', '0.995 (high confidence)'),
    ('EC confidence', '0.997 (4.1.1.39 RuBisCO)'),
    ('Pfam', 'PF00016 (RuBisCO large subunit N-term)'),
    ('Fragment length', '194 / 459 aa (42% of full protein)'),
    ('BLAST top hit', 'P48687.1 — Holococcolithophora sphaeroidea'),
    ('BLAST identity', '89% (E=2×10⁻¹³⁰)'),
    ('Known coccolithophore Km', '~40 µM (Emiliania huxleyi)'),
    ('Fold difference', '3–7× lower than published values'),
    ('Priority', 'HIGH — warrants experimental validation'),
]
y = 0.97
for k, v in info:
    bold = k in ('Predicted Km', 'BLAST top hit', 'Priority')
    ax.text(0.02, y, f'{k}:', transform=ax.transAxes,
            fontsize=5.8, fontweight='bold' if bold else 'normal',
            va='top', color='#1A5276')
    ax.text(0.42, y, v, transform=ax.transAxes,
            fontsize=5.8, va='top',
            color='#E53935' if bold else '#1a1a1a',
            fontweight='bold' if bold else 'normal')
    y -= 0.065

from pathlib import Path
Path('figures').mkdir(exist_ok=True)
fig.savefig('figures/tara_metagenome_figure.pdf', format='pdf', bbox_inches='tight')
fig.savefig('figures/tara_metagenome_figure.png', format='png', bbox_inches='tight', dpi=300)
print('Saved: figures/tara_metagenome_figure.pdf/.png')
plt.close()
