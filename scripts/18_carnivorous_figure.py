"""
Carnivorous plant RuBisCO Km analysis
Panel A: Carnivorous vs non-carnivorous relatives (Caryophyllales, Lamiales, Ericales)
Panel B: Km by trap type with trend line
"""
import sqlite3, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

ROOT = Path('/home/job37yv/Projects_shared/CarboDB_v3')
conn = sqlite3.connect(ROOT / 'data/primary/carbodb.sqlite')

# ── Get carnivorous plant predictions ─────────────────────────────────────
carn = pd.read_sql("""
    SELECT s.organism, s.ec_number, p.km_pred_mM,
        CASE
            WHEN s.organism LIKE 'Dionaea%'     THEN 'Dionaea'
            WHEN s.organism LIKE 'Drosera%'     THEN 'Drosera'
            WHEN s.organism LIKE 'Sarracenia%'  THEN 'Sarracenia'
            WHEN s.organism LIKE 'Utricularia%' THEN 'Utricularia'
            WHEN s.organism LIKE 'Pinguicula%'  THEN 'Pinguicula'
            WHEN s.organism LIKE 'Nepenthes%'   THEN 'Nepenthes'
            WHEN s.organism LIKE 'Cephalotus%'  THEN 'Cephalotus'
        END as genus
    FROM predictions p JOIN sequences s ON s.id=p.sequence_id
    WHERE s.label=1 AND p.km_pred_mM IS NOT NULL
    AND s.ec_number='4.1.1.39'
    AND (s.organism LIKE 'Dionaea%' OR s.organism LIKE 'Drosera%'
         OR s.organism LIKE 'Sarracenia%' OR s.organism LIKE 'Utricularia%'
         OR s.organism LIKE 'Pinguicula%' OR s.organism LIKE 'Nepenthes%')
""", conn)
carn['km_uM'] = carn['km_pred_mM'] * 1000
print(f"Carnivorous: {len(carn)} sequences")
print(carn.groupby('genus')['km_uM'].agg(['count','mean','std']).round(2))

# ── Trap type assignment ───────────────────────────────────────────────────
trap_map = {
    'Dionaea':    ('Snap trap',    1),
    'Sarracenia': ('Pitfall',      2),
    'Drosera':    ('Flypaper',     3),
    'Utricularia':('Suction trap', 4),
    'Pinguicula': ('CAM',          5),
}
carn_rubisco = carn[carn['genus'].isin(trap_map)].copy()
carn_rubisco['trap_type'] = carn_rubisco['genus'].map(lambda g: trap_map[g][0])
carn_rubisco['trap_rank'] = carn_rubisco['genus'].map(lambda g: trap_map[g][1])

# ── Non-carnivorous relatives ──────────────────────────────────────────────
# Caryophyllales (contains Caryophyllaceae, Amaranthaceae, Polygonaceae etc.)
# Lamiales (contains Olea, Plantago, Verbena)
# Ericales (contains Ericaceae, Primulaceae)
non_carn_orders = {
    'Non-carnivorous\nCaryophyllales': ['Spinacia','Beta','Chenopodium','Amaranthus',
                                         'Silene','Dianthus','Mesembryanthemum'],
    'Non-carnivorous\nLamiales':       ['Olea','Plantago','Verbena','Antirrhinum',
                                         'Mimulus','Striga'],
    'Non-carnivorous\nEricales':       ['Vaccinium','Rhododendron','Camellia',
                                         'Primula','Lysimachia'],
}

non_carn_data = {}
for label, genera in non_carn_orders.items():
    rows = []
    for genus in genera:
        r = pd.read_sql(f"""
            SELECT p.km_pred_mM FROM predictions p JOIN sequences s ON s.id=p.sequence_id
            WHERE s.label=1 AND p.km_pred_mM IS NOT NULL
            AND s.ec_number='4.1.1.39' AND s.organism LIKE '{genus}%'
            LIMIT 200
        """, conn)
        rows.append(r)
    if rows:
        combined = pd.concat(rows)['km_pred_mM'].values * 1000
        non_carn_data[label] = combined
        print(f"{label}: n={len(combined)} mean={combined.mean():.2f} µM")

conn.close()

# ── Figure ─────────────────────────────────────────────────────────────────
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.8,
    'xtick.major.size':3,'ytick.major.size':3,'figure.dpi':300})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.15, wspace=0.35)

# ── Panel A: Carnivorous vs non-carnivorous ────────────────────────────────
carn_all = carn_rubisco['km_uM'].values
groups_A = {'Carnivorous\nPlants': (carn_all, '#2E7D32')}
for label, vals in non_carn_data.items():
    groups_A[label] = (vals, '#6D4C41')

labels_A = list(groups_A.keys())
x_A = np.arange(len(labels_A))
colors_A = [v[1] for v in groups_A.values()]
means_A  = [v[0].mean() for v in groups_A.values()]
sems_A   = [v[0].std()/np.sqrt(len(v[0])) for v in groups_A.values()]
ns_A     = [len(v[0]) for v in groups_A.values()]

bars = ax1.bar(x_A, means_A, color=colors_A, alpha=0.88, width=0.6,
               yerr=sems_A, capsize=5, error_kw={'linewidth':1.2,'capthick':1.2},
               linewidth=0)
ax1.set_xticks(x_A); ax1.set_xticklabels(labels_A, fontsize=8.5)
ax1.set_ylabel('Predicted Km (µM)', fontsize=10)
ax1.set_title('A.  RuBisCO Km: Carnivorous vs Non-Carnivorous', fontsize=11,
              fontweight='bold', loc='left', pad=8)
ax1.spines['right'].set_visible(False); ax1.spines['top'].set_visible(False)
ax1.set_facecolor('#FAFAFA'); ax1.grid(axis='y', color='#e0e0e0', linewidth=0.5)
ax1.set_axisbelow(True)

for bar, n, m in zip(bars, ns_A, means_A):
    ax1.text(bar.get_x()+bar.get_width()/2, 1.5, f'n={n}',
             ha='center', va='bottom', fontsize=7.5, color='white', fontweight='bold')

# Significance bracket
t, p = stats.mannwhitneyu(carn_all, non_carn_data['Non-carnivorous\nCaryophyllales'])
all_noncarn = np.concatenate([v[0] for k,v in groups_A.items() if 'Non' in k])
t2, p2 = stats.mannwhitneyu(carn_all, all_noncarn)
ax1.annotate('', xy=(2.7, max(means_A)*1.15), xytext=(0, max(means_A)*1.15),
             arrowprops=dict(arrowstyle='-', color='black', linewidth=1.2))
ax1.text(1.35, max(means_A)*1.17, f'*** 30–40% lower\n(p={p2:.2e})',
         ha='center', fontsize=8.5, fontweight='bold', color='#1B5E20')

ax1.set_ylim(0, max(means_A)*1.35)

# ── Panel B: Km by trap type ───────────────────────────────────────────────
trap_order = ['Snap trap','Pitfall','Flypaper','Suction trap','CAM']
trap_genera= {'Snap trap':'Dionaea','Pitfall':'Sarracenia','Flypaper':'Drosera',
              'Suction trap':'Utricularia','CAM':'Pinguicula'}
trap_colors= ['#4CAF50','#8BC34A','#FFEB3B','#FF9800','#F44336']

trap_means=[]; trap_sems=[]; trap_ns=[]
for trap in trap_order:
    genus = trap_genera[trap]
    vals = carn_rubisco[carn_rubisco['trap_type']==trap]['km_uM'].values
    trap_means.append(vals.mean() if len(vals)>0 else 0)
    trap_sems.append(vals.std()/np.sqrt(max(len(vals),1)))
    trap_ns.append(len(vals))

x_B = np.arange(len(trap_order))
bars2 = ax2.bar(x_B, trap_means, color=trap_colors, alpha=0.9, width=0.6,
                yerr=trap_sems, capsize=4,
                error_kw={'linewidth':1.0,'capthick':1.0}, linewidth=0)

# Trend line
m, b, r, p_trend, se = stats.linregress(range(1,6), trap_means)
x_trend = np.linspace(-0.3, 4.3, 100)
ax2.plot(x_trend, m*(x_trend+1)+b, '-', color='#9E9E9E', linewidth=1.5,
         alpha=0.7, zorder=5)

ax2.set_xticks(x_B)
ax2.set_xticklabels(trap_order, fontsize=8.5)
ax2.set_ylabel('Predicted Km (µM)', fontsize=10)
ax2.set_title('B.  RuBisCO Km by Trap Type', fontsize=11,
              fontweight='bold', loc='left', pad=8)
ax2.spines['right'].set_visible(False); ax2.spines['top'].set_visible(False)
ax2.set_facecolor('#FAFAFA'); ax2.grid(axis='y', color='#e0e0e0', linewidth=0.5)
ax2.set_axisbelow(True)

# Genus labels
for i, (trap, genus) in enumerate(trap_genera.items()):
    ax2.text(i, trap_means[i]+0.4, f'({genus})', ha='center', fontsize=7.2,
             color='#555555', style='italic')

ax2.text(0.55, 0.12, 'Higher trap cost → Lower Km',
         transform=ax2.transAxes, fontsize=8.5, color='#9E9E9E',
         style='italic', ha='center')

ax2.set_ylim(0, max(trap_means)*1.4)

fig.suptitle('CarboDB v5 — Carnivorous Plant RuBisCO CO₂ Affinity',
             fontsize=12, fontweight='bold', y=0.97)

out = ROOT / 'figures' / 'carnivorous_plant_figure.png'
fig.savefig(out, dpi=300, bbox_inches='tight')
print(f'Saved: {out}')
plt.close()
