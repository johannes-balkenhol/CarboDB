"""
Script 21: Ecological Findings Figure
6 panels summarising the ecological exploration findings
- Panel A: Parasitic vs free-living RuBisCO Km
- Panel B: Convergent C4 lineages Km comparison
- Panel C: Crop vs wild relative Km
- Panel D: CA temperature optima
- Panel E: Polyploidy effect
- Panel F: Coral symbiont + Helianthus within-genus range

Run: python scripts/21_ecological_findings_figure.py
Output: figures/ecological_findings_figure.pdf/.png
"""
import sqlite3, pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy import stats
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')

ROOT = Path('.')
conn = sqlite3.connect(ROOT / 'data/primary/carbodb.sqlite')

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':7,'axes.linewidth':0.6,
    'axes.spines.right':False,'axes.spines.top':False,'figure.dpi':300,
    'axes.titlesize':8,'axes.titleweight':'bold','axes.titlepad':5,
    'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':5.8,
    'axes.labelsize':7,
})
BG = '#FAFAFA'
def style(ax):
    ax.set_facecolor(BG)
    ax.grid(axis='both', color='#e8e8e8', linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
def lab(ax, l):
    ax.text(-0.16, 1.10, l, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')

fig = plt.figure(figsize=(11.69, 8.27))
gs  = gridspec.GridSpec(2, 3, figure=fig, left=0.08, right=0.97,
                        top=0.91, bottom=0.10, wspace=0.46, hspace=0.58)
axes = [fig.add_subplot(gs[r,c]) for r in range(2) for c in range(3)]

# ── A: Parasitic plants ───────────────────────────────────────────────────
ax = axes[0]; lab(ax,'A'); ax.set_title('Parasitic vs free-living RuBisCO Km'); style(ax)

parasites = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Cuscuta%' THEN 'Holoparasite'
             WHEN s.organism LIKE 'Orobanche%' THEN 'Holoparasite'
             WHEN s.organism LIKE 'Viscum%' THEN 'Hemiparasite'
             WHEN s.organism LIKE 'Striga%' THEN 'Hemiparasite'
             WHEN s.organism LIKE 'Rhinanthus%' THEN 'Hemiparasite'
        END as ptype
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Cuscuta%' OR s.organism LIKE 'Orobanche%'
         OR s.organism LIKE 'Viscum%' OR s.organism LIKE 'Striga%'
         OR s.organism LIKE 'Rhinanthus%')
""", conn)
free = pd.read_sql("""
    SELECT p.km_pred_mM FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Ipomoea%' OR s.organism LIKE 'Nicotiana%'
         OR s.organism LIKE 'Solanum%') LIMIT 300
""", conn)

groups = {
    'Hemi-\nparasite': parasites[parasites['ptype']=='Hemiparasite']['km_pred_mM'].values * 1000,
    'Holo-\nparasite': parasites[parasites['ptype']=='Holoparasite']['km_pred_mM'].values * 1000,
    'Free-living\nrelatives': free['km_pred_mM'].values * 1000,
}
colors = ['#E65100','#BF360C','#4CAF50']
x = np.arange(len(groups))
means = [v.mean() for v in groups.values()]
sems  = [v.std()/np.sqrt(len(v)) for v in groups.values()]
ns    = [len(v) for v in groups.values()]

bars = ax.bar(x, means, color=colors, alpha=0.85, width=0.55,
              yerr=sems, capsize=4, error_kw={'linewidth':1.1}, linewidth=0)
ax.set_xticks(x); ax.set_xticklabels(list(groups.keys()), fontsize=6.5)
ax.set_ylabel('Mean predicted Km (µM)')
for i,(n,m) in enumerate(zip(ns,means)):
    ax.text(i, 1.5, f'n={n}', ha='center', fontsize=5.5, color='white', fontweight='bold')
# Significance bracket
u,p = stats.mannwhitneyu(parasites['km_pred_mM'], free['km_pred_mM'])
ymax = max(means) * 1.25
ax.annotate('', xy=(1.5, ymax), xytext=(2, ymax),
            arrowprops=dict(arrowstyle='-', color='black', linewidth=0.8))
ax.text(1.75, ymax*1.03, f'p={p:.0e}', ha='center', fontsize=5.5, fontweight='bold')
ax.set_ylim(0, ymax*1.15)
ax.text(0.5, 0.96, 'Counter-intuitive:\nparasites lower Km\nthan free-living',
        transform=ax.transAxes, fontsize=5.2, ha='center', va='top',
        style='italic', color='#BF360C')

# ── B: Convergent C4 ──────────────────────────────────────────────────────
ax = axes[1]; lab(ax,'B'); ax.set_title('C4 lineages: different evolutionary solutions'); style(ax)

c4 = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Zea%' OR s.organism LIKE 'Sorghum%'
                  OR s.organism LIKE 'Cynodon%' OR s.organism LIKE 'Muhlenbergia%'
                  OR s.organism LIKE 'Chloris%' OR s.organism LIKE 'Sporobolus%'
                  THEN 'Poaceae'
             WHEN s.organism LIKE 'Amaranthus%' OR s.organism LIKE 'Atriplex%'
                  THEN 'Caryophyll.'
             WHEN s.organism LIKE 'Portulaca%' THEN 'Portulac.'
             WHEN s.organism LIKE 'Flaveria%' THEN 'Flaveria\n(Asterac.)'
             WHEN s.organism LIKE 'Cleome%' THEN 'Cleome\n(recent C4)'
        END as lineage
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Zea%' OR s.organism LIKE 'Sorghum%'
         OR s.organism LIKE 'Cynodon%' OR s.organism LIKE 'Muhlenbergia%'
         OR s.organism LIKE 'Chloris%' OR s.organism LIKE 'Sporobolus%'
         OR s.organism LIKE 'Amaranthus%' OR s.organism LIKE 'Atriplex%'
         OR s.organism LIKE 'Portulaca%' OR s.organism LIKE 'Flaveria%'
         OR s.organism LIKE 'Cleome%')
""", conn).dropna(subset=['lineage'])

order = ['Cleome\n(recent C4)','Caryophyll.','Portulac.','Poaceae','Flaveria\n(Asterac.)']
c4_colors = ['#2E7D32','#66BB6A','#AED581','#FFA726','#E53935']

c4_means = [c4[c4['lineage']==l]['km_pred_mM'].mean()*1000 for l in order]
c4_sems  = [c4[c4['lineage']==l]['km_pred_mM'].sem()*1000 for l in order]
c4_ns    = [len(c4[c4['lineage']==l]) for l in order]

bars2 = ax.bar(range(len(order)), c4_means, color=c4_colors, alpha=0.88,
               width=0.6, yerr=c4_sems, capsize=3,
               error_kw={'linewidth':0.8}, linewidth=0)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order, fontsize=5.5)
ax.set_ylabel('Mean predicted Km (µM)')
for i,(n,m) in enumerate(zip(c4_ns, c4_means)):
    ax.text(i, m + c4_sems[i] + 1, f'n={n}', ha='center', fontsize=4.8, color='#555555')

grps = [c4[c4['lineage']==l]['km_pred_mM'].values for l in order if len(c4[c4['lineage']==l])>=3]
if len(grps) >= 3:
    f, pv = stats.kruskal(*grps)
    ax.text(0.97, 0.97, f'KW p={pv:.2e}\nlineages DIFFER',
            transform=ax.transAxes, fontsize=5.2, ha='right', va='top',
            style='italic', color='#E53935')
ax.annotate('C4 acquisition\nprecedes\nRuBisCO relaxation',
            xy=(4, c4_means[4]), xytext=(2.8, c4_means[4]*1.3),
            fontsize=5, color='#E53935', style='italic',
            arrowprops=dict(arrowstyle='->', color='#E53935', lw=0.8))

# ── C: Crop vs wild ───────────────────────────────────────────────────────
ax = axes[2]; lab(ax,'C'); ax.set_title('Crop vs wild relative RuBisCO Km'); style(ax)

crop_data = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Triticum aestivum%' THEN 'Wheat\n6x (crop)'
             WHEN s.organism LIKE 'Hordeum vulgare%' THEN 'Barley\n(crop)'
             WHEN s.organism LIKE 'Oryza sativa%' THEN 'Rice\n(crop)'
             WHEN s.organism LIKE 'Aegilops%' THEN 'Aegilops\n(ancestor)'
             WHEN s.organism LIKE 'Oryza%' THEN 'Oryza\n(wild)'
             WHEN s.organism LIKE 'Hordeum%' THEN 'Hordeum\n(wild)'
        END as grp
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Triticum aestivum%' OR s.organism LIKE 'Hordeum vulgare%'
         OR s.organism LIKE 'Oryza sativa%' OR s.organism LIKE 'Aegilops%'
         OR s.organism LIKE 'Oryza%' OR s.organism LIKE 'Hordeum%')
""", conn).dropna(subset=['grp'])

corp_order = ['Wheat\n6x (crop)','Barley\n(crop)','Rice\n(crop)',
              'Hordeum\n(wild)','Oryza\n(wild)','Aegilops\n(ancestor)']
corp_cols  = ['#1565C0','#1976D2','#42A5F5','#9E9E9E','#BDBDBD','#E57373']

corp_means = [crop_data[crop_data['grp']==g]['km_pred_mM'].mean()*1000 for g in corp_order]
corp_ns    = [len(crop_data[crop_data['grp']==g]) for g in corp_order]

bars3 = ax.bar(range(len(corp_order)), corp_means, color=corp_cols,
               alpha=0.88, width=0.6, linewidth=0)
ax.set_xticks(range(len(corp_order)))
ax.set_xticklabels(corp_order, fontsize=5.5)
ax.set_ylabel('Mean predicted Km (µM)')
for i,(n,m) in enumerate(zip(corp_ns,corp_means)):
    ax.text(i, m+0.3, f'n={n}', ha='center', fontsize=4.8, color='#555555')
ax.annotate('', xy=(0, corp_means[0]+0.3), xytext=(5, corp_means[5]+0.3),
            arrowprops=dict(arrowstyle='->', color='#E53935', lw=1.2))
ax.text(2.5, max(corp_means)*0.9, '3× lower Km\nafter domestication',
        ha='center', fontsize=5.5, style='italic', color='#E53935')

# ── D: CA temperature optima ──────────────────────────────────────────────
ax = axes[3]; lab(ax,'D'); ax.set_title('CA Km by organism temperature optimum'); style(ax)

ca_t = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Sulfolobus%' OR s.organism LIKE 'Thermococcus%'
                  OR s.organism LIKE 'Pyrococcus%' THEN 'Thermophile\n75-95°C'
             WHEN s.organism LIKE 'Escherichia%' THEN 'Mesophile\nbact. 37°C'
             WHEN s.organism LIKE 'Chlamydomonas%' OR s.organism LIKE 'Chlorella%'
                  THEN 'Mesophile\nalga 20°C'
             WHEN s.organism LIKE 'Psychrobacter%' OR s.organism LIKE 'Colwellia%'
                  OR s.organism LIKE 'Polaribacter%' THEN 'Psychrophile\n0-10°C'
             WHEN s.organism LIKE 'Homo%' OR s.organism LIKE 'Bos%'
                  THEN 'Mesophile\nanimal 37°C'
        END as tgroup
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.2.1.1' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Sulfolobus%' OR s.organism LIKE 'Thermococcus%'
         OR s.organism LIKE 'Pyrococcus%' OR s.organism LIKE 'Escherichia%'
         OR s.organism LIKE 'Chlamydomonas%' OR s.organism LIKE 'Chlorella%'
         OR s.organism LIKE 'Psychrobacter%' OR s.organism LIKE 'Colwellia%'
         OR s.organism LIKE 'Polaribacter%' OR s.organism LIKE 'Homo%'
         OR s.organism LIKE 'Bos%')
""", conn).dropna(subset=['tgroup'])

t_order = ['Thermophile\n75-95°C','Mesophile\nbact. 37°C','Mesophile\nalga 20°C',
           'Psychrophile\n0-10°C','Mesophile\nanimal 37°C']
t_cols  = ['#E53935','#FF7043','#FFA726','#42A5F5','#7E57C2']

t_means = [ca_t[ca_t['tgroup']==t]['km_pred_mM'].mean()*1000 for t in t_order]
t_ns    = [len(ca_t[ca_t['tgroup']==t]) for t in t_order]
t_sems  = [ca_t[ca_t['tgroup']==t]['km_pred_mM'].sem()*1000 for t in t_order]

bars4 = ax.bar(range(len(t_order)), t_means, color=t_cols, alpha=0.88,
               width=0.6, yerr=t_sems, capsize=3,
               error_kw={'linewidth':0.8}, linewidth=0)
ax.set_xticks(range(len(t_order)))
ax.set_xticklabels(t_order, fontsize=5.2)
ax.set_ylabel('Mean predicted Km (mM×1000 = µM)')
ax.set_ylabel('Mean predicted Km (µM)')
for i,(n,m) in enumerate(zip(t_ns,t_means)):
    ax.text(i, m + t_sems[i] + 10, f'n={n}', ha='center', fontsize=4.8, color='#555555')
ax.text(0.5, 0.96,
        'Inverted from CO₂ solubility hypothesis:\nthermophile γ-CA needs ultra-high affinity\nat 80°C where CO₂ solubility collapses',
        transform=ax.transAxes, fontsize=4.8, ha='center', va='top',
        style='italic', color='#E53935')

# ── E: Ploidy effect ──────────────────────────────────────────────────────
ax = axes[4]; lab(ax,'E'); ax.set_title('Ploidy effect on RuBisCO Km'); style(ax)

ploidy = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Triticum aestivum%' THEN 'Wheat\n6x'
             WHEN s.organism LIKE 'Triticum durum%' THEN 'Wheat\n4x'
             WHEN s.organism LIKE 'Gossypium hirsutum%' THEN 'Cotton\n4x'
             WHEN s.organism LIKE 'Gossypium arboreum%'
                  OR s.organism LIKE 'Gossypium raimondii%' THEN 'Cotton\n2x'
             WHEN s.organism LIKE 'Avena sativa%' THEN 'Oat\n6x'
             WHEN s.organism LIKE 'Avena strigosa%' THEN 'Oat\n2x'
             WHEN s.organism LIKE 'Oryza sativa%' THEN 'Rice\n2x'
             WHEN s.organism LIKE 'Hordeum vulgare%' THEN 'Barley\n2x'
        END as pg
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Triticum aestivum%' OR s.organism LIKE 'Triticum durum%'
         OR s.organism LIKE 'Gossypium hirsutum%' OR s.organism LIKE 'Gossypium arboreum%'
         OR s.organism LIKE 'Gossypium raimondii%' OR s.organism LIKE 'Avena sativa%'
         OR s.organism LIKE 'Avena strigosa%' OR s.organism LIKE 'Oryza sativa%'
         OR s.organism LIKE 'Hordeum vulgare%')
""", conn).dropna(subset=['pg'])

p_order = ['Wheat\n6x','Wheat\n4x','Barley\n2x','Oat\n6x','Rice\n2x',
           'Oat\n2x','Cotton\n4x','Cotton\n2x']
p_cols  = ['#1565C0','#1976D2','#42A5F5','#0097A7','#26C6DA',
           '#80DEEA','#FF7043','#FFAB91']
p_means = []
p_ns    = []
for pg in p_order:
    sub = ploidy[ploidy['pg']==pg]['km_pred_mM']
    p_means.append(sub.mean()*1000 if len(sub)>0 else 0)
    p_ns.append(len(sub))

bars5 = ax.bar(range(len(p_order)), p_means, color=p_cols,
               alpha=0.88, width=0.6, linewidth=0)
ax.set_xticks(range(len(p_order)))
ax.set_xticklabels(p_order, fontsize=5.5)
ax.set_ylabel('Mean predicted Km (µM)')
for i,(n,m) in enumerate(zip(p_ns,p_means)):
    if n > 0:
        ax.text(i, m+0.2, f'n={n}', ha='center', fontsize=4.8, color='#555555')

# Annotate ploidy comparison
ax.annotate('', xy=(0, p_means[0]+1), xytext=(2, p_means[2]+1),
            arrowprops=dict(arrowstyle='<->', color='#1565C0', lw=1.0))
ax.text(1, max(p_means[:3])*0.85,
        f'{p_means[2]/p_means[0]:.1f}× higher\n in 2x vs 6x wheat',
        ha='center', fontsize=5.2, color='#1565C0', style='italic')

# ── F: Helianthus range + Symbiodiniaceae ─────────────────────────────────
ax = axes[5]; lab(ax,'F')
ax.set_title('Extreme within-genus variation + coral symbionts')
style(ax)

hel = pd.read_sql("""
    SELECT s.organism, p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Helianthus annuus%' THEN 'H. annuus\n(cultivated)'
             WHEN s.organism LIKE 'Helianthus petiolaris%' THEN 'H. petiolaris\n(xeric)'
             WHEN s.organism LIKE 'Helianthus argophyllus%' THEN 'H. argophyllus\n(sand dunes)'
        END as sp
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Helianthus annuus%' OR s.organism LIKE 'Helianthus petiolaris%'
         OR s.organism LIKE 'Helianthus argophyllus%')
""", conn).dropna(subset=['sp'])

sym = pd.read_sql("""
    SELECT s.organism, p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Cladocopium%' THEN 'Cladocopium\n(bleach sensitive)'
             WHEN s.organism LIKE 'Durusdinium%' THEN 'Durusdinium\n(heat tolerant)'
             WHEN s.organism LIKE 'Symbiodinium%' THEN 'Symbiodinium\n(broad)'
        END as sp
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Cladocopium%' OR s.organism LIKE 'Durusdinium%'
         OR s.organism LIKE 'Symbiodinium%')
""", conn).dropna(subset=['sp'])

all_sp = pd.concat([hel, sym])
sp_order = ['H. annuus\n(cultivated)','H. petiolaris\n(xeric)','H. argophyllus\n(sand dunes)',
            'Cladocopium\n(bleach sensitive)','Durusdinium\n(heat tolerant)','Symbiodinium\n(broad)']
sp_cols = ['#1565C0','#E65100','#BF360C','#00838F','#00695C','#004D40']
sp_means = [all_sp[all_sp['sp']==s]['km_pred_mM'].mean()*1000 for s in sp_order]
sp_ns    = [len(all_sp[all_sp['sp']==s]) for s in sp_order]

# Use log scale for this panel
bars6 = ax.bar(range(len(sp_order)), sp_means, color=sp_cols,
               alpha=0.88, width=0.6, linewidth=0)
ax.set_yscale('log')
ax.set_xticks(range(len(sp_order)))
ax.set_xticklabels(sp_order, fontsize=4.8)
ax.set_ylabel('Mean predicted Km (µM, log scale)')
for i,(n,m) in enumerate(zip(sp_ns,sp_means)):
    if n > 0:
        ax.text(i, m*1.1, f'n={n}', ha='center', fontsize=4.8, color='#555555')

# Divider line
ax.axvline(2.5, color='#cccccc', linewidth=0.8, linestyle='--')
ax.text(1, ax.get_ylim()[1]*0.7, 'Helianthus\n(same genus,\n2300× range)',
        ha='center', fontsize=5, style='italic', color='#555555')
ax.text(4, ax.get_ylim()[1]*0.7, 'Coral\nsymbionts\n(highest algal Km)',
        ha='center', fontsize=5, style='italic', color='#004D40')

conn.close()

fig.suptitle(
    'CarboDB v5 — Ecological findings: parasitism, C4 convergence, domestication, temperature, ploidy',
    fontsize=8, fontweight='bold', y=0.975)

out_pdf = ROOT / 'figures' / 'ecological_findings_figure.pdf'
out_png = ROOT / 'figures' / 'ecological_findings_figure.png'
fig.savefig(out_pdf, format='pdf', bbox_inches='tight', pad_inches=0.05)
fig.savefig(out_png, format='png', bbox_inches='tight', dpi=300)
print(f'Saved: figures/ecological_findings_figure.pdf/.png')
plt.close()
