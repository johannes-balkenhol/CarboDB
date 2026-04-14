"""
Script 19: Taxonomic Deep Dive
- Kingdom-level Km distributions for all major EC classes
- Top genera/families with extreme Km (high and low)
- Convergent low-Km lineages across kingdoms
- Supplementary taxonomic tables for Paper 1

Output:
  data/taxonomy/kingdom_km_summary.json
  data/taxonomy/top_genera_low_km.csv
  data/taxonomy/top_genera_high_km.csv
  data/taxonomy/convergent_low_km.csv
  data/taxonomy/supplementary_taxonomic_table.csv
  figures/taxonomy_figure.pdf/.png

Run: python scripts/19_taxonomic_deepdive.py
"""
import sqlite3, pandas as pd, numpy as np, json, re
from pathlib import Path
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings; warnings.filterwarnings('ignore')

ROOT = Path('.')
CONN = sqlite3.connect(ROOT / 'data/primary/carbodb.sqlite')
(ROOT / 'data/taxonomy').mkdir(exist_ok=True)
(ROOT / 'figures').mkdir(exist_ok=True)

EC_NAMES = {
    '4.1.1.39': 'RuBisCO', '4.2.1.1': 'Carbonic anhydrase',
    '4.1.1.49': 'PEPC',    '4.1.1.31': 'PEPCK',
    '6.3.4.14': 'Pyr. carboxylase', '6.4.1.1': 'ACC',
    '4.1.1.32': 'PEPCK-GTP',
}

# ── Kingdom assignment from organism string ────────────────────────────────
KINGDOM_PATTERNS = {
    'Bacteria': [
        'bacteria','Escherichia','Bacillus','Streptomyces','Pseudomonas',
        'Synechococcus','Cyanobacterium','Rhodobacter','Chlorobium',
        'Thermodesulfovibrio','Aquifex','Helicobacter','Clostridium',
        'Staphylococcus','Mycobacterium','Salmonella','Vibrio',
        'Prochlorococcus','Nitrosomonas','Thiobacillus','Acidithiobacillus',
        'Magnetococcus','Rhodospirillum','Thiomicrospira','Halothiobacillus',
    ],
    'Archaea': [
        'archaea','Methan','Sulfolobus','Thermococcus','Archaeoglobus',
        'Haloferax','Pyrococcus','Thermoplasma','Crenarchaeota',
        'Euryarchaeota','Nanoarchaeota','Thaumarchaeota',
    ],
    'Plants': [
        'Arabidopsis','Oryza','Zea','Vitis','Solanum','Glycine','Nicotiana',
        'Spinacia','Beta','Triticum','Hordeum','Sorghum','Populus','Medicago',
        'Physcomitrella','Marchantia','Selaginella','Pinus','Picea',
        'Chlamydomonas','Chlorella','Volvox','Ostreococcus','Micromonas',
        'Galdieria','Cyanidium','Porphyra','Gracilaria',
        'Dionaea','Drosera','Nepenthes','Utricularia','Sarracenia','Pinguicula',
        'Cynodon','Chloris','Zoysia','Sorghum','Panicum',
    ],
    'Animals': [
        'Homo','Mus','Rattus','Sus','Bos','Gallus','Danio','Xenopus',
        'Drosophila','Caenorhabditis','Strongylocentrotus','Ciona',
        'Equus','Ovis','Canis','Felis','Macaca',
    ],
    'Fungi': [
        'Saccharomyces','Aspergillus','Candida','Neurospora','Fusarium',
        'Cryptococcus','Schizosaccharomyces','Yarrowia','Pichia',
    ],
}

def assign_kingdom(organism):
    if not organism: return 'Unknown'
    for kingdom, patterns in KINGDOM_PATTERNS.items():
        for p in patterns:
            if p.lower() in organism.lower():
                return kingdom
    return 'Unknown'

def extract_genus(organism):
    if not organism: return 'Unknown'
    parts = organism.strip().split()
    return parts[0] if parts else 'Unknown'

print("Loading predictions + sequences...")
df = pd.read_sql("""
    SELECT s.id, s.uniprot_id, s.ec_number, s.organism, s.source,
           s.reviewed, s.length,
           p.km_pred_mM, p.km_pred_log10, p.co2_prob, p.ec_pred
    FROM sequences s
    JOIN predictions p ON p.sequence_id = s.id
    WHERE s.label=1 AND p.km_pred_mM IS NOT NULL
    AND s.ec_number IN ('4.1.1.39','4.2.1.1','4.1.1.49','4.1.1.31',
                        '6.3.4.14','6.4.1.1','4.1.1.32','6.3.5.5',
                        '6.3.4.16','6.4.1.3','6.4.1.2')
""", CONN)

print(f"Loaded {len(df):,} sequences with Km predictions")

df['kingdom'] = df['organism'].apply(assign_kingdom)
df['genus']   = df['organism'].apply(extract_genus)
df['km_uM']   = df['km_pred_mM'] * 1000

# ── 1. Kingdom-level Km summary ───────────────────────────────────────────
print("\n=== KINGDOM-LEVEL Km SUMMARY ===")
kingdom_summary = df.groupby(['ec_number','kingdom']).agg(
    n=('km_pred_log10','count'),
    mean_log10=('km_pred_log10','mean'),
    std_log10=('km_pred_log10','std'),
    median_log10=('km_pred_log10','median'),
    mean_mM=('km_pred_mM','mean'),
    median_mM=('km_pred_mM','median'),
).reset_index()
kingdom_summary['mean_uM'] = kingdom_summary['mean_mM'] * 1000
print(kingdom_summary[kingdom_summary['ec_number']=='4.1.1.39'].sort_values('mean_log10').to_string(index=False))
kingdom_summary.to_csv('data/taxonomy/kingdom_km_summary.csv', index=False)
json.dump(kingdom_summary.to_dict(orient='records'),
          open('data/taxonomy/kingdom_km_summary.json','w'), indent=2, default=str)

# Kingdom ANOVA for RuBisCO
rubisco = df[df['ec_number']=='4.1.1.39']
kingdoms_rubisco = [rubisco[rubisco['kingdom']==k]['km_pred_log10'].values
                    for k in ['Plants','Bacteria','Archaea','Animals','Fungi']
                    if len(rubisco[rubisco['kingdom']==k]) >= 5]
f, p = stats.f_oneway(*kingdoms_rubisco)
print(f"\nRuBisCO kingdom ANOVA: F={f:.2f} p={p:.2e}")

# ── 2. Top genera by low Km ───────────────────────────────────────────────
print("\n=== TOP GENERA LOW Km (RuBisCO) ===")
rubisco_genera = rubisco.groupby('genus').agg(
    n=('km_pred_mM','count'),
    mean_km_mM=('km_pred_mM','mean'),
    median_km_mM=('km_pred_mM','median'),
    std_km_mM=('km_pred_mM','std'),
    kingdom=('kingdom','first'),
    example_organism=('organism','first'),
).reset_index()
rubisco_genera = rubisco_genera[rubisco_genera['n'] >= 3]
rubisco_genera['mean_km_uM'] = rubisco_genera['mean_km_mM'] * 1000

low_km = rubisco_genera.nsmallest(30, 'mean_km_mM')
high_km = rubisco_genera.nlargest(20, 'mean_km_mM')
print("Top 15 lowest Km genera (RuBisCO):")
print(low_km[['genus','kingdom','n','mean_km_uM','example_organism']].head(15).to_string(index=False))
print("\nTop 10 highest Km genera (RuBisCO):")
print(high_km[['genus','kingdom','n','mean_km_uM','example_organism']].head(10).to_string(index=False))
low_km.to_csv('data/taxonomy/top_genera_low_km.csv', index=False)
high_km.to_csv('data/taxonomy/top_genera_high_km.csv', index=False)

# ── 3. Convergent low-Km analysis ─────────────────────────────────────────
# Identify genera with mean Km < 0.005 mM across ALL EC classes
print("\n=== CONVERGENT LOW Km ACROSS EC CLASSES ===")
convergent = df.groupby(['ec_number','genus','kingdom']).agg(
    n=('km_pred_mM','count'),
    mean_km_mM=('km_pred_mM','mean'),
).reset_index()
convergent = convergent[(convergent['mean_km_mM'] < 0.005) & (convergent['n'] >= 2)]
print(convergent.sort_values('mean_km_mM').head(30).to_string(index=False))
convergent.to_csv('data/taxonomy/convergent_low_km.csv', index=False)

# ── 4. Full supplementary table ───────────────────────────────────────────
supp = df.groupby(['ec_number','genus','kingdom']).agg(
    n=('km_pred_mM','count'),
    mean_km_mM=('km_pred_mM','mean'),
    median_km_mM=('km_pred_mM','median'),
    std_km_mM=('km_pred_mM','std'),
    min_km_mM=('km_pred_mM','min'),
    max_km_mM=('km_pred_mM','max'),
    n_reviewed=('reviewed','sum'),
).reset_index()
supp = supp[supp['n'] >= 2].sort_values(['ec_number','mean_km_mM'])
supp.to_csv('data/taxonomy/supplementary_taxonomic_table.csv', index=False)
print(f"\nSupplementary table: {len(supp):,} genus×EC combinations")

# ── 5. Figure ─────────────────────────────────────────────────────────────
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.linewidth':0.6,
    'axes.spines.right':False,'axes.spines.top':False,'figure.dpi':300,
    'axes.titlesize':8,'axes.titleweight':'bold'})

KCOLS = {'Plants':'#2E7D32','Bacteria':'#1565C0','Archaea':'#BF360C',
          'Animals':'#6A1B9A','Fungi':'#E65100','Unknown':'#9E9E9E'}

fig = plt.figure(figsize=(11.69, 8.27))
gs  = gridspec.GridSpec(2, 3, figure=fig, left=0.08, right=0.97,
                        top=0.91, bottom=0.10, wspace=0.44, hspace=0.58)
axes = [fig.add_subplot(gs[r,c]) for r in range(2) for c in range(3)]

BG = '#FAFAFA'
def style(ax): ax.set_facecolor(BG); ax.grid(axis='both',color='#e8e8e8',linewidth=0.4,zorder=0); ax.set_axisbelow(True)
def lab(ax,l): ax.text(-0.16,1.10,l,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top',ha='left')

# Panel A: RuBisCO Km by kingdom violin
ax = axes[0]; lab(ax,'A'); ax.set_title('RuBisCO Km by kingdom'); style(ax)
k_order = ['Plants','Bacteria','Archaea','Fungi','Animals']
k_data  = [rubisco[rubisco['kingdom']==k]['km_pred_log10'].values for k in k_order]
k_data  = [d for d in k_data if len(d) >= 5]
k_labels= [k for k,d in zip(k_order,[rubisco[rubisco['kingdom']==k]['km_pred_log10'].values for k in k_order]) if len(d)>=5]
vp = ax.violinplot(k_data, positions=range(len(k_labels)), widths=0.7,
                   showmedians=True, showextrema=False)
for body,k in zip(vp['bodies'],k_labels):
    body.set_facecolor(KCOLS.get(k,'#888888')); body.set_alpha(0.75); body.set_edgecolor('none')
vp['cmedians'].set_color('#222222'); vp['cmedians'].set_linewidth(1.2)
ax.set_xticks(range(len(k_labels))); ax.set_xticklabels(k_labels, fontsize=6, rotation=25, ha='right')
ax.set_ylabel('Predicted log₁₀(Km / mM)')
for i,(k,d) in enumerate(zip(k_labels,k_data)):
    ax.text(i, ax.get_ylim()[0]+0.1 if ax.get_ylim()[0]<0 else 0.05,
            f'n={len(d):,}', ha='center', fontsize=5, color='#555555')

# Panel B: Kingdom Km means across all EC classes (heatmap)
ax = axes[1]; lab(ax,'B'); ax.set_title('Mean predicted Km (log₁₀ mM) by kingdom × EC')
ec_plot = ['4.1.1.39','4.2.1.1','4.1.1.49','4.1.1.31','6.3.4.14','6.4.1.1']
k_plot  = ['Plants','Bacteria','Archaea','Fungi','Animals']
heat = np.full((len(k_plot),len(ec_plot)), np.nan)
for ri,k in enumerate(k_plot):
    for ci,ec in enumerate(ec_plot):
        sub = df[(df['kingdom']==k)&(df['ec_number']==ec)]
        if len(sub)>=3: heat[ri,ci] = sub['km_pred_log10'].mean()
im = ax.imshow(heat, cmap='RdYlGn_r', aspect='auto', vmin=-2.5, vmax=1.5)
ax.set_xticks(range(len(ec_plot)))
ax.set_xticklabels([EC_NAMES.get(e,e) for e in ec_plot], rotation=30, ha='right', fontsize=6)
ax.set_yticks(range(len(k_plot))); ax.set_yticklabels(k_plot, fontsize=6.5)
plt.colorbar(im, ax=ax, fraction=0.038, pad=0.04, label='log₁₀(Km/mM)')
for ri in range(len(k_plot)):
    for ci in range(len(ec_plot)):
        v = heat[ri,ci]
        if not np.isnan(v):
            ax.text(ci,ri,f'{v:.1f}', ha='center', va='center', fontsize=5.5,
                    color='white' if abs(v)>1.5 else '#222222')

# Panel C: Top 15 lowest-Km genera (all EC classes, dot plot)
ax = axes[2]; lab(ax,'C'); ax.set_title('Top 15 lowest-Km genera (RuBisCO)'); style(ax)
top15 = low_km.head(15)
y = np.arange(len(top15))
cols = [KCOLS.get(k,'#888888') for k in top15['kingdom']]
ax.scatter(top15['mean_km_uM'], y[::-1], c=cols, s=50, zorder=4, alpha=0.88)
ax.set_yticks(y)
ax.set_yticklabels([f"{r['genus']} ({r['kingdom'][:3]})" for _,r in top15.iloc[::-1].iterrows()], fontsize=5.8)
ax.set_xlabel('Mean predicted Km (µM)')
from matplotlib.patches import Patch
patches = [Patch(color=c,label=k) for k,c in KCOLS.items() if k!='Unknown']
ax.legend(handles=patches, fontsize=4.5, loc='lower right')

# Panel D: Km gradient across kingdoms for RuBisCO only (strip plot)
ax = axes[3]; lab(ax,'D'); ax.set_title('RuBisCO Km gradient across kingdoms'); style(ax)
np.random.seed(42)
yp = 0; yticks=[]; ylabels=[]
for k in ['Plants','Fungi','Archaea','Bacteria']:
    sub = rubisco[rubisco['kingdom']==k]['km_pred_log10'].values
    if len(sub) < 5: continue
    jit = np.random.normal(0, 0.06, len(sub))
    ax.scatter(sub, np.full(len(sub),yp)+jit, c=KCOLS[k], s=4, alpha=0.3,
               edgecolors='none', zorder=3, rasterized=True)
    ax.plot([np.median(sub)]*2, [yp-0.38,yp+0.38], color=KCOLS[k], linewidth=2, zorder=4)
    ax.text(np.median(sub)+0.08, yp, f'{10**np.median(sub)*1000:.2f} µM',
            va='center', fontsize=5.5, color=KCOLS[k])
    yticks.append(yp); ylabels.append(f'{k}\n(n={len(sub):,})'); yp+=1
ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=6)
ax.set_xlabel('Predicted log₁₀(Km / mM)')
ax.axvline(-1.52, color='#aaaaaa', linewidth=0.6, linestyle='--', alpha=0.7)

# Panel E: Within-kingdom Km variance decomposition (horizontal bars)
ax = axes[4]; lab(ax,'E'); ax.set_title('Km variance: between vs within kingdom'); style(ax)
ec_var = []
for ec in ['4.1.1.39','4.2.1.1','4.1.1.49','4.1.1.31']:
    sub = df[df['ec_number']==ec]
    total_var = sub['km_pred_log10'].var()
    group_means = sub.groupby('kingdom')['km_pred_log10'].mean()
    group_ns    = sub.groupby('kingdom')['km_pred_log10'].count()
    between_var = np.average((group_means - sub['km_pred_log10'].mean())**2, weights=group_ns)
    within_var  = total_var - between_var
    ec_var.append({'ec': EC_NAMES[ec], 'between': between_var/total_var*100,
                   'within': within_var/total_var*100})
ec_var_df = pd.DataFrame(ec_var)
x = np.arange(len(ec_var_df))
ax.bar(x, ec_var_df['between'], color='#1B4F8A', alpha=0.85, label='Between-kingdom', width=0.55)
ax.bar(x, ec_var_df['within'],  bottom=ec_var_df['between'], color='#C47900', alpha=0.85,
       label='Within-kingdom', width=0.55)
ax.set_xticks(x); ax.set_xticklabels(ec_var_df['ec'], fontsize=6.5)
ax.set_ylabel('% of total Km variance'); ax.set_ylim(0,105)
ax.legend(fontsize=5.5)
for i,row in ec_var_df.iterrows():
    ax.text(i, row['between']/2, f"{row['between']:.0f}%", ha='center',
            va='center', fontsize=5.5, color='white', fontweight='bold')

# Panel F: Novel top genera/families summary table as text
ax = axes[5]; lab(ax,'F'); ax.set_title('Top convergent low-Km lineages')
ax.axis('off')
col_labels = ['Genus','Kingdom','EC','n','Mean Km (µM)']
table_data = []
for _,row in convergent.sort_values('mean_km_mM').head(12).iterrows():
    table_data.append([row['genus'], row['kingdom'][:4], row['ec_number'],
                       str(int(row['n'])), f"{row['mean_km_mM']*1000:.3f}"])
t = ax.table(cellText=table_data, colLabels=col_labels,
             loc='center', cellLoc='center')
t.auto_set_font_size(False); t.set_fontsize(6)
t.scale(1, 1.3)
for (r,c), cell in t.get_celld().items():
    if r == 0: cell.set_facecolor('#1B4F8A'); cell.set_text_props(color='white', fontweight='bold')
    elif r % 2 == 0: cell.set_facecolor('#F5F5F5')
    cell.set_edgecolor('#CCCCCC')

fig.suptitle('CarboDB v5 — Taxonomic deep dive: Km distributions across kingdoms and lineages',
             fontsize=8, fontweight='bold', y=0.975)
fig.savefig('figures/taxonomy_figure.pdf', format='pdf', bbox_inches='tight')
fig.savefig('figures/taxonomy_figure.png', format='png', bbox_inches='tight', dpi=300)
print('\nSaved: figures/taxonomy_figure.pdf/.png')
plt.close()

print("\n=== DONE ===")
print("Output files:")
for f in sorted(Path('data/taxonomy').glob('*')): print(f" {f}")
print(" figures/taxonomy_figure.pdf/.png")

CONN.close()
