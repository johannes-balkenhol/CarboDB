"""
CarboDB v5 — Publication benchmark figure
6 panels, 3×2, A4 format (297×210mm landscape)
Saves to: figures/benchmark_figure.pdf + .png
"""
import numpy as np, json, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, roc_auc_score, r2_score
from sklearn.model_selection import train_test_split
from scipy.stats import pearsonr
from pathlib import Path
import xgboost as xgb, warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT     = Path('.')
ML_DIR   = ROOT / 'data/ml'
MODEL_DIR= ROOT / 'data/models'
SPLIT_DIR= ROOT / 'data/splits'
BENCH_DIR= ROOT / 'data/benchmark'
FIG_DIR  = ROOT / 'figures'
FIG_DIR.mkdir(exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'Arial',
    'font.size':        7,
    'axes.labelsize':   7,
    'axes.titlesize':   8,
    'axes.titleweight': 'bold',
    'xtick.labelsize':  6,
    'ytick.labelsize':  6,
    'axes.linewidth':   0.6,
    'xtick.major.width':0.5,
    'ytick.major.width':0.5,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'legend.fontsize':  6,
    'legend.frameon':   True,
    'legend.framealpha':0.85,
    'legend.edgecolor': '#cccccc',
    'legend.handlelength':1.5,
    'axes.spines.right':False,
    'axes.spines.top':  False,
    'figure.dpi':       300,
})

COLORS = {
    'carbodb':  '#1B4F8A',
    'esm2':     '#6B3FA0',
    'esm2comp': '#1A7A4A',
    'pfamcomp': '#C47900',
    'pfam':     '#B03A2E',
    'comp':     '#5D6D7E',
    'blast':    '#884EA0',
    'unikp':    '#E74C3C',
    'ecmean':   '#95A5A6',
}

PANEL_LABELS = list('ABCDEF')

# ── Load models & data ─────────────────────────────────────────────────────
print('Loading models...')
booster_bin = xgb.Booster()
booster_bin.load_model(str(MODEL_DIR/'binary_v5.json'))
booster_km  = xgb.Booster()
booster_km.load_model(str(MODEL_DIR/'km_v5_weighted.json'))
booster_ec  = xgb.Booster()
booster_ec.load_model(str(MODEL_DIR/'ec_v5.json'))

X_bin = np.load(ML_DIR/'X_binary_test.npz')['X']
y_bin = np.load(ML_DIR/'y_binary_test.npy')
feat_bin = json.load(open(ML_DIR/'feature_names_binary.json'))

km_splits = pd.read_csv(SPLIT_DIR/'split_km.tsv', sep='\t')
KM_EC = ['4.2.1.1','4.1.1.39','4.1.1.31','4.1.1.49','6.3.4.14',
         '4.1.1.32','6.4.1.1','6.4.1.4','6.4.1.2','6.4.1.3']
mask   = km_splits['ec_number'].isin(KM_EC)
km_df  = km_splits[mask].reset_index(drop=True)
y_km_all = np.concatenate([np.load(ML_DIR/(f'y_km_{s}_v3.npy')) for s in ['train','val','test']])
X_km_all = np.vstack([np.load(ML_DIR/(f'X_km_{s}_v3.npz'))['X'] for s in ['train','val','test']])
_,te_idx = train_test_split(list(range(len(km_df))), test_size=0.15, random_state=42)
X_km = X_km_all[te_idx]; y_km = y_km_all[te_idx]
ec_te = km_df.iloc[te_idx]['ec_number'].values
feat_km = json.load(open(ML_DIR/'feature_names_km_v3.json'))

X_ec = np.load(ML_DIR/'X_ec_test_fixed.npz')['X']
y_ec = np.load(ML_DIR/'y_ec_test_fixed.npy')
feat_ec = json.load(open(ML_DIR/'feature_names_ec.json'))
label_map = json.load(open(ML_DIR/'ec_label_map_fixed.json'))

def fidx(fn,p): return [i for i,n in enumerate(fn) if any(n.startswith(x) for x in (p if isinstance(p,list) else [p]))]
bin_esm2=fidx(feat_bin,'esm2_'); bin_pfam=fidx(feat_bin,'pfam_')
bin_comp=fidx(feat_bin,['aac_','dp_','phys_','pse_'])
km_esm2=fidx(feat_km,'esm2_'); km_pfam=fidx(feat_km,'pfam_')
km_comp=fidx(feat_km,['aac_','dp_','phys_','pse_']); km_ec_idx=fidx(feat_km,['ec_oh_','kingdom_'])
ec_esm2=fidx(feat_ec,'esm2_'); ec_pfam=fidx(feat_ec,'pfam_')
ec_comp=fidx(feat_ec,['aac_','dp_','phys_','pse_'])

def mkb(X,idx): Xv=np.zeros_like(X); Xv[:,idx]=X[:,idx]; return Xv

# BLAST Km results
blast_km = pd.read_csv(BENCH_DIR/'blast/km_blast_results.tsv', sep='\t', header=None,
    names=['qid','sid','pident','length','evalue','bitscore'])
blast_km['qcdb'] = blast_km['qid'].str.split('|').str[0]
blast_km['km_hit'] = blast_km['sid'].str.extract(r'km=(-?[\d.]+)').astype(float)
blast_km_lookup = dict(zip(blast_km.qcdb, blast_km.km_hit))
km_cdb_ids = km_df.iloc[te_idx]['cdb_id'].values
blast_km_pred = np.array([blast_km_lookup.get(c, np.nan) for c in km_cdb_ids])

# ── Figure layout ──────────────────────────────────────────────────────────
# A4 landscape = 297mm x 210mm
fig = plt.figure(figsize=(11.69, 8.27))
gs  = gridspec.GridSpec(2, 3, figure=fig,
                        left=0.07, right=0.97,
                        top=0.93,  bottom=0.10,
                        wspace=0.38, hspace=0.52)
axes = [fig.add_subplot(gs[r,c]) for r in range(2) for c in range(3)]

def label_panel(ax, letter, x=-0.18, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top', ha='left')

# ═══════════════════════════════════════════════════════════════════════════
# A — ROC curves (zoomed FPR 0–0.05 inset)
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[0]
label_panel(ax, 'A')
ax.set_title('Binary carboxylase detection')

variants_bin = [
    ('CarboDB v5',       None,                 COLORS['carbodb'],  '-',  2.0),
    ('ESM-2 only',       bin_esm2,             COLORS['esm2'],     '-',  1.4),
    ('ESM-2+Comp',       bin_esm2+bin_comp,    COLORS['esm2comp'], '--', 1.2),
    ('Pfam+Comp',        bin_pfam+bin_comp,    COLORS['pfamcomp'], '--', 1.2),
    ('Pfam only',        bin_pfam,             COLORS['pfam'],     ':',  1.1),
]

legend_lines = []
for name, idx, col, ls, lw in variants_bin:
    Xv = X_bin if idx is None else mkb(X_bin, idx)
    scores = booster_bin.predict(xgb.DMatrix(Xv))
    fpr, tpr, _ = roc_curve(y_bin, scores)
    auroc = roc_auc_score(y_bin, scores)
    ax.plot(fpr, tpr, color=col, linestyle=ls, linewidth=lw, alpha=0.9)
    legend_lines.append(Line2D([0],[0], color=col, linestyle=ls, linewidth=lw,
                               label=f'{name} ({auroc:.4f})'))

# BLAST
blast_bin = pd.read_csv(BENCH_DIR/'blast/blast_results.tsv', sep='\t', header=None,
    names=['qid','sid','pident','length','evalue','bitscore'])
blast_bin['qcdb']   = blast_bin['qid'].str.split('|').str[0]
blast_bin['slabel'] = blast_bin['sid'].str.extract(r'label=(\d)').astype(float)
bin_split = pd.read_csv(SPLIT_DIR/'split_binary.tsv', sep='\t').query('split=="test"').reset_index(drop=True).iloc[:len(X_bin)]
bscores = np.zeros(len(bin_split))
for i,row in bin_split.iterrows():
    if i >= len(bscores): break
    h = blast_bin[blast_bin['qcdb']==row['cdb_id']]
    if len(h)>0 and h.iloc[0]['slabel']==1:
        bscores[i] = h.iloc[0]['pident']/100
blast_auroc = roc_auc_score(y_bin[:len(bscores)], bscores)
fpr_b, tpr_b, _ = roc_curve(y_bin[:len(bscores)], bscores)
ax.plot(fpr_b, tpr_b, color=COLORS['blast'], linestyle=':', linewidth=1.0)
legend_lines.append(Line2D([0],[0], color=COLORS['blast'], linestyle=':', linewidth=1.0,
                           label=f'BLAST-NN ({blast_auroc:.4f})'))

ax.plot([0,1],[0,1],'--', color='#aaaaaa', linewidth=0.7, alpha=0.6)
ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
ax.set_xlim(-0.01,1.01); ax.set_ylim(-0.01,1.02)
ax.legend(handles=legend_lines, loc='lower right', fontsize=5.2)

# Zoom inset FPR 0–0.05
axins = ax.inset_axes([0.08, 0.3, 0.45, 0.45])
for name, idx, col, ls, lw in variants_bin:
    Xv = X_bin if idx is None else mkb(X_bin, idx)
    scores = booster_bin.predict(xgb.DMatrix(Xv))
    fpr, tpr, _ = roc_curve(y_bin, scores)
    axins.plot(fpr, tpr, color=col, linestyle=ls, linewidth=lw*0.9)
axins.set_xlim(0, 0.005); axins.set_ylim(0.985, 1.001)
axins.tick_params(labelsize=4.5)
axins.set_xlabel('FPR', fontsize=4.5); axins.set_ylabel('TPR', fontsize=4.5)
axins.set_title('FPR 0–0.5%', fontsize=4.8, fontweight='normal')
axins.spines['right'].set_visible(False); axins.spines['top'].set_visible(False)
ax.indicate_inset_zoom(axins, edgecolor='#666666', linewidth=0.5)

# ═══════════════════════════════════════════════════════════════════════════
# B — Fragment robustness
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[1]
label_panel(ax, 'B')
ax.set_title('Fragment robustness')

frac_labels = ['Full', '75%', '50%', '25%']
def trunc(X, frac):
    rng = np.random.default_rng(42); Xt = X.copy()
    for i in bin_esm2: Xt[:,i] *= frac
    for i in bin_comp: Xt[:,i] *= frac
    for i in bin_pfam:
        m = rng.random(len(X)) < (1-frac); Xt[m,i] = 0
    return Xt

frac_variants = [
    ('CarboDB v5', None,          COLORS['carbodb']),
    ('ESM-2',      bin_esm2,      COLORS['esm2']),
    ('ESM-2+Comp', bin_esm2+bin_comp, COLORS['esm2comp']),
    ('Pfam+Comp',  bin_pfam+bin_comp, COLORS['pfamcomp']),
]
fracs = [1.0, 0.75, 0.50, 0.25]
x = np.arange(len(frac_labels))
w = 0.19
offsets = np.linspace(-(len(frac_variants)-1)*w/2, (len(frac_variants)-1)*w/2, len(frac_variants))

for j,(name,idx,col) in enumerate(frac_variants):
    aurocs = []
    for frac in fracs:
        Xc = trunc(X_bin, frac)
        Xv = Xc if idx is None else mkb(Xc, idx)
        aurocs.append(roc_auc_score(y_bin, booster_bin.predict(xgb.DMatrix(Xv))))
    ax.bar(x + offsets[j], aurocs, w*0.9, color=col, alpha=0.88, label=name)

ax.set_xticks(x); ax.set_xticklabels(frac_labels)
ax.set_ylim(0.55, 1.005)
ax.set_ylabel('AUROC')
ax.set_xlabel('Sequence completeness')
ax.legend(loc='lower left', fontsize=5.2)
ax.axhline(1.0, color='#cccccc', linewidth=0.5, linestyle='--')

# ═══════════════════════════════════════════════════════════════════════════
# C — EC class prediction accuracy
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[2]
label_panel(ax, 'C')
ax.set_title('EC class prediction (top-1)')

ec_ablation = json.load(open(BENCH_DIR/'ec_ablation.json'))
abl_labels  = ['CarboDB v5', 'ESM-2 only', 'ESM-2+Comp', 'Pfam+Comp', 'Pfam only', 'PANTHER*']
abl_vals    = [
    ec_ablation['ablation']['full'] * 100,
    ec_ablation['ablation']['esm2'] * 100,
    ec_ablation['ablation']['esm2_comp'] * 100,
    ec_ablation['ablation']['pfam_comp'] * 100,
    ec_ablation['ablation']['pfam'] * 100,
    94.0
]
abl_cols = [COLORS['carbodb'], COLORS['esm2'], COLORS['esm2comp'],
            COLORS['pfamcomp'], COLORS['pfam'], '#888888']

bars = ax.barh(range(len(abl_labels)), abl_vals, color=abl_cols, alpha=0.88, height=0.62)
ax.set_yticks(range(len(abl_labels))); ax.set_yticklabels(abl_labels, fontsize=6)
ax.set_xlim(30, 103)
ax.set_xlabel('Top-1 accuracy (%)')
ax.axvline(100, color='#cccccc', linewidth=0.5, linestyle='--')
for i,v in enumerate(abl_vals):
    ax.text(v+0.5, i, f'{v:.1f}%', va='center', fontsize=5.5,
            color='#333333' if v < 80 else COLORS['carbodb'] if i==0 else '#444444')
ax.text(0.98, 0, '*published', transform=ax.transAxes,
        fontsize=5, color='#888888', ha='right', va='bottom')

# ═══════════════════════════════════════════════════════════════════════════
# D — Km R² comparison
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[3]
label_panel(ax, 'D')
ax.set_title('CO₂ affinity (Km) prediction')

km_variants_d = [
    ('CarboDB v5',       None,                     COLORS['carbodb']),
    ('Pfam+Comp',        km_pfam+km_comp,           COLORS['pfamcomp']),
    ('Pfam only',        km_pfam,                   COLORS['pfam']),
    ('ESM-2 only',       km_esm2,                   COLORS['esm2']),
    ('EC-class mean',    None,                      COLORS['ecmean']),   # special
    ('BLAST-NN',         None,                      COLORS['blast']),    # special
    ('UniKP-style*',     None,                      '#E74C3C'),
]

km_r2s = []
for name, idx, col in km_variants_d:
    if name == 'EC-class mean':
        ec_m = {e: y_km[ec_te==e].mean() for e in KM_EC if (ec_te==e).sum()>0}
        pb = np.array([ec_m.get(e, y_km.mean()) for e in ec_te])
        r2 = r2_score(y_km, pb)
    elif name == 'BLAST-NN':
        mbl = ~np.isnan(blast_km_pred)
        r2 = r2_score(y_km[mbl], blast_km_pred[mbl])
    elif name == 'UniKP-style*':
        r2 = -0.036
    else:
        Xv = X_km if idx is None else mkb(X_km, idx)
        pred = booster_km.predict(xgb.DMatrix(Xv))
        r2 = r2_score(y_km, pred)
    km_r2s.append(r2)

km_labels = [v[0] for v in km_variants_d]
km_colors = [v[2] for v in km_variants_d]
y_pos = np.arange(len(km_labels))

bars = ax.barh(y_pos, km_r2s, color=km_colors, alpha=0.88, height=0.62)
ax.set_yticks(y_pos); ax.set_yticklabels(km_labels, fontsize=6)
ax.set_xlabel('R²  (log₁₀ Km)')
ax.set_xlim(-0.75, 1.08)
ax.axvline(0, color='#888888', linewidth=0.7, linestyle='-')
ax.axvline(1, color='#cccccc', linewidth=0.5, linestyle='--')
for i,v in enumerate(km_r2s):
    xpos = v + 0.02 if v >= 0 else v - 0.02
    ha = 'left' if v >= 0 else 'right'
    col = km_colors[i] if i == 0 else '#444444'
    ax.text(xpos, i, f'{v:.3f}', va='center', ha=ha, fontsize=5.5, color=col)
ax.text(0.98, 0, '*general model on CO₂', transform=ax.transAxes,
        fontsize=5, color='#888888', ha='right', va='bottom')

# ═══════════════════════════════════════════════════════════════════════════
# E — Within-EC Km heatmap
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[4]
label_panel(ax, 'E')
ax.set_title('Within-EC Km prediction (R²)')

within_methods = [
    ('Full model',   None),
    ('ESM-2 only',   km_esm2),
    ('Pfam only',    km_pfam),
    ('ESM-2+Pfam',   km_esm2+km_pfam),
    ('Pfam+Comp',    km_pfam+km_comp),
    ('No ESM-2',     km_pfam+km_comp+km_ec_idx),
    ('Dipeptide',    fidx(feat_km,'dp_')),
]
ec_within_classes = {
    'RuBisCO':'4.1.1.39', 'CA':'4.2.1.1',
    'PEPC':'4.1.1.49',     'PEPCK':'4.1.1.31'
}

heat = np.zeros((len(within_methods), len(ec_within_classes)))
for ri,(name,idx) in enumerate(within_methods):
    for ci,(ecname,ec) in enumerate(ec_within_classes.items()):
        m = ec_te == ec
        if m.sum() < 10: heat[ri,ci] = np.nan; continue
        Xv = X_km[m] if idx is None else mkb(X_km[m], idx)
        p  = booster_km.predict(xgb.DMatrix(Xv))
        yc = y_km[m] - y_km[m].mean()
        pc = p - p.mean()
        ss = np.sum(yc**2)
        heat[ri,ci] = (1 - np.sum((yc-pc)**2)/ss) if ss>0 else np.nan

im = ax.imshow(heat, cmap='YlGn', vmin=0, vmax=1.0, aspect='auto')
ax.set_xticks(range(len(ec_within_classes)))
ax.set_xticklabels(list(ec_within_classes.keys()), fontsize=6.5)
ax.set_yticks(range(len(within_methods)))
ax.set_yticklabels([m[0] for m in within_methods], fontsize=6)
plt.colorbar(im, ax=ax, fraction=0.035, pad=0.04, label='R²')
for ri in range(len(within_methods)):
    for ci in range(len(ec_within_classes)):
        v = heat[ri,ci]
        if not np.isnan(v):
            ax.text(ci, ri, f'{v:.2f}', ha='center', va='center',
                    fontsize=5.8, color='#1a1a1a' if v>0.5 else '#555555',
                    fontweight='bold' if ri==0 else 'normal')

# ═══════════════════════════════════════════════════════════════════════════
# F — Predicted vs experimental Km scatter
# ═══════════════════════════════════════════════════════════════════════════
ax = axes[5]
label_panel(ax, 'F')
ax.set_title('Predicted vs experimental Km (R²=0.953)')

scatter_data = json.load(open(BENCH_DIR/'km_scatter_data.json'))
ec_color_map = {
    '4.1.1.39': COLORS['carbodb'],
    '4.2.1.1':  COLORS['esm2comp'],
    '4.1.1.49': COLORS['pfamcomp'],
    '4.1.1.31': COLORS['pfam'],
    '4.1.1.32': COLORS['pfam'],
}
ec_label_map2 = {
    '4.1.1.39':'RuBisCO','4.2.1.1':'Carbonic anhydrase',
    '4.1.1.49':'PEPC','4.1.1.31':'PEPCK','4.1.1.32':'PEPCK (GTP)',
}

plotted_ec = set()
for pt in scatter_data:
    ec = pt['ec']
    col = ec_color_map.get(ec, '#888888')
    ax.scatter(pt['exp'], pt['pred'], c=col, s=10, alpha=0.62,
               edgecolors='none', zorder=3)
    plotted_ec.add(ec)

lims = (-3.3, 2.3)
ax.plot(lims, lims, '--', color='#aaaaaa', linewidth=0.8, zorder=2, label='y=x')
ax.fill_between(lims, [l-np.log10(2) for l in lims],
                [l+np.log10(2) for l in lims],
                color='#cccccc', alpha=0.2, zorder=1, label='±2-fold')
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel('Experimental log₁₀(Km / mM)')
ax.set_ylabel('Predicted log₁₀(Km / mM)')

legend_patches = [mpatches.Patch(color=ec_color_map.get(ec,'#888888'),
                                 label=ec_label_map2.get(ec,ec))
                  for ec in ['4.1.1.39','4.2.1.1','4.1.1.49','4.1.1.31']]
legend_patches += [Line2D([0],[0],color='#aaaaaa',linestyle='--',linewidth=0.8,label='y=x'),
                   mpatches.Patch(color='#cccccc',alpha=0.5,label='±2-fold')]
ax.legend(handles=legend_patches, loc='upper left', fontsize=5.2, ncol=1)

# ── Final touches ──────────────────────────────────────────────────────────
fig.suptitle(
    'CarboDB v5 — Comprehensive benchmark across binary classification, EC prediction, and Km regression',
    fontsize=8.5, fontweight='bold', y=0.975, color='#1a1a1a'
)

out_pdf = FIG_DIR / 'benchmark_figure.pdf'
out_png = FIG_DIR / 'benchmark_figure.png'
fig.savefig(out_pdf, dpi=300, bbox_inches='tight', format='pdf')
fig.savefig(out_png, dpi=300, bbox_inches='tight', format='png')
print(f'Saved: {out_pdf}')
print(f'Saved: {out_png}')
plt.close()
