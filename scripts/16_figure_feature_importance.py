"""
CarboDB v5 — Figure 2: Feature Importance
6 panels, A4 landscape
Run: python scripts/16_figure_feature_importance.py
"""
import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from pathlib import Path
import xgboost as xgb, warnings
warnings.filterwarnings('ignore')

ROOT=Path('.'); ML_DIR=ROOT/'data/ml'; MODEL_DIR=ROOT/'data/models'
SPLIT_DIR=ROOT/'data/splits'; SHAP_DIR=ROOT/'data/shap'; FIG_DIR=ROOT/'figures'
FIG_DIR.mkdir(exist_ok=True)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,
    'axes.titlesize':8,'axes.titleweight':'bold','axes.titlepad':5,
    'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':0.6,
    'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2.5,
    'ytick.major.size':2.5,'legend.fontsize':5.8,'legend.frameon':True,
    'legend.framealpha':0.9,'legend.edgecolor':'#dddddd',
    'axes.spines.right':False,'axes.spines.top':False,'figure.dpi':300})

C={'binary':'#1B4F8A','km':'#1A6B3A','pfam':'#B03A2E','esm2':'#6B3FA0',
   'comp':'#C47900','dp':'#D97706','ec':'#888888','blast':'#7B4EA0'}
BG='#FAFAFA'

# Load SHAP data
bin_g = json.load(open(SHAP_DIR/'shap_binary_global.json'))
km_g  = json.load(open(SHAP_DIR/'shap_km_global.json'))
km_ec = json.load(open(SHAP_DIR/'shap_km_per_ec.json'))

# Load ML data for ablation curves
booster_bin=xgb.Booster(); booster_bin.load_model(str(MODEL_DIR/'binary_v5.json'))
booster_km =xgb.Booster(); booster_km.load_model(str(MODEL_DIR/'km_v5_weighted.json'))
X_bin=np.load(ML_DIR/'X_binary_test.npz')['X']
y_bin=np.load(ML_DIR/'y_binary_test.npy')
feat_bin=json.load(open(ML_DIR/'feature_names_binary.json'))
km_splits=__import__('pandas').read_csv(SPLIT_DIR/'split_km.tsv',sep='\t')
KM_EC=['4.2.1.1','4.1.1.39','4.1.1.31','4.1.1.49','6.3.4.14','4.1.1.32',
       '6.4.1.1','6.4.1.4','6.4.1.2','6.4.1.3']
mask=km_splits['ec_number'].isin(KM_EC)
km_df=km_splits[mask].reset_index(drop=True)
y_km_all=np.concatenate([np.load(ML_DIR/f'y_km_{s}_v3.npy') for s in ['train','val','test']])
X_km_all=np.vstack([np.load(ML_DIR/f'X_km_{s}_v3.npz')['X'] for s in ['train','val','test']])
_,te_idx=train_test_split(list(range(len(km_df))),test_size=0.15,random_state=42)
X_km=X_km_all[te_idx]; y_km=y_km_all[te_idx]
ec_te=km_df.iloc[te_idx]['ec_number'].values
feat_km=json.load(open(ML_DIR/'feature_names_km_v3.json'))

def fidx(fn,p): return [i for i,n in enumerate(fn) if any(n.startswith(x) for x in (p if isinstance(p,list) else [p]))]
bin_esm2=fidx(feat_bin,'esm2_'); bin_pfam=fidx(feat_bin,'pfam_')
bin_comp=fidx(feat_bin,['aac_','dp_','phys_','pse_'])
km_esm2=fidx(feat_km,'esm2_'); km_pfam=fidx(feat_km,'pfam_')
km_comp=fidx(feat_km,['aac_','dp_','phys_','pse_']); km_dp=fidx(feat_km,'dp_')
km_ec_idx=fidx(feat_km,['ec_oh_','kingdom_'])

from sklearn.metrics import roc_auc_score
def mkb(X,idx): Xv=np.zeros_like(X); Xv[:,idx]=X[:,idx]; return Xv

fig=plt.figure(figsize=(11.69,8.27))
gs=gridspec.GridSpec(2,3,figure=fig,left=0.08,right=0.97,top=0.91,bottom=0.09,
                     wspace=0.44,hspace=0.58)
axes=[fig.add_subplot(gs[r,c]) for r in range(2) for c in range(3)]

def style(ax):
    ax.set_facecolor(BG)
    ax.grid(axis='both',color='#e8e8e8',linewidth=0.4,zorder=0)
    ax.set_axisbelow(True)

def lab(ax,l,x=-0.16,y=1.10):
    ax.text(x,y,l,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top',ha='left')

# ── A: SHAP group importance — binary vs Km (diverging bars) ──────────────
ax=axes[0]; lab(ax,'A'); ax.set_title('SHAP group importance'); style(ax)
groups=['ESM-2\nembedding','Pfam\ndomains','Dipeptide\ncomp.','EC one-hot']
bin_vals=[69.26,30.74,0.0,0.0]
km_vals=[13.81,78.33,6.23,1.63]
x=np.arange(len(groups)); w=0.35
ax.bar(x-w/2,bin_vals,w,color=C['binary'],alpha=0.88,label='Binary model',linewidth=0)
ax.bar(x+w/2,km_vals, w,color=C['km'],   alpha=0.88,label='Km model',linewidth=0)
ax.set_xticks(x); ax.set_xticklabels(groups,fontsize=6.5)
ax.set_ylabel('SHAP importance (%)')
ax.set_ylim(0,85)
ax.legend(fontsize=5.5)
for i,(bv,kv) in enumerate(zip(bin_vals,km_vals)):
    if bv>2: ax.text(i-w/2,bv+1,f'{bv:.0f}%',ha='center',va='bottom',fontsize=5.5,color=C['binary'])
    if kv>2: ax.text(i+w/2,kv+1,f'{kv:.0f}%',ha='center',va='bottom',fontsize=5.5,color=C['km'])

# ── B: Top 12 individual features — binary (horizontal SHAP bar) ──────────
ax=axes[1]; lab(ax,'B'); ax.set_title('Top features — binary model'); style(ax)
top_bin=bin_g['top_global'][:12]
feat_labels=[f['feature'].replace('pfam_','').replace('esm2_','ESM2-') for f in top_bin]
shap_vals=[f['mean_abs_shap'] for f in top_bin]
feat_cols=[C['pfam'] if f['group']=='Pfam domains' else C['esm2'] for f in top_bin]
y_pos=np.arange(len(feat_labels))
ax.barh(y_pos[::-1],shap_vals,color=feat_cols,alpha=0.88,height=0.72,linewidth=0)
ax.set_yticks(y_pos); ax.set_yticklabels(feat_labels[::-1],fontsize=5.8)
ax.set_xlabel('Mean |SHAP value|')
patches=[mpatches.Patch(color=C['pfam'],label='Pfam domains',alpha=0.88),
         mpatches.Patch(color=C['esm2'],label='ESM-2 dims',alpha=0.88)]
ax.legend(handles=patches,fontsize=5.0,loc='lower right')

# ── C: Top 12 individual features — Km (horizontal SHAP bar) ─────────────
ax=axes[2]; lab(ax,'C'); ax.set_title('Top features — Km model'); style(ax)
top_km=km_g['top_global'][:12]
feat_labels_km=[f['feature'].replace('pfam_','').replace('esm2_','ESM2-').replace('dp_','dp-') for f in top_km]
shap_km_v=[f['mean_abs_shap'] for f in top_km]
km_feat_cols=[]
for f in top_km:
    g=f['group']
    km_feat_cols.append(C['pfam'] if g=='Pfam domains' else C['esm2'] if 'ESM' in g else C['dp'])
y_pos=np.arange(len(feat_labels_km))
ax.barh(y_pos[::-1],shap_km_v,color=km_feat_cols,alpha=0.88,height=0.72,linewidth=0)
ax.set_yticks(y_pos); ax.set_yticklabels(feat_labels_km[::-1],fontsize=5.8)
ax.set_xlabel('Mean |SHAP value|')
patches2=[mpatches.Patch(color=C['pfam'],label='Pfam domains',alpha=0.88),
          mpatches.Patch(color=C['esm2'],label='ESM-2 dims',alpha=0.88),
          mpatches.Patch(color=C['dp'],label='Dipeptide',alpha=0.88)]
ax.legend(handles=patches2,fontsize=5.0,loc='lower right')

# ── D: Ablation curve — binary (remove one feature group at a time) ───────
ax=axes[3]; lab(ax,'D'); ax.set_title('Feature ablation — binary AUROC'); style(ax)
abl_labels=['Full\nmodel','No\nESM-2','No\nPfam','No\nComp','ESM-2\nonly','Pfam\nonly','Comp\nonly']
abl_sets=[None, km_pfam+bin_comp, bin_esm2+bin_comp, bin_esm2+bin_pfam,
          bin_esm2, bin_pfam, bin_comp]
abl_auroc=[]
for idx in abl_sets:
    Xv=X_bin if idx is None else mkb(X_bin,idx)
    abl_auroc.append(roc_auc_score(y_bin,booster_bin.predict(xgb.DMatrix(Xv))))
abl_colors=[C['binary'],'#e74c3c','#e67e22','#888888',C['esm2'],C['pfam'],C['comp']]
x=np.arange(len(abl_labels))
bars=ax.bar(x,abl_auroc,color=abl_colors,alpha=0.88,width=0.65,linewidth=0)
ax.set_xticks(x); ax.set_xticklabels(abl_labels,fontsize=6)
ax.set_ylim(0.45,1.02); ax.set_ylabel('AUROC')
ax.axhline(1.0,color='#cccccc',linewidth=0.5,linestyle='--')
for i,v in enumerate(abl_auroc):
    ax.text(i,v+0.005,f'{v:.3f}',ha='center',va='bottom',fontsize=5.2,
            color='white' if v>0.75 else '#333333')

# ── E: Ablation — Km R² ───────────────────────────────────────────────────
ax=axes[4]; lab(ax,'E'); ax.set_title('Feature ablation — Km R²'); style(ax)
abl_km_labels=['Full\nmodel','No\nESM-2','No\nPfam','No\nComp','ESM-2\nonly',
               'Pfam\nonly','Dipep.\nonly']
abl_km_sets=[None,
             km_pfam+km_comp+km_ec_idx,
             km_esm2+km_comp+km_ec_idx,
             km_esm2+km_pfam+km_ec_idx,
             km_esm2, km_pfam, km_dp]
abl_km_r2=[]
for idx in abl_km_sets:
    Xv=X_km if idx is None else mkb(X_km,idx)
    p=booster_km.predict(xgb.DMatrix(Xv))
    abl_km_r2.append(r2_score(y_km,p))
abl_km_cols=[C['km'],'#e74c3c','#e67e22','#888888',C['esm2'],C['pfam'],C['dp']]
x=np.arange(len(abl_km_labels))
ax.bar(x,abl_km_r2,color=abl_km_cols,alpha=0.88,width=0.65,linewidth=0)
ax.set_xticks(x); ax.set_xticklabels(abl_km_labels,fontsize=6)
ax.set_ylabel('R²')
ax.axhline(0,color='#777777',linewidth=0.8)
ax.axhline(1,color='#cccccc',linewidth=0.5,linestyle='--')
for i,v in enumerate(abl_km_r2):
    ypos=max(v,0)+0.02
    ax.text(i,ypos,f'{v:.3f}',ha='center',va='bottom',fontsize=5.2,
            color='white' if v>0.5 else '#333333')

# ── F: Within-EC heatmap + per-EC top feature ─────────────────────────────
ax=axes[5]; lab(ax,'F'); ax.set_title('Within-EC Km R² by feature group')
within_methods=[
    ('Full model',None),('Pfam only',km_pfam),('ESM-2 only',km_esm2),
    ('ESM-2+Pfam',km_esm2+km_pfam),('Pfam+Comp',km_pfam+km_comp),
    ('Dipeptide',km_dp),('AAC only',fidx(feat_km,'aac_'))]
ec_within={'RuBisCO':'4.1.1.39','CA':'4.2.1.1','PEPC':'4.1.1.49','PEPCK':'4.1.1.31'}
heat=np.full((len(within_methods),len(ec_within)),np.nan)
for ri,(name,idx) in enumerate(within_methods):
    for ci,(ecname,ec) in enumerate(ec_within.items()):
        m=ec_te==ec
        if m.sum()<10: continue
        Xv=X_km[m] if idx is None else mkb(X_km[m],idx)
        p=booster_km.predict(xgb.DMatrix(Xv))
        yc=y_km[m]-y_km[m].mean(); pc=p-p.mean()
        ss=np.sum(yc**2)
        heat[ri,ci]=(1-np.sum((yc-pc)**2)/ss) if ss>0 else np.nan
im=ax.imshow(heat,cmap='YlGn',vmin=0,vmax=1.0,aspect='auto')
ax.set_xticks(range(len(ec_within))); ax.set_xticklabels(list(ec_within.keys()),fontsize=7)
ax.set_yticks(range(len(within_methods))); ax.set_yticklabels([m[0] for m in within_methods],fontsize=6)
plt.colorbar(im,ax=ax,fraction=0.038,pad=0.04,label='R²')
for ri in range(len(within_methods)):
    for ci in range(len(ec_within)):
        v=heat[ri,ci]
        if not np.isnan(v):
            fw='bold' if ri==0 else 'normal'
            ax.text(ci,ri,f'{v:.2f}',ha='center',va='center',fontsize=5.8,
                    color='#1a1a1a' if v>0.4 else '#777',fontweight=fw)

fig.suptitle('CarboDB v5 — Feature importance: SHAP analysis, ablation, and within-EC drivers',
             fontsize=8,fontweight='bold',y=0.975)
fig.savefig(FIG_DIR/'feature_importance_figure.pdf',format='pdf',bbox_inches='tight',pad_inches=0.05)
fig.savefig(FIG_DIR/'feature_importance_figure.png',format='png',bbox_inches='tight',dpi=300)
print('Saved: figures/feature_importance_figure.pdf/.png')
plt.close()
