"""
CarboDB v5 — Figure 3: Biological Findings
6 panels, A4 landscape
Run: python scripts/17_figure_biological.py
"""
import numpy as np, json, pandas as pd, sqlite3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')

ROOT=Path('.'); BENCH_DIR=ROOT/'data/benchmark'; BIO_DIR=ROOT/'data/biological'
FIG_DIR=ROOT/'figures'; FIG_DIR.mkdir(exist_ok=True)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,
    'axes.titlesize':8,'axes.titleweight':'bold','axes.titlepad':5,
    'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':0.6,
    'xtick.major.width':0.5,'ytick.major.width':0.5,'xtick.major.size':2.5,
    'ytick.major.size':2.5,'legend.fontsize':5.8,'legend.frameon':True,
    'legend.framealpha':0.9,'legend.edgecolor':'#dddddd',
    'axes.spines.right':False,'axes.spines.top':False,'figure.dpi':300})

BG='#FAFAFA'
EC_COLS={'4.1.1.39':'#1B4F8A','4.2.1.1':'#1A7A4A','4.1.1.49':'#C47900',
         '4.1.1.31':'#B03A2E','4.1.1.32':'#7B2D8B','6.3.4.14':'#5D6D7E',
         '6.4.1.1':'#884EA0','6.3.5.5':'#666666'}
EC_NAMES={'4.1.1.39':'RuBisCO','4.2.1.1':'Carb. anhydrase','4.1.1.49':'PEPC',
          '4.1.1.31':'PEPCK','4.1.1.32':'PEPCK-GTP','6.3.4.14':'Pyr. carbox.',
          '6.4.1.1':'ACC','6.3.5.5':'CPS'}

fig=plt.figure(figsize=(11.69,8.27))
gs=gridspec.GridSpec(2,3,figure=fig,left=0.08,right=0.97,top=0.91,bottom=0.10,
                     wspace=0.46,hspace=0.58)
axes=[fig.add_subplot(gs[r,c]) for r in range(2) for c in range(3)]

def style(ax):
    ax.set_facecolor(BG)
    ax.grid(axis='both',color='#e8e8e8',linewidth=0.4,zorder=0)
    ax.set_axisbelow(True)

def lab(ax,l,x=-0.16,y=1.10):
    ax.text(x,y,l,transform=ax.transAxes,fontsize=10,fontweight='bold',va='top',ha='left')

# ── A: Predicted Km distribution by EC class (violin) ─────────────────────
ax=axes[0]; lab(ax,'A'); ax.set_title('Predicted Km by EC class')
violin_df=pd.read_csv(BENCH_DIR/'km_violin_data.csv')
ecs=['4.1.1.39','4.2.1.1','4.1.1.49','4.1.1.31','4.1.1.32','6.3.4.14','6.4.1.1','6.3.5.5']
ecs=[e for e in ecs if e in violin_df['ec_number'].values]
data_v=[violin_df[violin_df['ec_number']==ec]['km_pred_log10'].values for ec in ecs]
vp=ax.violinplot(data_v,positions=range(len(ecs)),widths=0.7,
                  showmedians=True,showextrema=False)
for i,(body,ec) in enumerate(zip(vp['bodies'],ecs)):
    body.set_facecolor(EC_COLS.get(ec,'#888888'))
    body.set_alpha(0.75); body.set_edgecolor('none')
vp['cmedians'].set_color('#222222'); vp['cmedians'].set_linewidth(1.2)
ax.set_xticks(range(len(ecs)))
ax.set_xticklabels([EC_NAMES.get(e,e) for e in ecs],rotation=35,ha='right',fontsize=5.8)
ax.set_ylabel('Predicted log₁₀(Km / mM)')
ax.axhline(0,color='#aaaaaa',linewidth=0.6,linestyle='--')
style(ax)

# ── B: Database composition (stacked bar) ────────────────────────────────
ax=axes[1]; lab(ax,'B'); ax.set_title('Database composition')
style(ax)
categories=['CO₂-interacting\n(label=1)','Evol. related\n(label=2)','Non-CO₂\n(label=0)']
counts=[503841,121512,1755093]
colors=['#1B4F8A','#C47900','#888888']
bars=ax.bar(categories,counts,color=colors,alpha=0.88,width=0.55,linewidth=0)
ax.set_ylabel('Number of sequences')
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,_:f'{x/1e6:.1f}M'))
for bar,n in zip(bars,counts):
    pct=n/sum(counts)*100
    ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+8000,
            f'{n:,}\n({pct:.1f}%)',ha='center',va='bottom',fontsize=6,fontweight='bold')
ax.set_ylim(0,max(counts)*1.22)

# ── C: Novel low-Km candidates dot plot ───────────────────────────────────
ax=axes[2]; lab(ax,'C'); ax.set_title('Novel low-Km candidates (predicted)')
style(ax)
novel=json.load(open(BIO_DIR/'novel_km_candidates.json'))
low5=novel['low_km_top5']
high5=novel['high_km_top5']

labels_low=[f"{c['organism'].split()[0][:12]}\n{c['organism'].split()[-1][:8]}"
            for c in low5]
vals_low=[np.log10(c['km_pred_mM']) for c in low5]
labels_high=[f"{c['organism'].split()[0][:12]}\n{c['organism'].split()[-1][:8]}"
             for c in high5]
vals_high=[np.log10(c['km_pred_mM']) for c in high5]

y_low=np.arange(len(low5))
y_high=np.arange(len(low5),len(low5)+len(high5))

ax.scatter(vals_low,y_low,c='#1B4F8A',s=55,zorder=4,alpha=0.9)
ax.scatter(vals_high,y_high,c='#B03A2E',s=55,zorder=4,alpha=0.9,marker='D')
ax.set_yticks(list(y_low)+list(y_high))
ax.set_yticklabels(labels_low+labels_high,fontsize=5.5)
ax.set_xlabel('Predicted log₁₀(Km / mM)')
ax.axvline(-1.5,color='#cccccc',linewidth=0.5,linestyle=':')
patches=[mpatches.Patch(color='#1B4F8A',label='Low Km (novel)'),
         mpatches.Patch(color='#B03A2E',label='High Km (novel)',alpha=0.8)]
ax.legend(handles=patches,fontsize=5.2,loc='lower right')
for v,y in zip(vals_low,y_low):
    ax.text(v-0.04,y,f'{10**v*1000:.2f} µM',ha='right',va='center',fontsize=4.8,color='#1B4F8A')

# ── D: Carnivorous plants predicted Km ────────────────────────────────────
ax=axes[3]; lab(ax,'D'); ax.set_title('Carnivorous plants — predicted Km')
carn=pd.read_csv(BENCH_DIR/'carnivorous_km.csv')
genera_order=['Sarracenia','Drosera','Dionaea','Utricularia','Pinguicula','Nepenthes','Cephalotus']
genera_order=[g for g in genera_order if g in carn['genus'].values]
carn_colors={'Sarracenia':'#2196F3','Drosera':'#4CAF50','Dionaea':'#F44336',
             'Utricularia':'#9C27B0','Pinguicula':'#FF9800','Nepenthes':'#00BCD4',
             'Cephalotus':'#795548'}

# Strip/jitter plot
np.random.seed(42)
for i,genus in enumerate(genera_order):
    sub=carn[carn['genus']==genus]['km_pred_log10'].values
    jitter=np.random.normal(0,0.08,len(sub))
    col=carn_colors.get(genus,'#888888')
    ax.scatter(sub,np.full(len(sub),i)+jitter,c=col,s=8,alpha=0.55,
               edgecolors='none',zorder=3,rasterized=True)
    med=np.median(sub)
    ax.plot([med,med],[i-0.35,i+0.35],color=col,linewidth=2.0,zorder=4)

# Reference lines: typical C3 (−1.52), C4 (−2.7)
ax.axvline(-1.52,color='#1B4F8A',linewidth=0.8,linestyle='--',alpha=0.6,label='Typical C3 RuBisCO')
ax.axvline(-2.70,color='#1A7A4A',linewidth=0.8,linestyle=':',alpha=0.6,label='Typical C4 RuBisCO')
ax.set_yticks(range(len(genera_order))); ax.set_yticklabels(genera_order)
ax.set_xlabel('Predicted log₁₀(Km / mM)')
ax.legend(fontsize=4.8,loc='lower right')
style(ax)
ax.text(0.98,0.02,'median shown as line',transform=ax.transAxes,
        fontsize=4.8,color='#888888',ha='right',va='bottom')

# ── E: Ecological patterns — low-Km groups dot/strip ─────────────────────
ax=axes[4]; lab(ax,'E'); ax.set_title('Ecological Km patterns')
style(ax)

eco_groups={
    'C4 Chloridoideae\ngrasses':('4.1.1.39',['Cynodon','Chloris','Zoysia','Bouteloua','Dactyloctenium'],'#2E7D32'),
    'Galdieria/Cyanidiales\n(hot springs)':('4.1.1.39',['Galdieria','Cyanidium'],'#BF360C'),
    'Marine deep-reef\nmacroalgae':('4.1.1.39',['Lobophora','Desmarestia','Ishige'],'#1565C0'),
    'Helicobacter pylori\n(gastric)':('4.2.1.1',['Helicobacter'],'#6A1B9A'),
    'Sus scrofa\n(mammalian CA)':('4.2.1.1',['Sus scrofa'],'#E65100'),
}

conn=sqlite3.connect('data/primary/carbodb.sqlite')
np.random.seed(99)
y_pos=0
yticks=[]; ytick_labels=[]
for label,(ec,genera,col) in eco_groups.items():
    org_filter=" OR ".join([f"s.organism LIKE '{g}%'" for g in genera])
    df=pd.read_sql(f'''SELECT p.km_pred_log10 FROM predictions p
        JOIN sequences s ON s.id=p.sequence_id
        WHERE s.ec_number='{ec}' AND ({org_filter})
        AND p.km_pred_log10 IS NOT NULL LIMIT 500''', conn)
    if len(df)==0: continue
    vals=df['km_pred_log10'].values
    jitter=np.random.normal(0,0.07,len(vals))
    ax.scatter(vals,np.full(len(vals),y_pos)+jitter,c=col,s=8,alpha=0.55,
               edgecolors='none',zorder=3,rasterized=True)
    med=np.median(vals)
    ax.plot([med,med],[y_pos-0.38,y_pos+0.38],color=col,linewidth=2.2,zorder=4)
    ax.text(med+0.05,y_pos,f'{10**med:.3f} mM',va='center',fontsize=5.2,color=col)
    yticks.append(y_pos); ytick_labels.append(label); y_pos+=1

conn.close()
ax.set_yticks(yticks); ax.set_yticklabels(ytick_labels,fontsize=5.8)
ax.set_xlabel('Predicted log₁₀(Km / mM)')
ax.axvline(-1.52,color='#aaaaaa',linewidth=0.6,linestyle='--')

# ── F: Km prediction validation — experimental vs predicted ───────────────
ax=axes[5]; lab(ax,'F'); ax.set_title('Model validation — 7 known proteins')
style(ax)

val_data=[
    ('P00875','Spinach RuBisCO','4.1.1.39',-2.00,-1.979,'✓'),
    ('P00880','R. rubrum RuBisCO\n(Form II)','4.1.1.39',0.00,-1.774,'✗'),
    ('P00918','Human CA-II','4.2.1.1',0.903,0.117,'✗'),
    ('P00864','E. coli PEPC','4.1.1.49',-0.046,-0.021,'~'),
    ('P11498','Human pyr. carbox.','6.4.1.1',-0.398,0.277,'~'),
    ('P05165','Human propionyl-CoA C.','6.4.1.3',-0.301,0.572,'~'),
    ('P0ABD5','E. coli acetyl-CoA C.','6.4.1.2',-0.699,-0.291,'✓'),
]
exp_v=[v[3] for v in val_data]
pred_v=[v[4] for v in val_data]
names=[v[1] for v in val_data]
marks=[v[5] for v in val_data]
ec_c=[EC_COLS.get(v[2],'#888888') for v in val_data]
mark_col={'✓':'#1A7A4A','~':'#C47900','✗':'#B03A2E'}

lims=(-2.5,1.3)
ax.plot(lims,lims,'--',color='#aaaaaa',linewidth=0.9,zorder=2)
ax.fill_between(lims,[l-np.log10(2) for l in lims],[l+np.log10(2) for l in lims],
                color='#cccccc',alpha=0.22,zorder=1)
for i,(e,p,n,m,col) in enumerate(zip(exp_v,pred_v,names,marks,ec_c)):
    ax.scatter(e,p,c=col,s=60,zorder=4,edgecolors='white',linewidth=0.5)
    ax.text(e+0.04,p,n.split('\n')[0],fontsize=4.8,va='center',color='#333333')

ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel('Experimental log₁₀(Km / mM)')
ax.set_ylabel('Predicted log₁₀(Km / mM)')
patches=[mpatches.Patch(color=mark_col['✓'],label='Within 2-fold ✓'),
         mpatches.Patch(color=mark_col['~'],label='Within 5-fold ~'),
         mpatches.Patch(color=mark_col['✗'],label='>5-fold error ✗')]
ax.legend(handles=patches,fontsize=5.0,loc='upper left')
ax.text(0.03,0.97,'n=7 validation proteins',transform=ax.transAxes,
        fontsize=5.5,va='top',color='#555555')

fig.suptitle('CarboDB v5 — Biological findings: database landscape, novel candidates, and ecological patterns',
             fontsize=8,fontweight='bold',y=0.975)
fig.savefig(FIG_DIR/'biological_findings_figure.pdf',format='pdf',bbox_inches='tight',pad_inches=0.05)
fig.savefig(FIG_DIR/'biological_findings_figure.png',format='png',bbox_inches='tight',dpi=300)
print('Saved: figures/biological_findings_figure.pdf/.png')
plt.close()
