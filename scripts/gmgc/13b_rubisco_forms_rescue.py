#!/usr/bin/env python3
"""
RuBisCO Form classifier v2 — Pfam-architecture-based, NOT gated on v5 EC routing.

Rationale: script 13 gated on (ec_pred == '4.1.1.39'), but v5 routes many small-
subunit sequences (PF00101) to non-RuBisCO EC classes because PF00101 alone
lacks catalytic activity. This script ignores v5's EC routing and classifies
any sequence with carboxylase Pfam architecture matching a RuBisCO form.

Classification rules (Pfam architecture only):
  Form I  : has PF00101 (RuBisCO_small)  → cyanobacterial / canonical Calvin
  Form II : has PF00016+PF02788 (large+N), NO PF00101, alpha/beta-proteobacteria
  Form III: has PF00016, archaeal taxonomy
  Form IV : has PF00016 only (no small, no PF02788) → RLP candidate

Reads:  data/gmgc/results/predictions_joined_full.parquet
        data/gmgc/interim_full/candidates_full_pfam.tsv
Writes: data/gmgc/results/rubisco_forms_v2.tsv
        data/gmgc/figures/13b_rubisco_forms_v2.png
"""
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path.home() / "Projects_shared/CarboDB_v3"
PRED = ROOT / "data/gmgc/results/predictions_joined_full.parquet"
PFAM = ROOT / "data/gmgc/interim_full/candidates_full_pfam.tsv"
OUT_T = ROOT / "data/gmgc/results/rubisco_forms_v2.tsv"
OUT_F = ROOT / "data/gmgc/figures/13b_rubisco_forms_v2.png"

print(f"Loading predictions...", file=sys.stderr)
df = pd.read_parquet(PRED)
df_carb = df[df['is_carboxylase']].copy()
print(f"  Total carboxylase predictions: {len(df_carb):,}", file=sys.stderr)

print(f"\nLoading Pfam architecture from {PFAM}", file=sys.stderr)
pfam_arch = defaultdict(set)
with open(PFAM) as f:
    next(f)
    for line in f:
        sid, acc, ev, bs = line.strip().split('\t')
        pfam_arch[sid].add(acc)
print(f"  Sequences with Pfam hits: {len(pfam_arch):,}", file=sys.stderr)

# Find ALL sequences with any RuBisCO-family Pfam
RUBISCO_PFAMS = {'PF00016', 'PF02788', 'PF00101'}
candidates = set()
for sid, hits in pfam_arch.items():
    if hits & RUBISCO_PFAMS:
        candidates.add(sid)
print(f"  Sequences with any RuBisCO-family Pfam: {len(candidates):,}", file=sys.stderr)

# Restrict to predictions joined to taxonomy
df_rbc = df_carb[df_carb['seq_id'].isin(candidates)].copy()
print(f"  Sequences also in joined predictions: {len(df_rbc):,}", file=sys.stderr)

def classify(seq_id, tax_name):
    p = pfam_arch.get(seq_id, set())
    has_large  = 'PF00016' in p
    has_largeN = 'PF02788' in p
    has_small  = 'PF00101' in p
    tn = str(tax_name) if tax_name else ''

    # Form I: has the small subunit (regardless of taxonomy — cyanos, plants, some chemolithos)
    if has_small:
        return "I (cyanobacterial/canonical)"

    # Archaeal taxonomy → Form III
    if any(k in tn for k in ['Archaea','archaeon','Methano','Thermo',
                              'Euryarchaeota','Crenarchaeota','Thaumarchaeota',
                              'Nitrosopelagicus','Nitrosopumilus']):
        return "III (archaeal)"

    # Has large + N domain → Form II if alpha/beta proteobacteria, otherwise II/IV ambiguous
    if has_large and has_largeN:
        if any(k in tn for k in ['proteobacteri','Rhodo','Rhizob','alpha','beta',
                                   'Bradyrhizob','Methylobacter']):
            return "II (alpha-proteobacterial)"
        return "II or IV (multi-domain, non-cyanobacterial)"

    # Only PF00016 → Form IV (RLP)
    if has_large:
        return "IV (RuBisCO-like protein candidate)"

    # Only small subunit, no large — orphan small subunit (probably Form I fragment)
    if has_small:
        return "I (orphan small subunit)"

    return "unclassified"

df_rbc['rubisco_form'] = df_rbc.apply(
    lambda r: classify(r['seq_id'], r.get('tax_name')), axis=1
)

print(f"\n=== RuBisCO Form distribution (v2, Pfam-only) ===", file=sys.stderr)
print(df_rbc['rubisco_form'].value_counts(), file=sys.stderr)

print(f"\n=== Median Km (mM) per form ===", file=sys.stderr)
print(df_rbc.groupby('rubisco_form')['km_pred_mM'].agg(['count','median','mean']),
      file=sys.stderr)

# Form I breakdown by taxonomy
form_i = df_rbc[df_rbc['rubisco_form'].str.startswith('I')]
if len(form_i) > 0:
    print(f"\n=== Form I taxonomic breakdown (top 15) ===", file=sys.stderr)
    print(form_i[form_i['tax_rank']=='genus']['tax_name'].value_counts().head(15),
          file=sys.stderr)

# Write output
df_rbc[['seq_id','tax_name','tax_rank','seq_len','km_pred_mM','ec_pred',
        'ec_conf','n_carboxy_pfam','rubisco_form']].to_csv(OUT_T, sep='\t', index=False)
print(f"\nWrote {OUT_T}", file=sys.stderr)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
form_counts = df_rbc['rubisco_form'].value_counts()
colors = plt.cm.Set2.colors[:len(form_counts)]
axes[0].barh(form_counts.index[::-1], form_counts.values[::-1], color=colors)
axes[0].set_xlabel('Sequence count')
axes[0].set_title(f'RuBisCO form distribution (Pfam-architecture classification, n={len(df_rbc):,})')

forms = df_rbc['rubisco_form'].value_counts().index.tolist()
data = [df_rbc[df_rbc['rubisco_form']==f]['km_pred_mM'].values for f in forms]
bp = axes[1].boxplot(data, tick_labels=[f.split('(')[0].strip() for f in forms],
                       showmeans=True, patch_artist=True)
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
axes[1].set_yscale('log')
axes[1].set_ylabel('Predicted Km (mM, log-scale)')
axes[1].set_title('Km distribution per RuBisCO form')
axes[1].tick_params(axis='x', rotation=20)
plt.tight_layout()
plt.savefig(OUT_F, dpi=120, bbox_inches='tight')
print(f"Saved {OUT_F}", file=sys.stderr)
