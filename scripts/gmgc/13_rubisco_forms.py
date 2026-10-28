#!/usr/bin/env python3
"""
Classify predicted RuBisCO (EC 4.1.1.39) sequences into Forms I-IV based on
Pfam architecture + taxonomy + sequence length.

Form I:   PF00016 + PF02788 + PF00101 present, cyanobacteria/eukaryotic taxonomy
Form II:  PF00016 + PF02788, NO PF00101, alpha/beta proteobacteria
Form III: PF00016, archaeal taxonomy
Form IV:  PF00016 only, weak conservation, often called "RuBisCO-like protein" (RLP)

Reads:  data/gmgc/results/predictions_joined{,_full}.parquet  (chooses _full if present)
        data/gmgc/interim_full/candidates_full_pfam.tsv  (or pilot version as fallback)
Writes: data/gmgc/results/rubisco_forms.tsv
        data/gmgc/figures/13_rubisco_forms_by_taxonomy.png
"""
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path.home() / "Projects_shared/CarboDB_v3"
# Prefer full if it exists, else pilot
PRED_FULL = ROOT / "data/gmgc/results/predictions_joined_full.parquet"
PRED_PILOT = ROOT / "data/gmgc/results/predictions_joined.parquet"
PFAM_FULL = ROOT / "data/gmgc/interim_full/candidates_full_pfam.tsv"
PFAM_PILOT = ROOT / "data/gmgc/interim/candidates_pfam.tsv"
OUT_T = ROOT / "data/gmgc/results/rubisco_forms.tsv"
OUT_F = ROOT / "data/gmgc/figures/13_rubisco_forms.png"

PRED = PRED_FULL if PRED_FULL.exists() else PRED_PILOT
PFAM = PFAM_FULL if PFAM_FULL.exists() else PFAM_PILOT
print(f"Using predictions:  {PRED}", file=sys.stderr)
print(f"Using Pfam table:   {PFAM}", file=sys.stderr)

df = pd.read_parquet(PRED)
rbc = df[df['ec_pred'] == '4.1.1.39'].copy()
print(f"\nRuBisCO predictions: {len(rbc)}", file=sys.stderr)

if len(rbc) == 0:
    print("No RuBisCO predictions in dataset. Exiting.", file=sys.stderr)
    sys.exit(0)

# Load Pfam architecture per sequence
pfam_hits = defaultdict(set)
with open(PFAM) as f:
    next(f)
    for line in f:
        sid, acc, ev, bs = line.strip().split('\t')
        pfam_hits[sid].add(acc)

# Classify each RuBisCO
def classify_form(seq_id, tax_name, tax_rank, seq_len):
    pfams = pfam_hits.get(seq_id, set())
    has_large = 'PF00016' in pfams
    has_large_N = 'PF02788' in pfams
    has_small = 'PF00101' in pfams

    # Archaea -> Form III
    if tax_name and ('Archaea' in str(tax_name) or 'archaeon' in str(tax_name).lower()
                       or 'Euryarchaeota' in str(tax_name) or 'Methano' in str(tax_name)):
        return "III (archaeal)"

    # Has small subunit -> Form I (canonical Calvin cycle)
    if has_small:
        return "I (cyanobacterial/canonical)"

    # Large + N domain present, alpha/beta proteobacteria -> Form II
    if has_large and has_large_N:
        if tax_name and ('proteobacteri' in str(tax_name).lower()
                          or 'Rhodo' in str(tax_name) or 'Rhizob' in str(tax_name)):
            return "II (alpha-proteobacterial)"
        return "II or IV (multi-domain, non-cyanobacterial)"

    # Only PF00016 (no N domain, no small)
    return "IV (RuBisCO-like protein candidate)"

rbc['rubisco_form'] = rbc.apply(
    lambda r: classify_form(r['seq_id'], r.get('tax_name'),
                              r.get('tax_rank'), r.get('seq_len')),
    axis=1
)

# Summary
print("\n=== RuBisCO form distribution ===", file=sys.stderr)
print(rbc['rubisco_form'].value_counts(), file=sys.stderr)
print(f"\n=== Median predicted Km by form (mM) ===", file=sys.stderr)
print(rbc.groupby('rubisco_form')['km_pred_mM'].agg(['count', 'median', 'mean']),
      file=sys.stderr)

rbc[['seq_id','tax_name','tax_rank','seq_len','km_pred_mM','ec_conf',
     'n_carboxy_pfam','rubisco_form']].to_csv(OUT_T, sep='\t', index=False)
print(f"\nWrote {OUT_T}", file=sys.stderr)

# Plot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
form_counts = rbc['rubisco_form'].value_counts()
axes[0].barh(form_counts.index[::-1], form_counts.values[::-1], color='steelblue')
axes[0].set_xlabel('Count')
axes[0].set_title(f'RuBisCO form distribution (n={len(rbc)})')

# Km per form boxplot
forms = rbc['rubisco_form'].value_counts().index.tolist()
data = [rbc[rbc['rubisco_form']==f]['km_pred_mM'].values for f in forms]
axes[1].boxplot(data, labels=[f.split('(')[0] for f in forms], showmeans=True)
axes[1].set_yscale('log')
axes[1].set_ylabel('Predicted Km (mM, log)')
axes[1].set_title('Km distribution by RuBisCO form')
axes[1].tick_params(axis='x', rotation=30)

plt.tight_layout()
plt.savefig(OUT_F, dpi=120)
plt.close()
print(f"Saved {OUT_F}", file=sys.stderr)
