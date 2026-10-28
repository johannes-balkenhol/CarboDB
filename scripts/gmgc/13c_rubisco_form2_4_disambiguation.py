#!/usr/bin/env python3
"""Disambiguate the 471 'Form II or IV' sequences."""
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd

ROOT = Path.home() / "Projects_shared/CarboDB_v3"

rbc = pd.read_csv(ROOT / "data/gmgc/results/rubisco_forms_v2.tsv", sep='\t')
ambig_mask = rbc['rubisco_form'] == 'II or IV (multi-domain, non-cyanobacterial)'
ambig = rbc[ambig_mask].copy()
print(f"Ambiguous Form II or IV sequences: {len(ambig)}", file=sys.stderr)

pfam_arch = defaultdict(list)
with open(ROOT / "data/gmgc/interim_full/candidates_full_pfam.tsv") as f:
    next(f)
    for line in f:
        sid, acc, ev, bs = line.strip().split('\t')
        pfam_arch[sid].append((acc, float(bs), float(ev)))

preds = pd.read_parquet(ROOT / "data/gmgc/results/predictions_joined_full.parquet")
seq_len_map = dict(zip(preds['seq_id'], preds['seq_len']))

def reclassify(sid, current_form):
    pfams = pfam_arch.get(sid, [])
    pfam_accs = set(p[0] for p in pfams)
    seq_len = seq_len_map.get(sid, 0)
    if 'PF02788' in pfam_accs:
        return 'II (alpha-proteobacterial)', 'PF02788 present'
    if seq_len >= 450:
        return 'II (likely, length-based)', f'len={seq_len}>=450'
    elif seq_len <= 420:
        return 'IV (likely RLP, length-based)', f'len={seq_len}<=420'
    else:
        return current_form, f'len={seq_len} in 420-450 boundary'

new_classifications = []
for _, row in ambig.iterrows():
    new_form, reason = reclassify(row['seq_id'], row['rubisco_form'])
    new_classifications.append({
        'seq_id': row['seq_id'],
        'old_form': row['rubisco_form'],
        'new_form': new_form,
        'reason': reason,
    })

reclass = pd.DataFrame(new_classifications)
print("\n=== Reclassification of 471 ambiguous sequences ===")
print(reclass['new_form'].value_counts())
print(f"\nStill ambiguous: {(reclass['new_form']==reclass['old_form']).sum()}")
print(f"Reassigned to Form II: {reclass['new_form'].str.startswith('II ').sum()}")
print(f"Reassigned to Form IV: {reclass['new_form'].str.startswith('IV ').sum()}")

rbc_v3 = rbc.copy()
reclass_map = dict(zip(reclass['seq_id'], reclass['new_form']))
rbc_v3['rubisco_form'] = rbc_v3.apply(
    lambda r: reclass_map.get(r['seq_id'], r['rubisco_form']),
    axis=1
)

print("\n=== Final Form distribution (v3) ===")
print(rbc_v3['rubisco_form'].value_counts())

rbc_v3.to_csv(ROOT / "data/gmgc/results/rubisco_forms_v3.tsv", sep='\t', index=False)
reclass.to_csv(ROOT / "data/gmgc/results/rubisco_form_disambiguation_log.tsv",
                sep='\t', index=False)

km_data = preds[preds['is_carboxylase']][['seq_id','km_pred_mM']]
rbc_v3 = rbc_v3.merge(km_data, on='seq_id', how='left')
print("\n=== Median Km per form (v3 after disambiguation) ===")
print(rbc_v3.groupby('rubisco_form')['km_pred_mM'].agg(['count','median','mean']).round(3))
