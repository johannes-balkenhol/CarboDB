#!/usr/bin/env python3
"""
Literature benchmark: compare v5 Km predictions to published
measured Km values for marine carboxylases from the same genera.

This is the honest validation: do our predictions match measurements?
"""
import pandas as pd
import sys

# Published Km values from peer-reviewed literature
# Format: (genus, EC, Km_mM, source, condition)
LIT_KM_MEASURED = [
    # RuBisCO Form I — Prochlorococcus
    ("Prochlorococcus", "4.1.1.39", 0.750, "Scott et al. 2007, MIT9313 purified pH 7.5"),
    ("Prochlorococcus", "4.1.1.39", 0.295, "Roberts et al. 2012, carboxysomes pH 8"),
    ("Prochlorococcus", "4.1.1.39", 0.263, "Hopkinson et al. 2014, MED4 crude extract pH 8"),
    # RuBisCO Form I — Synechococcus
    ("Synechococcus",   "4.1.1.39", 0.180, "Badger 1980, S. elongatus Km(CO2)"),
    # RuBisCO Form II
    ("Rhodospirillum",  "4.1.1.39", 0.090, "Form II Robinson 2003"),
    # Carbonic anhydrase
    ("Escherichia coli (b-CA)", "4.2.1.1", 8.4, "Cronk et al. (E. coli CynT)"),
    # Acetyl-CoA carboxylase
    ("Escherichia coli", "6.4.1.2", 0.030, "BRENDA E. coli ACC for acetyl-CoA"),
    # PEP carboxylase
    ("Synechococcus",   "4.1.1.31", 0.140, "BRENDA cyanobacterial PEPC Km(PEP)"),
    ("Escherichia coli", "4.1.1.31", 0.150, "BRENDA E. coli PEPC"),
]

# Load our predictions
df = pd.read_parquet('data/gmgc/results/predictions_joined_full.parquet')
carb = df[df['is_carboxylase'] & (df['ec_conf']>=0.9)].copy()

print('═══════════════════════════════════════════════════════════════')
print(' LITERATURE BENCHMARK — v5 predictions vs measured Km')
print('═══════════════════════════════════════════════════════════════\n')

print(f'{"Genus":<25} {"EC":<10} {"Lit Km":<10} {"Our median":<12} {"Our 95% range":<20} {"Ratio":<10}')
print('-' * 90)

for genus, ec, lit_km, source in LIT_KM_MEASURED:
    if genus == 'Escherichia coli (b-CA)':
        sub = carb[(carb['ec_pred']==ec) & 
                    (carb['tax_name'].str.contains('Escherichia', na=False))]
    elif genus == 'Rhodospirillum':
        sub = carb[(carb['ec_pred']==ec) &
                    (carb['tax_name'].str.contains('Rhodospir', na=False))]
    else:
        sub = carb[(carb['ec_pred']==ec) &
                    (carb['tax_name'].str.contains(genus, na=False))]
    
    if len(sub) == 0:
        print(f'{genus:<25} {ec:<10} {lit_km:<10.4f}  -- no predictions for this genus -- ')
        continue
    
    med = sub['km_pred_mM'].median()
    p5 = sub['km_pred_mM'].quantile(0.05)
    p95 = sub['km_pred_mM'].quantile(0.95)
    ratio = med / lit_km
    range_str = f'{p5:.3f}-{p95:.3f}'
    
    flag = ""
    if 0.5 <= ratio <= 2.0:
        flag = "  ✓ within 2x"
    elif 0.3 <= ratio <= 3.0:
        flag = "  ~ within 3x"
    elif ratio < 0.1 or ratio > 10:
        flag = "  ✗ off by >10x"
    else:
        flag = "  ! off by 3-10x"
    
    print(f'{genus:<25} {ec:<10} {lit_km:<10.4f} {med:<12.4f} {range_str:<20} {ratio:<10.2f}{flag}')
    print(f'  source: {source}')
    print(f'  our n={len(sub):,}')
    print()

print('═══ Overall benchmark interpretation ═══')
print('Within 2x of literature: predictions are quantitatively reliable')
print('Within 3x:                qualitatively useful, exact value imprecise')
print('Off by >10x:              prediction-vs-measurement mismatch — needs investigation')
