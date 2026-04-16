"""
Script 22: TARA Oceans Metagenome Scan
- Load 5 TARA protist-fraction assemblies (MGYS00006491)
- Run CarboDB v5 binary + EC + Km models on all predicted CDS
- Summarise: how many carboxylases, which EC classes, Km distribution
- Compare Km predictions to CarboDB reference database
- Output: data/metagenome/tara_scan_results.tsv + figures/metagenome_figure.pdf

Run: python scripts/22_metagenome_scan.py
"""
import gzip, json, re, sqlite3, os
import numpy as np
import pandas as pd
import xgboost as xgb
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')

ROOT    = Path('.')
MGDIR   = ROOT / 'data/metagenome/tara_srf'
ML_DIR  = ROOT / 'data/ml'
MOD_DIR = ROOT / 'data/models'
OUT_DIR = ROOT / 'data/metagenome'

# ── 1. LOAD MODELS ────────────────────────────────────────────────────────
print("Loading models...")
bin_model = xgb.Booster(); bin_model.load_model(MOD_DIR/'binary_v5.json')
ec_model  = xgb.Booster(); ec_model.load_model(MOD_DIR/'ec_v5.json')
km_model  = xgb.Booster(); km_model.load_model(MOD_DIR/'km_v5_weighted.json')

feat_names_bin = json.load(open(ML_DIR/'feature_names_binary.json'))
feat_names_km  = json.load(open(ML_DIR/'feature_names_km.json'))
ec_label_map = json.load(open(ML_DIR/'ec_label_map_fixed.json'))
ec_classes = [ec for ec,_ in sorted(ec_label_map.items(), key=lambda x: x[1])]
print(f"  Binary features: {len(feat_names_bin)}")
print(f"  Km features: {len(feat_names_km)}")
print(f"  EC classes: {len(ec_classes)}")

# ── 2. FEATURE FUNCTIONS (from script 11) ─────────────────────────────────
AAs = list('ACDEFGHIKLMNPQRSTVWY')

def composition_features(seq):
    seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','',seq.upper())
    n   = len(seq) or 1
    # AAC
    aac = {f'aac_{a}': seq.count(a)/n for a in AAs}
    # Dipeptide (400 features — use only the ones in feat_names)
    dp = {}
    for i in range(len(seq)-1):
        k = f'dp_{seq[i]}{seq[i+1]}'
        dp[k] = dp.get(k,0) + 1
    dp = {k: v/max(n-1,1) for k,v in dp.items()}
    # Physical
    mw_map = {'A':89,'C':121,'D':133,'E':147,'F':165,'G':75,'H':155,'I':131,
               'K':146,'L':131,'M':149,'N':132,'P':115,'Q':146,'R':174,'S':105,
               'T':119,'V':117,'W':204,'Y':181}
    charge_map = {'D':-1,'E':-1,'K':1,'R':1,'H':0.1}
    phys = {
        'phys_length': n,
        'phys_mw': sum(mw_map.get(a,110) for a in seq)/n,
        'phys_charge': sum(charge_map.get(a,0) for a in seq),
        'phys_gravy': sum({'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,
                           'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,'M':1.9,'N':-3.5,
                           'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,
                           'W':-0.9,'Y':-1.3}.get(a,0) for a in seq)/n,
    }
    return {**aac, **dp, **phys}

def pfam_features_from_tsv(pfam_tsv, seq_ids, feat_names):
    """Extract Pfam hit scores for sequences from pre-computed Pfam TSV."""
    pfam_cols = [f for f in feat_names if f.startswith('pfam_')]
    pfam_domains = [f.replace('pfam_','') for f in pfam_cols]
    
    # Parse pfam TSV: seq_id, pfam_id, ...
    pfam_hits = {}  # seq_id -> {pfam_id: best_score}
    try:
        with open(pfam_tsv) as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 3: continue
                sid = parts[0]
                pid = parts[1]  # e.g. PF00016
                try: score = float(parts[2]) if len(parts)>2 else 1.0
                except: score = 1.0
                if sid not in pfam_hits: pfam_hits[sid] = {}
                if pid not in pfam_hits[sid] or pfam_hits[sid][pid] < score:
                    pfam_hits[sid][pid] = score
    except Exception as e:
        print(f"  Pfam parse warning: {e}")
    
    # Build feature matrix
    result = {}
    for sid in seq_ids:
        hits = pfam_hits.get(sid, {})
        result[sid] = {f'pfam_{d}': hits.get(d, 0.0) for d in pfam_domains}
    return result

def build_feature_vector(seq, pfam_feats, feat_names):
    comp = composition_features(seq)
    feat = {}
    for f in feat_names:
        if f in comp:
            feat[f] = comp[f]
        elif f in pfam_feats:
            feat[f] = pfam_feats[f]
        elif f.startswith('esm2_'):
            feat[f] = 0.0  # ESM-2 not available — zero fill
        elif f.startswith('ec_oh_') or f.startswith('kingdom_'):
            feat[f] = 0.0
        else:
            feat[f] = 0.0
    return np.array([feat.get(f,0.0) for f in feat_names], dtype=np.float32)

# ── 3. SCAN EACH ASSEMBLY ─────────────────────────────────────────────────
SAMPLES = {
    'MGYA00679207': 'ERZ17499708',
    'MGYA00679210': 'ERZ17294090',
    'MGYA00679214': 'ERZ17294017',
    'MGYA00679222': 'ERZ17499738',
    'MGYA00679225': 'ERZ17499655',
}

BATCH = 2000  # process in batches for memory efficiency
BIN_THRESH = 0.5

all_results = []

for mgya, erz in SAMPLES.items():
    faa_gz  = MGDIR / f'{mgya}_{erz}_FASTA_predicted_cds.faa.gz'
    pfam_f  = MGDIR / f'{mgya}_{erz}_FASTA_pfam.tsv'
    tax_f   = MGDIR / f'{mgya}_{erz}_FASTA_SSU_OTU.tsv'
    
    print(f"\n{'='*60}")
    print(f"Processing {mgya} ({erz})...")
    
    # Parse FASTA
    seqs = {}
    with gzip.open(faa_gz, 'rt') as f:
        sid, buf = None, []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if sid: seqs[sid] = ''.join(buf)
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if sid: seqs[sid] = ''.join(buf)
    
    print(f"  Sequences loaded: {len(seqs):,}")
    
    # Load Pfam for this assembly
    print(f"  Loading Pfam annotations...")
    pfam_data = pfam_features_from_tsv(pfam_f, list(seqs.keys()), feat_names_bin)
    
    # Process in batches
    seq_ids = list(seqs.keys())
    n_total = len(seq_ids)
    n_co2 = 0
    batch_results = []
    
    print(f"  Scanning {n_total:,} sequences...")
    for i in range(0, n_total, BATCH):
        batch_ids  = seq_ids[i:i+BATCH]
        batch_seqs = [seqs[s] for s in batch_ids]
        batch_pfam = [pfam_data.get(s,{}) for s in batch_ids]
        
        # Build binary feature matrix
        X_bin = np.stack([
            build_feature_vector(seq, pf, feat_names_bin)
            for seq, pf in zip(batch_seqs, batch_pfam)
        ])
        
        # Binary prediction
        dbin = xgb.DMatrix(X_bin, feature_names=feat_names_bin)
        bin_probs = bin_model.predict(dbin)
        
        # Filter CO2-interacting
        co2_mask = bin_probs >= BIN_THRESH
        if co2_mask.sum() == 0:
            continue
        
        co2_ids   = [batch_ids[j]  for j in range(len(batch_ids))  if co2_mask[j]]
        co2_seqs  = [batch_seqs[j] for j in range(len(batch_seqs)) if co2_mask[j]]
        co2_pfam  = [batch_pfam[j] for j in range(len(batch_pfam)) if co2_mask[j]]
        co2_probs = bin_probs[co2_mask]
        n_co2    += len(co2_ids)
        
        # EC prediction
        X_km = np.stack([
            build_feature_vector(seq, pf, feat_names_km)
            for seq, pf in zip(co2_seqs, co2_pfam)
        ])
        dec = xgb.DMatrix(X_km, feature_names=feat_names_km)
        ec_probs_mat = ec_model.predict(dec)
        ec_preds  = [ec_classes[np.argmax(row)] for row in ec_probs_mat]
        ec_confs  = [float(np.max(row)) for row in ec_probs_mat]
        
        # Km prediction
        km_preds_log = km_model.predict(dec)
        km_preds_mM  = [float(10**v) for v in km_preds_log]
        
        for j in range(len(co2_ids)):
            batch_results.append({
                'sample':    mgya,
                'erz':       erz,
                'seq_id':    co2_ids[j],
                'bin_prob':  float(co2_probs[j]),
                'ec_pred':   ec_preds[j],
                'ec_conf':   ec_confs[j],
                'km_pred_mM': km_preds_mM[j],
                'km_pred_log10': float(km_preds_log[j]),
                'seq_len':   len(co2_seqs[j]),
            })
        
        if (i // BATCH) % 10 == 0:
            pct = min(100, 100*(i+BATCH)/n_total)
            print(f"    {pct:.0f}% — {n_co2} CO2-interacting found so far")
    
    print(f"  DONE: {n_co2:,} / {n_total:,} predicted CO2-interacting ({100*n_co2/n_total:.1f}%)")
    all_results.extend(batch_results)

# ── 4. SAVE RESULTS ───────────────────────────────────────────────────────
df = pd.DataFrame(all_results)
out_tsv = OUT_DIR / 'tara_scan_results.tsv'
df.to_csv(out_tsv, sep='\t', index=False)
print(f"\nSaved {len(df):,} predictions → {out_tsv}")

# ── 5. SUMMARY STATS ──────────────────────────────────────────────────────
print("\n" + "="*60)
print("SUMMARY")
print("="*60)

total_seqs = sum(
    sum(1 for line in gzip.open(MGDIR/f'{m}_{e}_FASTA_predicted_cds.faa.gz','rt')
        if line.startswith('>'))
    for m,e in SAMPLES.items()
)
print(f"Total sequences scanned: {total_seqs:,}")
print(f"CO2-interacting predicted: {len(df):,} ({100*len(df)/total_seqs:.2f}%)")
print(f"\nEC class distribution:")
print(df['ec_pred'].value_counts().head(15).to_string())
print(f"\nKm statistics (all CO2-interacting, mM):")
print(df['km_pred_mM'].describe())
print(f"\nKm by EC class (µM):")
print(df.groupby('ec_pred')['km_pred_mM'].agg(
    n='count',
    mean_uM=lambda x: round(x.mean()*1000,1),
    median_uM=lambda x: round(x.median()*1000,1)
).sort_values('n', ascending=False).head(15).to_string())
print(f"\nNovel low-Km candidates (< 0.005 mM = 5 µM):")
novel = df[df['km_pred_mM'] < 0.005].sort_values('km_pred_mM')
print(f"  {len(novel):,} sequences with predicted Km < 5 µM")
print(novel[['seq_id','ec_pred','km_pred_mM','sample']].head(20).to_string(index=False))
print("\n=== DONE ===")
