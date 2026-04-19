"""
pipeline/predict.py
Core prediction logic — feature extraction + XGBoost inference.
Uses ModelStore singleton — never reloads models.
"""
import os
import re
import json
import time
import subprocess
import tempfile
import logging
import numpy as np
import xgboost as xgb
from typing import Optional

from ..startup import ModelStore, EC_NAMES, KM_EC_CLASSES

log = logging.getLogger(__name__)

AA_LIST = 'ACDEFGHIKLMNPQRSTVWY'
AA_SET = set(AA_LIST)

# ── Sequence cleaning ──────────────────────────────────────────────────────

def clean_sequence(seq: str) -> str:
    """Remove whitespace, headers, non-AA characters."""
    seq = re.sub(r'>.*\n', '', seq)  # remove FASTA header
    seq = ''.join(c for c in seq.upper() if c in AA_SET)
    return seq


# ── Composition features (494 dims) ───────────────────────────────────────

def compute_composition_features(seq: str) -> dict:
    """Amino acid composition + dipeptide + physicochemical features."""
    n = len(seq)
    feats = {}

    # AAC — 20 features
    for aa in AA_LIST:
        feats[f'aac_{aa}'] = seq.count(aa) / n if n > 0 else 0.0

    # Dipeptide — 400 features
    dp_counts = {}
    for i in range(len(seq) - 1):
        dp = seq[i:i+2]
        if all(c in AA_SET for c in dp):
            dp_counts[dp] = dp_counts.get(dp, 0) + 1
    total_dp = sum(dp_counts.values()) or 1
    for a1 in AA_LIST:
        for a2 in AA_LIST:
            feats[f'dpc_{a1}{a2}'] = dp_counts.get(f'{a1}{a2}', 0) / total_dp

    # Physicochemical — 74 features (length, charge, hydrophobicity, etc.)
    hydrophobic = set('VILMFYWC')
    charged_pos = set('KRH')
    charged_neg = set('DE')
    polar = set('STNQ')
    aromatic = set('FYW')
    small = set('AGSV')

    feats['length'] = n
    feats['log_length'] = np.log10(n) if n > 0 else 0
    feats['frac_hydrophobic'] = sum(1 for c in seq if c in hydrophobic) / n if n > 0 else 0
    feats['frac_charged_pos'] = sum(1 for c in seq if c in charged_pos) / n if n > 0 else 0
    feats['frac_charged_neg'] = sum(1 for c in seq if c in charged_neg) / n if n > 0 else 0
    feats['frac_polar'] = sum(1 for c in seq if c in polar) / n if n > 0 else 0
    feats['frac_aromatic'] = sum(1 for c in seq if c in aromatic) / n if n > 0 else 0
    feats['net_charge'] = feats['frac_charged_pos'] - feats['frac_charged_neg']

    return feats


# ── HMMER Pfam features (19 dims) ─────────────────────────────────────────

def compute_pfam_features(seq: str, seq_id: str = "query") -> dict:
    """Run HMMER against Pfam-A.hmm, return binary hit features."""
    pfam_hmm = os.environ.get("PFAM_HMM", "data/Pfam-A.hmm")
    feats = {f: 0.0 for f in ModelStore.pfam_features}
    pfam_hits = []

    if not os.path.exists(pfam_hmm):
        log.warning(f"Pfam-A.hmm not found at {pfam_hmm} — Pfam features will be zero")
        return feats, pfam_hits

    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as f:
            f.write(f'>{seq_id}\n{seq}\n')
            fasta_path = f.name

        with tempfile.NamedTemporaryFile(suffix='.tbl', delete=False) as f:
            tbl_path = f.name

        result = subprocess.run(
            ['hmmscan', '--domtblout', tbl_path, '--noali', '-E', '0.01',
             '--cpu', '4', pfam_hmm, fasta_path],
            capture_output=True, timeout=60
        )

        # Parse hits
        with open(tbl_path) as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.split()
                if len(parts) < 5: continue
                pfam_id = parts[1].split('.')[0]  # PF00016
                feat_name = f'pfam_{pfam_id}'
                if feat_name in feats:
                    feats[feat_name] = 1.0
                if pfam_id not in pfam_hits:
                    pfam_hits.append(pfam_id)

    except Exception as e:
        log.warning(f"HMMER failed: {e}")
    finally:
        for p in [fasta_path, tbl_path]:
            try: os.unlink(p)
            except: pass

    return feats, pfam_hits


# ── ESM-2 features (1280 dims) ────────────────────────────────────────────

def compute_esm2_features(seq: str, seq_id: str = "query") -> np.ndarray:
    """Compute ESM-2 embeddings using loaded ModelStore.esm_model."""
    if ModelStore.esm_model is None:
        log.warning("ESM-2 not loaded — returning zeros")
        return np.zeros(1280)

    try:
        import torch
        batch_converter = ModelStore.esm_alphabet.get_batch_converter()
        data = [(seq_id, seq)]
        _, _, tokens = batch_converter(data)
        tokens = tokens.to(ModelStore.esm_device)

        with torch.no_grad():
            results = ModelStore.esm_model(tokens, repr_layers=[33])
        emb = results['representations'][33][0, 1:len(seq)+1].mean(0)
        return emb.cpu().numpy()
    except Exception as e:
        log.warning(f"ESM-2 inference failed: {e}")
        return np.zeros(1280)


# ── Assemble feature vector ────────────────────────────────────────────────

def build_feature_vector(comp_feats: dict, pfam_feats: dict, esm2_emb: np.ndarray) -> np.ndarray:
    """Assemble features in exact order matching training data."""
    all_feats = {}
    all_feats.update(comp_feats)
    all_feats.update(pfam_feats)
    
    # ESM-2 features
    for i, v in enumerate(esm2_emb):
        all_feats[f'esm2_{i}'] = float(v)

    # Build vector in exact training order
    vec = np.array([all_feats.get(f, 0.0) for f in ModelStore.feature_names], dtype=np.float32)
    return vec


def build_feature_vector_no_esm2(comp_feats: dict, pfam_feats: dict) -> np.ndarray:
    """Feature vector with ESM-2 set to zeros (fast mode)."""
    return build_feature_vector(comp_feats, pfam_feats, np.zeros(1280))


# ── XGBoost inference ─────────────────────────────────────────────────────

def run_xgboost_predict(feature_vec: np.ndarray, kingdom: str = "plant") -> dict:
    """Run all three XGBoost models on a feature vector."""
    dmat = xgb.DMatrix(feature_vec.reshape(1, -1))

    # Binary prediction
    bin_prob = float(ModelStore.xgb_binary.predict(dmat)[0])
    is_carboxylase = bin_prob >= 0.5

    # EC prediction
    ec_probs = ModelStore.xgb_ec.predict(dmat)[0]
    ec_idx = int(np.argmax(ec_probs))
    ec_pred = ModelStore.ec_inv_map.get(ec_idx, 'unknown')
    ec_conf = float(ec_probs[ec_idx])
    ec_probs_dict = {ModelStore.ec_inv_map[i]: float(p)
                     for i, p in enumerate(ec_probs) if i in ModelStore.ec_inv_map}

    # Km prediction (only if carboxylase or high EC confidence)
    km_mM = None
    km_uM = None
    if (is_carboxylase or ec_conf > 0.8) and ec_pred in KM_EC_CLASSES:
        try:
            km_log = float(ModelStore.xgb_km.predict(dmat)[0])
            km_mM = float(10 ** km_log)
            km_uM = km_mM * 1000
        except Exception as e:
            log.warning(f"Km prediction failed: {e}")

    return {
        'is_carboxylase': is_carboxylase,
        'carboxylase_probability': bin_prob,
        'ec_predicted': ec_pred,
        'ec_name': EC_NAMES.get(ec_pred, ec_pred),
        'ec_confidence': ec_conf,
        'ec_probabilities': ec_probs_dict,
        'km_predicted_mM': km_mM,
        'km_predicted_uM': km_uM,
    }


# ── Novelty flag ──────────────────────────────────────────────────────────

def compute_novelty_flag(pfam_hits: list, bin_prob: float, ec_conf: float) -> str:
    """
    Simple novelty estimate without BLAST.
    Green: confirmed Pfam hit + high confidence
    Yellow: some evidence
    Red: no Pfam, low confidence — novel/uncertain
    """
    if pfam_hits and bin_prob > 0.8 and ec_conf > 0.8:
        return 'known'
    elif pfam_hits or bin_prob > 0.5:
        return 'borderline'
    else:
        return 'novel'


# ── Main predict function ─────────────────────────────────────────────────

def predict_sequence(sequence: str, mode: str = 'fast',
                     kingdom: str = 'plant', seq_id: str = 'query') -> dict:
    """
    Full prediction pipeline.
    mode: 'fast' | 'standard' | 'pfam' | 'composite'
    """
    t_start = time.time()
    seq = clean_sequence(sequence)

    if len(seq) < 10:
        raise ValueError(f"Sequence too short after cleaning: {len(seq)} aa")
    if len(seq) > 5000:
        raise ValueError(f"Sequence too long: {len(seq)} aa (max 5000)")

    # Composition features (always)
    comp_feats = compute_composition_features(seq)

    # Pfam features
    pfam_feats = {f: 0.0 for f in ModelStore.pfam_features}
    pfam_hits = []

    if mode in ('fast', 'standard', 'composite', 'pfam'):
        pfam_feats_computed, pfam_hits = compute_pfam_features(seq, seq_id)
        pfam_feats.update(pfam_feats_computed)

    if mode == 'pfam':
        # Pfam-only: no composition, no ESM-2
        vec = build_feature_vector_no_esm2({}, pfam_feats)
        result = run_xgboost_predict(vec, kingdom)

    elif mode == 'fast':
        # Composition + Pfam, no ESM-2
        vec = build_feature_vector_no_esm2(comp_feats, pfam_feats)
        result = run_xgboost_predict(vec, kingdom)

    elif mode == 'standard':
        # Full pipeline with ESM-2
        esm2_emb = compute_esm2_features(seq, seq_id)
        vec = build_feature_vector(comp_feats, pfam_feats, esm2_emb)
        result = run_xgboost_predict(vec, kingdom)

    elif mode == 'composite':
        # Run all three, return ensemble
        # Fast
        vec_fast = build_feature_vector_no_esm2(comp_feats, pfam_feats)
        r_fast = run_xgboost_predict(vec_fast, kingdom)

        # Pfam-only
        vec_pfam = build_feature_vector_no_esm2({}, pfam_feats)
        r_pfam = run_xgboost_predict(vec_pfam, kingdom)

        # Standard
        esm2_emb = compute_esm2_features(seq, seq_id)
        vec_std = build_feature_vector(comp_feats, pfam_feats, esm2_emb)
        r_std = run_xgboost_predict(vec_std, kingdom)

        # Ensemble: weighted average (standard has highest weight)
        result = r_std.copy()
        if r_std['km_predicted_uM'] and r_fast['km_predicted_uM']:
            kms = [r_fast['km_predicted_uM'], r_pfam.get('km_predicted_uM') or r_fast['km_predicted_uM'], r_std['km_predicted_uM']]
            weights = [0.2, 0.1, 0.7]
            ens_km_uM = sum(w * k for w, k in zip(weights, kms))
            result['composite_results'] = {
                'fast':  {'km_uM': r_fast['km_predicted_uM'],  'ec': r_fast['ec_predicted']},
                'pfam':  {'km_uM': r_pfam['km_predicted_uM'],  'ec': r_pfam['ec_predicted']},
                'standard': {'km_uM': r_std['km_predicted_uM'], 'ec': r_std['ec_predicted']},
                'ensemble_km_uM': round(ens_km_uM, 2),
                'agreement': 'high' if abs(r_fast['km_predicted_uM'] - r_std['km_predicted_uM']) / r_std['km_predicted_uM'] < 0.5 else 'low',
            }
    else:
        raise ValueError(f"Unknown mode: {mode}")

    # Add metadata
    result['sequence_length'] = len(seq)
    result['pfam_hits'] = pfam_hits
    result['novelty_flag'] = compute_novelty_flag(
        pfam_hits, result['carboxylase_probability'], result['ec_confidence'])
    result['features_used'] = ['composition']
    if mode != 'pfam':
        result['features_used'].append('pfam')
    if mode in ('standard', 'composite') and ModelStore.esm_model is not None:
        result['features_used'].append('esm2')
    result['mode'] = mode
    result['kingdom'] = kingdom
    result['runtime_seconds'] = round(time.time() - t_start, 2)

    return result
