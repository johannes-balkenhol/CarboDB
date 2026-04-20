#!/usr/bin/env python3
"""
11_annotate_sequence.py
=======================
CarboDB — Step 11: Annotate a new protein sequence.

Pipeline:
  1. Parse & validate FASTA input
  2. Compute amino acid composition + physicochemical features (489 dims)
  3. Run HMMER against Pfam-A.hmm  → 19 Pfam features
  4. Run HMMER against InterPro representative HMMs → 5 InterPro features
  5. Compute ESM-2 embedding (esm2_t33_650M_UR50D) → 1280 dims
  6. XGBoost predict: binary (carboxylase?), EC class, Km (if carboxylase)
  7. Output JSON to stdout (or --out FILE)

Usage:
  python scripts/11_annotate_sequence.py --fasta my_protein.fasta
  python scripts/11_annotate_sequence.py --fasta my_protein.fasta --out result.json
  python scripts/11_annotate_sequence.py --fasta my_protein.fasta --no-esm2  # fast mode
  python scripts/11_annotate_sequence.py --sequence MAKT...   # inline sequence

Output JSON schema:
  {
    "cdb_query_id": "query_001",
    "sequence_length": 475,
    "is_carboxylase": true,
    "carboxylase_probability": 0.9987,
    "confidence": "high",
    "ec_predicted": "4.1.1.39",
    "ec_name": "ribulose-bisphosphate carboxylase",
    "ec_probabilities": {"4.1.1.39": 0.987, ...},
    "km_predicted_mM": 0.18,
    "km_predicted_log10": -0.74,
    "km_ec_used": "4.1.1.39",
    "features_used": ["composition", "pfam", "interpro", "esm2"],
    "pfam_hits": ["PF00016", "PF02788"],
    "top_pfam_by_shap": [...],
    "warnings": [],
    "runtime_seconds": 12.4
  }
"""

import argparse
import json
import logging
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, ROOT, setup_logging

log = setup_logging("11_annotate_sequence")

# ── Constants ──────────────────────────────────────────────────────────────────

MODEL_DIR  = ROOT / "data" / "models"
ML_DIR     = ROOT / "data" / "ml"

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# EC class → human-readable name
EC_NAMES = {
    "4.1.1.39":  "ribulose-bisphosphate carboxylase (RuBisCO)",
    "4.2.1.1":   "carbonic anhydrase",
    "6.3.4.16":  "acetyl-CoA carboxylase (biotin carboxylase subunit)",
    "6.3.4.14":  "pyruvate carboxylase",
    "6.3.5.5":   "carbamoyl-phosphate synthase",
    "6.3.4.18":  "3-methylcrotonyl-CoA carboxylase",
    "4.1.1.49":  "phosphoenolpyruvate carboxylase",
    "6.3.3.3":   "dethiobiotin synthase",
    "4.1.1.31":  "phosphoenolpyruvate carboxykinase",
    "4.1.1.112": "2-oxoglutarate carboxylase",
    "4.1.1.32":  "phosphoenolpyruvate carboxykinase (GTP)",
    "6.4.1.1":   "pyruvate carboxylase",
    "6.4.1.2":   "acetyl-CoA carboxylase",
    "6.4.1.3":   "propionyl-CoA carboxylase",
    "6.4.1.4":   "3-methylcrotonyl-CoA carboxylase",
    "4.1.1.38":  "phosphoenolpyruvate carboxykinase (pyrophosphate)",
}

# Trainable EC classes for Km (from script 09)
KM_TRAINABLE_EC = [
    "4.2.1.1", "4.1.1.39", "4.1.1.31", "4.1.1.49",
    "6.3.4.14", "4.1.1.32", "6.4.1.1", "6.4.1.4",
    "6.4.1.2", "6.4.1.3",
]

# Pfam IDs in feature order (must match training)
CARBOXY_PFAM_LIST = sorted(CFG.CARBOXY_PFAM)

# InterPro feature columns (must match training)
INTERPRO_COLS = ["n_pfam_hits", "n_panther_hits", "n_tigrfam_hits",
                 "n_cath_hits", "n_superfamily_hits"]


# ══════════════════════════════════════════════════════════════════════════════
# 1. FASTA parsing & validation
# ══════════════════════════════════════════════════════════════════════════════

def parse_fasta(fasta_path: Path = None, sequence: str = None):
    """
    Returns list of (seq_id, sequence) tuples.
    Accepts either a FASTA file or a raw sequence string.
    """
    records = []

    if sequence:
        seq = sequence.upper().strip()
        records.append(("query_001", seq))

    elif fasta_path:
        with open(fasta_path) as f:
            cur_id, cur_seq = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if cur_id:
                        records.append((cur_id, "".join(cur_seq)))
                    cur_id = line[1:].split()[0]
                    cur_seq = []
                elif line:
                    cur_seq.append(line.upper())
            if cur_id:
                records.append((cur_id, "".join(cur_seq)))
    else:
        raise ValueError("Provide either --fasta or --sequence")

    return records


def validate_sequence(seq_id: str, seq: str) -> list:
    """Returns list of warning strings; empty = valid."""
    warnings = []
    if len(seq) < CFG.SEQ_MIN_LEN:
        warnings.append(f"Sequence too short ({len(seq)} < {CFG.SEQ_MIN_LEN})")
    if len(seq) > CFG.SEQ_MAX_LEN:
        warnings.append(f"Sequence very long ({len(seq)} aa) — truncating ESM-2 to 1022")
    invalid = set(seq) - CFG.SEQ_VALID_AA
    if invalid:
        warnings.append(f"Non-standard amino acids removed: {sorted(invalid)}")
    n_unique = len(set(seq) & CFG.SEQ_VALID_AA)
    if n_unique < CFG.SEQ_MIN_UNIQUE:
        warnings.append(f"Low sequence complexity ({n_unique} unique AA)")
    return warnings


# ══════════════════════════════════════════════════════════════════════════════
# 2. Composition features (489 dims)
# ══════════════════════════════════════════════════════════════════════════════

def compute_composition(seq: str) -> dict:
    """
    Compute the same 489-dim composition vector used in training:
      - 20 AAC (amino acid composition)
      - 400 dipeptide composition
      - 20 pseudo-AAC (λ=1 correlation)
      - 20 physicochemical (hydrophobicity, charge, MW proxies)
      - 7 structural / global
      - 5 catalytic core motifs (inv_ features)
      - 17 EC motif hits (motif_ features) — set to 0 for unknown query

    Returns dict: feature_name → float
    """
    seq = "".join(c for c in seq.upper() if c in CFG.SEQ_VALID_AA)
    n = len(seq)
    if n == 0:
        return {}

    feats = {}

    # ── AAC ───────────────────────────────────────────────────────────────────
    for aa in AMINO_ACIDS:
        feats[f"aac_{aa}"] = seq.count(aa) / n

    # ── Dipeptide ─────────────────────────────────────────────────────────────
    dp_counts = {}
    for i in range(n - 1):
        dp = seq[i:i+2]
        if all(c in CFG.SEQ_VALID_AA for c in dp):
            dp_counts[dp] = dp_counts.get(dp, 0) + 1
    total_dp = max(1, n - 1)
    for a1 in AMINO_ACIDS:
        for a2 in AMINO_ACIDS:
            dp = a1 + a2
            feats[f"dp_{dp}"] = dp_counts.get(dp, 0) / total_dp

    # ── Pseudo-AAC (λ=1 sequence-order correlation) ───────────────────────────
    # Hydrophobicity scale (Kyte-Doolittle)
    kd = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
          "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
          "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
    if n > 1:
        corr = sum(kd.get(seq[i], 0) * kd.get(seq[i+1], 0)
                   for i in range(n-1)) / (n - 1)
    else:
        corr = 0.0
    for aa in AMINO_ACIDS:
        feats[f"pse_{aa}"] = feats[f"aac_{aa}"] / (1 + corr) if (1 + corr) != 0 else 0.0

    # ── Physicochemical ───────────────────────────────────────────────────────
    charged_pos = set("RKH")
    charged_neg = set("DE")
    aliphatic    = set("AVILM")
    aromatic     = set("FYW")
    polar        = set("STNQ")
    tiny         = set("AGCS")
    mw_approx = {"A":89,"R":174,"N":132,"D":133,"C":121,"Q":146,"E":147,
                 "G":75,"H":155,"I":131,"L":131,"K":128,"M":149,"F":165,
                 "P":115,"S":105,"T":119,"W":204,"Y":181,"V":117}

    feats["phys_charge_pos"]  = sum(seq.count(aa) for aa in charged_pos) / n
    feats["phys_charge_neg"]  = sum(seq.count(aa) for aa in charged_neg) / n
    feats["phys_net_charge"]  = feats["phys_charge_pos"] - feats["phys_charge_neg"]
    feats["phys_aliphatic"]   = sum(seq.count(aa) for aa in aliphatic) / n
    feats["phys_aromatic"]    = sum(seq.count(aa) for aa in aromatic) / n
    feats["phys_polar"]       = sum(seq.count(aa) for aa in polar) / n
    feats["phys_tiny"]        = sum(seq.count(aa) for aa in tiny) / n
    feats["phys_hydrophob"]   = sum(kd.get(aa, 0) * seq.count(aa)
                                    for aa in AMINO_ACIDS) / n
    feats["phys_mw_norm"]     = sum(mw_approx.get(aa, 110) * seq.count(aa)
                                    for aa in AMINO_ACIDS) / (n * 150)
    feats["phys_length"]      = np.log10(n)
    feats["phys_length_raw"]  = float(n)
    # Repeat composition (low complexity indicator)
    feats["phys_max_repeat"]  = max(seq.count(aa) for aa in AMINO_ACIDS) / n

    # Cys ratio (disulfide potential)
    feats["phys_cys_ratio"]   = seq.count("C") / n
    # Pro ratio (structural rigidity)
    feats["phys_pro_ratio"]   = seq.count("P") / n
    # Gly ratio (flexibility)
    feats["phys_gly_ratio"]   = seq.count("G") / n
    # His ratio (metal binding)
    feats["phys_his_ratio"]   = seq.count("H") / n
    # Asp+Glu (acidic)
    feats["phys_acidic"]      = (seq.count("D") + seq.count("E")) / n
    # Lys+Arg (basic, non-His)
    feats["phys_basic"]       = (seq.count("K") + seq.count("R")) / n
    # Instability index proxy
    feats["phys_instability"]  = (feats["phys_charge_pos"] +
                                  feats["phys_aromatic"]) / max(feats["phys_polar"] + 0.01, 0.01)

    # ── Catalytic core motifs (inv_ features) ─────────────────────────────────
    # Biotin-binding AMKM motif
    feats["inv_amkm"]  = 1.0 if "AMKM" in seq else 0.0
    # RuBisCO large subunit signature
    feats["inv_rubisco_sig"] = 1.0 if any(m in seq for m in ["GHYLNATAGTCE", "RIKFGETP"]) else 0.0
    # Phosphoenolpyruvate carboxylase signature
    feats["inv_pepc_sig"] = 1.0 if "QNTG" in seq else 0.0
    # Carbonic anhydrase zinc-binding
    feats["inv_ca_zinc"] = 1.0 if seq.count("H") >= 3 and n < 400 else 0.0
    # Generic ATP-grasp (biotin-dependent carboxylases)
    feats["inv_atp_grasp"] = 1.0 if any(m in seq for m in ["GRGREA", "GRGREP"]) else 0.0

    # ── EC motif placeholders (set 0 — unknown at query time) ─────────────────
    # These 17 features exist in the training matrix; we set 0 for new queries.
    for ec_short in ["4111", "4211", "6341", "6316", "6355", "6341b",
                     "4114", "6333", "4113", "41112", "4113b", "6411",
                     "6412", "6413", "6414", "4113c", "41149"]:
        feats[f"motif_{ec_short}"] = 0.0

    return feats


# ══════════════════════════════════════════════════════════════════════════════
# 3. Pfam features via HMMER
# ══════════════════════════════════════════════════════════════════════════════

def run_hmmer_pfam(seq_id: str, seq: str, tmp_dir: Path) -> dict:
    """
    Run hmmscan against Pfam-A.hmm.
    Returns dict: pfam_PFXXXXX → 1/0 for each domain in CARBOXY_PFAM_LIST,
    plus pfam_n_hits (count of any Pfam hits).
    """
    feats = {f"pfam_{pfam}": 0 for pfam in CARBOXY_PFAM_LIST}
    feats["pfam_n_hits"] = 0
    pfam_hits = []

    if not PATHS.PFAM_HMM.exists():
        log.warning("Pfam HMM not found at %s — skipping Pfam features", PATHS.PFAM_HMM)
        return feats, pfam_hits

    # Write temp FASTA
    fasta_tmp = tmp_dir / "query.fasta"
    with open(fasta_tmp, "w") as f:
        f.write(f">{seq_id}\n{seq}\n")

    out_tmp = tmp_dir / "hmmscan_pfam.tbl"

    try:
        cmd = [
            "hmmscan",
            "--domtblout", str(out_tmp),
            "-E", str(CFG.HMMER_EVALUE),
            "--cpu", "4",
            "--noali",
            str(PATHS.PFAM_HMM),
            str(fasta_tmp),
        ]
        result = subprocess.run(cmd, capture_output=True, timeout=120)
        if result.returncode != 0:
            log.warning("hmmscan failed: %s", result.stderr.decode()[:200])
            return feats, pfam_hits

        # Parse domtblout
        with open(out_tmp) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 4:
                    continue
                pfam_acc = parts[1].split(".")[0]  # e.g. PF00016
                evalue   = float(parts[12])
                if evalue > CFG.HMMER_EVALUE:
                    continue
                pfam_hits.append(pfam_acc)
                feats["pfam_n_hits"] += 1
                if pfam_acc in CFG.CARBOXY_PFAM:
                    feats[f"pfam_{pfam_acc}"] = 1

    except subprocess.TimeoutExpired:
        log.warning("hmmscan timed out")
    except FileNotFoundError:
        log.warning("hmmscan not found in PATH — skipping Pfam features")

    return feats, list(set(pfam_hits))


# ══════════════════════════════════════════════════════════════════════════════
# 4. InterPro features (simplified: count hits per database)
# ══════════════════════════════════════════════════════════════════════════════

def compute_interpro_features(pfam_hits: list) -> dict:
    """
    For Script 11 we derive InterPro features from the Pfam hits we already
    computed (same HMMER run), matching the 5 features used in training:
      n_pfam_hits, n_panther_hits, n_tigrfam_hits, n_cath_hits, n_superfamily_hits

    PANTHER/TIGRFAM/CATH/SUPERFAMILY require separate HMM databases not
    always present. We fill with 0 and flag a warning if needed.
    """
    feats = {col: 0 for col in INTERPRO_COLS}
    feats["n_pfam_hits"] = len(pfam_hits)
    # Others remain 0 unless full InterPro scan is available
    return feats


# ══════════════════════════════════════════════════════════════════════════════
# 5. ESM-2 embedding (1280 dims)
# ══════════════════════════════════════════════════════════════════════════════

def compute_esm2(seq_id: str, seq: str) -> np.ndarray:
    """
    Compute mean-pooled ESM-2 embedding (1280-dim).
    Truncates sequence to 1022 aa (ESM-2 limit).
    Returns np.ndarray shape (1280,).
    """
    try:
        import esm
        import torch
    except ImportError:
        log.error("ESM not installed. Run: pip install fair-esm")
        raise

    seq_trunc = seq[:1022]
    device = "cuda" if __import__("torch").cuda.is_available() else "cpu"

    log.info("  Loading ESM-2 model on %s ...", device)
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    data = [(seq_id, seq_trunc)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    with __import__("torch").no_grad():
        results = model(tokens, repr_layers=[33], return_contacts=False)

    # Mean pool over sequence positions (exclude BOS/EOS tokens)
    token_rep = results["representations"][33]
    seq_len   = len(seq_trunc)
    embedding = token_rep[0, 1:seq_len+1].mean(0).cpu().numpy()

    return embedding  # shape (1280,)


# ══════════════════════════════════════════════════════════════════════════════
# 6. Feature vector assembly
# ══════════════════════════════════════════════════════════════════════════════

def load_feature_names(task: str, suffix: str = "") -> list:
    p = ML_DIR / f"feature_names_{task}{suffix}.json"
    if p.exists():
        import json
        return json.load(open(p))
    log.warning("Feature names file not found: %s", p)
    return []


def assemble_feature_vector(comp_feats: dict,
                             pfam_feats: dict,
                             interpro_feats: dict,
                             esm2_emb: np.ndarray,
                             feat_names: list) -> np.ndarray:
    """
    Assemble features in exactly the order defined by feat_names
    (loaded from ML_DIR/feature_names_binary.json).
    Unknown features default to 0.
    """
    all_feats = {}
    all_feats.update(comp_feats)
    all_feats.update(pfam_feats)
    all_feats.update(interpro_feats)

    if esm2_emb is not None:
        for i, v in enumerate(esm2_emb):
            all_feats[f"esm2_{i}"] = float(v)
    else:
        # Fill zeros so vector shape is correct
        for i in range(CFG.ESM2_DIM):
            all_feats[f"esm2_{i}"] = 0.0

    if feat_names:
        vec = np.array([all_feats.get(f, 0.0) for f in feat_names],
                       dtype=np.float32)
    else:
        # Fallback: sort by key for reproducibility
        vec = np.array([all_feats[k] for k in sorted(all_feats)],
                       dtype=np.float32)

    return vec


def assemble_km_vector(binary_vec: np.ndarray,
                       ec_predicted: str,
                       kingdom: str,
                       feat_names_km: list) -> np.ndarray:
    """
    Build the Km feature vector = binary features + EC one-hot + kingdom.
    Matches the v3 Km matrix format used in training.
    """
    # EC one-hot (17 classes = KM_TRAINABLE_EC + 7 excluded → use training order)
    ec_oh_names = [f"ec_oh_{ec}" for ec in KM_TRAINABLE_EC]
    ec_oh = {name: 0.0 for name in ec_oh_names}
    ec_key = f"ec_oh_{ec_predicted}"
    if ec_key in ec_oh:
        ec_oh[ec_key] = 1.0

    # Kingdom one-hot (4 categories from training)
    kingdoms = ["bacteria", "plant", "archaea", "fungi"]
    kingdom_oh = {f"kingdom_{k}": 0.0 for k in kingdoms}
    kkey = f"kingdom_{kingdom.lower()}"
    if kkey in kingdom_oh:
        kingdom_oh[kkey] = 1.0

    # Combine with binary feature vector
    base_feats = {}
    if feat_names_km:
        binary_feat_names = load_feature_names("binary")
        for i, name in enumerate(binary_feat_names):
            if i < len(binary_vec):
                base_feats[name] = float(binary_vec[i])

    base_feats.update(ec_oh)
    base_feats.update(kingdom_oh)

    if feat_names_km:
        vec = np.array([base_feats.get(f, 0.0) for f in feat_names_km],
                       dtype=np.float32)
    else:
        vec = np.array(list(base_feats.values()), dtype=np.float32)

    return vec


# ══════════════════════════════════════════════════════════════════════════════
# 7. XGBoost prediction
# ══════════════════════════════════════════════════════════════════════════════

def load_models():
    import xgboost as xgb
    models = {}
    for name, fname in [("binary", "binary_v5.json"),
                        ("ec",     "ec_v5.json"),
                        ("km",     "km_v5_weighted.json")]:
        path = MODEL_DIR / fname
        if path.exists():
            b = xgb.Booster()
            b.load_model(str(path))
            models[name] = b
            log.info("  Loaded model: %s", fname)
        else:
            log.warning("  Model not found: %s", path)
            models[name] = None
    return models


def predict_binary(booster, vec: np.ndarray) -> tuple:
    import xgboost as xgb
    dmat = xgb.DMatrix(vec.reshape(1, -1))
    prob = float(booster.predict(dmat)[0])
    is_carb = prob >= 0.5
    # Confidence tiers matching script 10
    if prob >= 0.90:
        conf = "high"
    elif prob >= 0.70:
        conf = "medium"
    elif prob >= 0.50:
        conf = "low"
    else:
        conf = "non_carboxylase"
    return is_carb, prob, conf


def predict_ec(booster, vec: np.ndarray, ec_map_inv: dict) -> tuple:
    import xgboost as xgb
    dmat  = xgb.DMatrix(vec.reshape(1, -1))
    probs = booster.predict(dmat)[0]  # shape (n_classes,)
    top_idx   = int(np.argmax(probs))
    ec_pred   = ec_map_inv.get(top_idx, f"class_{top_idx}")
    ec_prob   = float(probs[top_idx])

    # Top-5 EC probabilities
    top5_idx = np.argsort(probs)[::-1][:5]
    ec_probs  = {ec_map_inv.get(int(i), f"class_{i}"): round(float(probs[i]), 4)
                 for i in top5_idx}

    return ec_pred, ec_prob, ec_probs


def predict_km(booster, vec: np.ndarray) -> tuple:
    import xgboost as xgb
    dmat    = xgb.DMatrix(vec.reshape(1, -1))
    log10_km = float(booster.predict(dmat)[0])
    km_mM    = float(10 ** log10_km)
    return km_mM, log10_km


# ══════════════════════════════════════════════════════════════════════════════
# 8. Main annotation pipeline
# ══════════════════════════════════════════════════════════════════════════════

def annotate_sequence(seq_id: str,
                      seq: str,
                      use_esm2: bool = True,
                      kingdom: str = "Bacteria") -> dict:
    """
    Full annotation pipeline for a single sequence.
    Returns result dict.
    """
    t0 = time.time()
    warnings = validate_sequence(seq_id, seq)
    seq_clean = "".join(c for c in seq.upper() if c in CFG.SEQ_VALID_AA)

    result = {
        "cdb_query_id":           seq_id,
        "sequence_length":        len(seq_clean),
        "is_carboxylase":         None,
        "carboxylase_probability": None,
        "confidence":             None,
        "ec_predicted":           None,
        "ec_name":                None,
        "ec_probabilities":       {},
        "km_predicted_mM":        None,
        "km_predicted_log10":     None,
        "km_ec_used":             None,
        "features_used":          [],
        "pfam_hits":              [],
        "warnings":               warnings,
        "runtime_seconds":        None,
    }

    # Load models & feature name lists
    log.info("Loading models...")
    models = load_models()
    if models["binary"] is None:
        result["warnings"].append("Binary model not found — cannot annotate")
        return result

    feat_names_binary = load_feature_names("binary")
    feat_names_ec     = load_feature_names("ec")
    feat_names_km     = load_feature_names("km", "_v3")

    # Load EC label map
    ec_map_path = ML_DIR / "ec_label_map_fixed.json"
    ec_map_inv  = {}
    if ec_map_path.exists():
        ec_map     = json.load(open(ec_map_path))
        ec_map_inv = {v: k for k, v in ec_map.items()}

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)

        # Step 1: Composition
        log.info("Computing composition features...")
        comp_feats = compute_composition(seq_clean)
        result["features_used"].append("composition")

        # Step 2: Pfam
        log.info("Running HMMER/Pfam scan...")
        pfam_feats, pfam_hits = run_hmmer_pfam(seq_id, seq_clean, tmp_path)
        result["pfam_hits"] = pfam_hits
        result["features_used"].append("pfam")

        # Step 3: InterPro
        interpro_feats = compute_interpro_features(pfam_hits)
        result["features_used"].append("interpro")

        # Step 4: ESM-2
        esm2_emb = None
        if use_esm2:
            log.info("Computing ESM-2 embedding (this takes ~10-30s)...")
            try:
                esm2_emb = compute_esm2(seq_id, seq_clean)
                result["features_used"].append("esm2")
            except Exception as e:
                log.warning("ESM-2 failed: %s — continuing without embedding", e)
                result["warnings"].append(f"ESM-2 failed: {e}")
        else:
            log.info("Skipping ESM-2 (--no-esm2 mode)")
            result["warnings"].append("ESM-2 skipped — predictions may be less accurate")

        # Step 5: Assemble binary feature vector
        log.info("Assembling feature vector...")
        binary_vec = assemble_feature_vector(
            comp_feats, pfam_feats, interpro_feats, esm2_emb, feat_names_binary)

        # Step 6: Binary prediction
        log.info("Predicting carboxylase activity...")
        is_carb, carb_prob, conf = predict_binary(models["binary"], binary_vec)

        result["is_carboxylase"]           = bool(is_carb)
        result["carboxylase_probability"]  = round(carb_prob, 4)
        result["confidence"]               = conf

        log.info("  is_carboxylase=%s  prob=%.4f  confidence=%s",
                 is_carb, carb_prob, conf)

        # Step 7: EC prediction (always run, even if not carboxylase)
        if models["ec"] is not None and ec_map_inv:
            # EC model uses same feature vector as binary
            ec_vec = assemble_feature_vector(
                comp_feats, pfam_feats, interpro_feats, esm2_emb, feat_names_ec)
            ec_pred, ec_prob, ec_probs = predict_ec(
                models["ec"], ec_vec, ec_map_inv)

            result["ec_predicted"]    = ec_pred
            result["ec_name"]         = EC_NAMES.get(ec_pred, ec_pred)
            result["ec_probabilities"] = ec_probs
            log.info("  ec_pred=%s  prob=%.4f", ec_pred, ec_prob)
        else:
            result["warnings"].append("EC model not available")

        # Step 8: Km prediction (only if carboxylase + EC is in trainable set)
        if is_carb and models["km"] is not None:
            ec_for_km = result["ec_predicted"]
            if ec_for_km in KM_TRAINABLE_EC:
                km_vec = assemble_km_vector(
                    binary_vec, ec_for_km, kingdom, feat_names_km)
                km_mM, km_log10 = predict_km(models["km"], km_vec)
                result["km_predicted_mM"]   = round(km_mM, 4)
                result["km_predicted_log10"] = round(km_log10, 4)
                result["km_ec_used"]         = ec_for_km
                log.info("  km_pred=%.4f mM (log10=%.3f)", km_mM, km_log10)
            else:
                result["warnings"].append(
                    f"Km prediction not available for EC {ec_for_km} "
                    f"(not in trainable set: {KM_TRAINABLE_EC})")

    result["runtime_seconds"] = round(time.time() - t0, 2)
    return result


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    import json

    ap = argparse.ArgumentParser(
        description="Annotate a protein sequence with CarboDB v5 models.")
    ap.add_argument("--fasta",    type=Path,
                    help="Input FASTA file (single or multi-sequence)")
    ap.add_argument("--sequence", type=str,
                    help="Raw amino acid sequence string")
    ap.add_argument("--out",      type=Path, default=None,
                    help="Output JSON file (default: stdout)")
    ap.add_argument("--kingdom",  type=str, default="Bacteria",
                    choices=["bacteria", "plant", "archaea", "fungi"],
                    help="Taxonomic kingdom for Km prediction (default: Bacteria)")
    ap.add_argument("--no-esm2",  action="store_true",
                    help="Skip ESM-2 embedding (faster but less accurate)")
    ap.add_argument("--pretty",   action="store_true",
                    help="Pretty-print JSON output")
    args = ap.parse_args()

    if not args.fasta and not args.sequence:
        ap.error("Provide --fasta FILE or --sequence SEQ")

    # Parse FASTA
    records = parse_fasta(fasta_path=args.fasta, sequence=args.sequence)
    log.info("Parsed %d sequence(s)", len(records))

    results = []
    for seq_id, seq in records:
        log.info("── Annotating: %s (len=%d) ──", seq_id, len(seq))
        res = annotate_sequence(
            seq_id, seq,
            use_esm2 = not args.no_esm2,
            kingdom  = args.kingdom,
        )
        results.append(res)

    # Output
    indent = 2 if args.pretty else None
    output = results[0] if len(results) == 1 else results
    out_str = json.dumps(output, indent=indent)

    if args.out:
        args.out.write_text(out_str)
        log.info("Saved: %s", args.out)
    else:
        print(out_str)


if __name__ == "__main__":
    main()
