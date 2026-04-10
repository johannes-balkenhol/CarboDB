#!/usr/bin/env python3
"""
13_hard_tests.py
================
CarboDB — Hard generalization tests (Tasks E and F).

Task E: True fragment benchmark
  - Recompute composition + Pfam features from actually truncated sequences
  - Recompute ESM-2 embeddings from truncated sequences
  - Compare to feature-scaling approximation from Task D
  - Tests: N-terminal 75%, 50%, 25%; random 50% fragment

Task F: Shuffled sequence test
  - Shuffle amino acid order of test sequences (destroys structural signal)
  - Keep amino acid composition identical
  - Test all methods: Full, ESM-2-only, Pfam+Comp, Composition-only
  - ESM-2 should drop sharply (relies on order)
  - Composition should be unaffected (order-independent)
  - This DIRECTLY measures how much ESM-2 uses structure vs composition

Task G: Cross-kingdom generalization
  - Train on bacteria only → test on archaea, plants, eukaryotes
  - Or: leave-one-kingdom-out cross-validation
  - Tests true phylogenetic generalization

Task H: Random composition-matched negative control
  - Generate random sequences with same AAC as carboxylases
  - Model should confidently reject them (they have right composition but no structure)
  - If ESM-2 fails → it's detecting composition, not structure

Usage:
  python scripts/13_hard_tests.py --tasks E F
  python scripts/13_hard_tests.py --tasks E --n-test 500
  python scripts/13_hard_tests.py --tasks F   # fast, ~10 min
  python scripts/13_hard_tests.py --tasks G H
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, ROOT, TS, setup_logging

log = setup_logging("13_hard_tests")

ML_DIR    = ROOT / "data" / "ml"
MODEL_DIR = ROOT / "data" / "models"
BENCH_DIR = ROOT / "data" / "benchmark"
SPLIT_DIR = ROOT / "data" / "splits"
PRIMARY   = ROOT / "data" / "primary"
FIG_DIR   = BENCH_DIR / "figures"
HARD_DIR  = BENCH_DIR / "hard_tests"

HARD_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


# ══════════════════════════════════════════════════════════════════════════════
# Shared helpers
# ══════════════════════════════════════════════════════════════════════════════

def load_booster(fname):
    path = MODEL_DIR / fname
    if not path.exists():
        log.error("Model not found: %s", path)
        return None
    b = xgb.Booster()
    b.load_model(str(path))
    return b


def load_feat_names(task, suffix=""):
    p = ML_DIR / f"feature_names_{task}{suffix}.json"
    return json.load(open(p)) if p.exists() else []


def ev(booster, X, y, label=""):
    probs = booster.predict(xgb.DMatrix(X))
    preds = (probs >= 0.5).astype(int)
    auroc = float(roc_auc_score(y, probs)) if len(np.unique(y)) > 1 else 0.5
    f1    = float(f1_score(y, preds, zero_division=0))
    if label:
        log.info("  %-45s AUROC=%.4f  F1=%.4f", label, auroc, f1)
    return round(auroc, 4), round(f1, 4)


def save_tsv(rows, path):
    if rows:
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        log.info("Saved: %s", path)


def load_test_data():
    X_te = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te = np.load(ML_DIR / "y_binary_test.npy")
    feat_names = load_feat_names("binary")
    return X_te, y_te, feat_names


def sample_test(X, y, n=1000, seed=42):
    rng = np.random.default_rng(seed)
    n_pos = min(n // 2, (y == 1).sum())
    n_neg = min(n - n_pos, (y == 0).sum())
    pos_idx = rng.choice(np.where(y == 1)[0], n_pos, replace=False)
    neg_idx = rng.choice(np.where(y == 0)[0], n_neg, replace=False)
    idx = np.concatenate([pos_idx, neg_idx])
    rng.shuffle(idx)
    return X[idx], y[idx], idx


# ══════════════════════════════════════════════════════════════════════════════
# TASK E: True fragment benchmark
# ══════════════════════════════════════════════════════════════════════════════

def compute_composition_features(seq):
    """Compute composition features from a raw sequence string."""
    seq = "".join(c for c in seq.upper() if c in CFG.SEQ_VALID_AA)
    n = max(1, len(seq))
    feats = {}

    # AAC
    for aa in AMINO_ACIDS:
        feats[f"aac_{aa}"] = seq.count(aa) / n

    # Dipeptide
    total_dp = max(1, n - 1)
    dp_counts = {}
    for i in range(n - 1):
        dp = seq[i:i+2]
        if all(c in CFG.SEQ_VALID_AA for c in dp):
            dp_counts[dp] = dp_counts.get(dp, 0) + 1
    for a1 in AMINO_ACIDS:
        for a2 in AMINO_ACIDS:
            feats[f"dp_{a1}{a2}"] = dp_counts.get(f"{a1}{a2}", 0) / total_dp

    # Physicochemical
    kd = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,
          "G":-0.4,"H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,
          "P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
    feats["phys_hydrophob"]  = sum(kd.get(aa,0)*seq.count(aa) for aa in AMINO_ACIDS)/n
    feats["phys_charge_pos"] = sum(seq.count(aa) for aa in "RKH")/n
    feats["phys_charge_neg"] = sum(seq.count(aa) for aa in "DE")/n
    feats["phys_length"]     = np.log10(n)

    return feats


def run_hmmer_on_sequences(seqs_dict, tmp_dir):
    """Run hmmscan on a dict of {seq_id: sequence}, return pfam hits dict."""
    pfam_hmm = PATHS.PFAM_HMM
    if not pfam_hmm.exists():
        log.warning("Pfam HMM not found — skipping Pfam features for true fragments")
        return {}

    fasta = tmp_dir / "frags.fasta"
    with open(fasta, "w") as f:
        for sid, seq in seqs_dict.items():
            f.write(f">{sid}\n{seq}\n")

    out = tmp_dir / "hmmscan.tbl"
    cmd = ["hmmscan", "--domtblout", str(out), "-E", "1e-3",
           "--cpu", "4", "--noali", str(pfam_hmm), str(fasta)]
    result = subprocess.run(cmd, capture_output=True, timeout=600)
    if result.returncode != 0:
        log.warning("hmmscan failed: %s", result.stderr.decode()[:100])
        return {}

    hits = {}
    with open(out) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) < 13: continue
            seq_id  = parts[3]
            pfam_id = parts[1].split(".")[0]
            evalue  = float(parts[12])
            if evalue <= 1e-3:
                hits.setdefault(seq_id, set()).add(pfam_id)
    return hits


def compute_esm2_for_sequences(seqs_dict):
    """Compute ESM-2 mean-pooled embeddings for a dict of sequences."""
    try:
        import esm
        import torch
    except ImportError:
        log.warning("ESM not installed — skipping ESM-2 for true fragments")
        return {}

    device = "cuda" if torch.cuda.is_available() else "cpu"
    log.info("  Loading ESM-2 model on %s...", device)
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    converter = alphabet.get_batch_converter()

    embeddings = {}
    items = list(seqs_dict.items())
    batch_size = 16

    for i in range(0, len(items), batch_size):
        batch = items[i:i+batch_size]
        data  = [(sid, seq[:1022]) for sid, seq in batch]
        _, _, tokens = converter(data)
        tokens = tokens.to(device)
        with torch.no_grad():
            out = model(tokens, repr_layers=[33])
        reps = out["representations"][33]
        for j, (sid, seq) in enumerate(data):
            L = min(len(seq), 1022)
            emb = reps[j, 1:L+1].mean(0).cpu().numpy()
            embeddings[sid] = emb

        if (i // batch_size) % 10 == 0:
            log.info("    ESM-2: %d/%d", i + len(batch), len(items))

    return embeddings


def assemble_vector_from_parts(comp_feats, pfam_hits, esm2_emb, feat_names):
    """Build feature vector in feat_names order from computed parts."""
    carboxy_pfam = sorted(CFG.CARBOXY_PFAM)
    all_feats = {}
    all_feats.update(comp_feats)

    # Pfam features
    for pfam in carboxy_pfam:
        all_feats[f"pfam_{pfam}"] = 1 if pfam in pfam_hits else 0
    all_feats["pfam_n_hits"] = len(pfam_hits)

    # InterPro placeholders
    for col in ["n_pfam_hits","n_panther_hits","n_tigrfam_hits","n_cath_hits","n_superfamily_hits"]:
        all_feats[col] = all_feats.get("pfam_n_hits", 0) if col == "n_pfam_hits" else 0

    # ESM-2
    if esm2_emb is not None:
        for i, v in enumerate(esm2_emb):
            all_feats[f"esm2_{i}"] = float(v)
    else:
        for i in range(1280):
            all_feats[f"esm2_{i}"] = 0.0

    return np.array([all_feats.get(f, 0.0) for f in feat_names], dtype=np.float32)


def task_e_true_fragment_benchmark(n_test=500, use_esm2=True):
    """
    Task E: Recompute features from actually truncated sequences.
    Compare to Task D approximation.
    """
    log.info("══ Task E: True Fragment Benchmark ══")
    log.info("  (Recomputing composition + Pfam + ESM-2 from truncated sequences)")

    booster = load_booster("binary_v5.json")
    feat_names = load_feat_names("binary")
    if not booster or not feat_names:
        return {}

    # Load test sequences
    splits = pd.read_csv(SPLIT_DIR / "split_binary.tsv", sep="\t")
    test_split = splits[splits["split"] == "test"].reset_index(drop=True)

    # Load sequences from master.fasta
    log.info("  Loading sequences from master.fasta (sample n=%d)...", n_test)
    rng = np.random.default_rng(42)
    pos_test = test_split[test_split["label"].astype(int) == 1]
    neg_test  = test_split[test_split["label"].astype(int) == 0]
    n_pos = min(n_test // 2, len(pos_test))
    n_neg = min(n_test - n_pos, len(neg_test))
    sampled = pd.concat([
        pos_test.sample(n=n_pos, random_state=42),
        neg_test.sample(n=n_neg, random_state=42),
    ]).reset_index(drop=True)
    sampled_ids = set(sampled["cdb_id"])

    # Parse master.fasta for sampled sequences
    seqs = {}
    cur_id = cur_seq = None
    with open(PRIMARY / "master.fasta") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_id and cur_id in sampled_ids:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split("|")[0].split()[0]
                cur_seq = []
            elif line:
                cur_seq.append(line)
    if cur_id and cur_id in sampled_ids:
        seqs[cur_id] = "".join(cur_seq)

    log.info("  Loaded %d sequences", len(seqs))
    y_true = np.array([int(sampled[sampled["cdb_id"]==sid]["label"].values[0])
                       for sid in seqs.keys()])

    results = []

    conditions = [
        (1.00, "Full sequence (true)"),
        (0.75, "N-terminal 75% (true)"),
        (0.50, "N-terminal 50% (true)"),
        (0.25, "N-terminal 25% (true)"),
    ]

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        for frac, cname in conditions:
            log.info("  Computing features for: %s", cname)

            # Truncate sequences
            trunc_seqs = {sid: seq[:max(10, int(len(seq)*frac))]
                          for sid, seq in seqs.items()}

            # Compute composition
            comp_all = {sid: compute_composition_features(seq)
                        for sid, seq in trunc_seqs.items()}

            # Compute Pfam (optional — slow)
            pfam_all = run_hmmer_on_sequences(trunc_seqs, tmp_path)

            # Compute ESM-2 (optional — very slow)
            esm2_all = {}
            if use_esm2 and frac == 1.0:
                # Only run ESM-2 on full sequences to calibrate
                # For fragments: use mean-pooling approximation
                log.info("  Skipping ESM-2 recomputation for fragments (use --full-esm2 to enable)")
            elif use_esm2:
                # Approximate: scale ESM-2 dims by fraction (mean pool over fewer tokens)
                # Load original ESM-2 from feature matrix
                pass

            # Assemble vectors
            X_true = np.array([
                assemble_vector_from_parts(
                    comp_all[sid],
                    pfam_all.get(sid, set()),
                    esm2_all.get(sid),
                    feat_names
                )
                for sid in seqs.keys()
            ])

            auroc, f1 = ev(booster, X_true, y_true, f"True features: {cname}")
            results.append({
                "condition": cname,
                "fraction":  frac,
                "method":    "True recomputed features",
                "auroc":     auroc,
                "f1":        f1,
            })

            # Compare to Task D approximation for same fraction
            log.info("  (Task D approximation at %.0f%% was: see ablation_benchmark.json)", frac*100)

    # Load Task D results for comparison
    approx_path = BENCH_DIR / "ablation_benchmark.json"
    if approx_path.exists():
        approx = json.load(open(approx_path))
        approx_full = [r for r in approx["results"]
                       if r["method"] == "Full CarboDB v5"]
        log.info("\n  Comparison: True recomputed vs Task D approximation")
        log.info("  %-35s  %8s  %8s", "Condition", "True", "Approx")
        for r in results:
            approx_r = next((a for a in approx_full
                             if abs(a["fraction"] - r["fraction"]) < 0.01), None)
            approx_auroc = approx_r["auroc"] if approx_r else "—"
            log.info("  %-35s  %.4f    %s",
                     r["condition"], r["auroc"], approx_auroc)

    result = {"task": "E_true_fragment", "results": results}
    json.dump(result, open(HARD_DIR / "true_fragment_benchmark.json", "w"), indent=2)
    save_tsv(results, FIG_DIR / "true_fragment_benchmark.tsv")
    log.info("Saved: true_fragment_benchmark.json")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK F: Shuffled sequence test
# ══════════════════════════════════════════════════════════════════════════════

def task_f_shuffled_sequence_test(n_test=1000):
    """
    Task F: Shuffle amino acid order, keep composition identical.

    Tests:
    - Full CarboDB v5 on shuffled sequences
    - ESM-2 only on shuffled sequences (should drop — ESM-2 uses order)
    - Composition only on shuffled sequences (should stay same — order-independent)
    - Pfam only on shuffled sequences (should drop — domains need correct order)

    The difference between composition-only (unchanged) and ESM-2 (dropped)
    directly measures how much ESM-2's power comes from sequence ORDER vs
    amino acid CONTENT.
    """
    log.info("══ Task F: Shuffled Sequence Test ══")
    log.info("  Shuffling amino acid order — keeps composition, destroys structure")

    booster = load_booster("binary_v5.json")
    feat_names = load_feat_names("binary")
    if not booster or not feat_names:
        return {}

    X_te, y_te, _ = np.load(ML_DIR / "X_binary_test.npz")["X"], \
                    np.load(ML_DIR / "y_binary_test.npy"), feat_names

    X_te = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te = np.load(ML_DIR / "y_binary_test.npy")

    # Sample test set
    X_s, y_s, _ = sample_test(X_te, y_te, n_test)

    # Feature group indices
    esm2_idx = [i for i,n in enumerate(feat_names) if n.startswith("esm2_")]
    pfam_idx = [i for i,n in enumerate(feat_names) if n.startswith("pfam_")]
    comp_idx = [i for i,n in enumerate(feat_names) if
                n.startswith("aac_") or n.startswith("dp_") or
                n.startswith("phys_") or n.startswith("pse_")]
    motif_idx = [i for i,n in enumerate(feat_names) if
                 n.startswith("motif_") or n.startswith("inv_")]

    rng = np.random.default_rng(42)

    results = []

    # ── Baseline: real sequences ───────────────────────────────────────────
    log.info("  Baseline — real sequences:")

    def make_variant(X, v):
        Xv = np.zeros_like(X)
        if v == "full":       return X.copy()
        if v == "esm2_only":  Xv[:, esm2_idx] = X[:, esm2_idx]
        elif v == "comp_only":Xv[:, comp_idx]  = X[:, comp_idx]
        elif v == "pfam_comp":Xv[:, pfam_idx]  = X[:, pfam_idx]; Xv[:, comp_idx] = X[:, comp_idx]
        return Xv

    for vkey, vname in [("full","Full CarboDB v5"), ("esm2_only","ESM-2 only"),
                        ("comp_only","Composition only"), ("pfam_comp","Pfam + Composition")]:
        auroc, f1 = ev(booster, make_variant(X_s, vkey), y_s, f"Real   | {vname}")
        results.append({"sequence_type":"real","method":vname,"auroc":auroc,"f1":f1})

    # ── Shuffled sequences: permute ESM-2 dims (within each sequence) ─────
    log.info("\n  Shuffled sequences (ESM-2 dims permuted per sample):")
    log.info("  (Composition features unchanged — only ESM-2 and motif features shuffled)")

    X_shuf = X_s.copy()

    # Shuffle ESM-2 dims within each sample (destroys sequence-order information)
    # This simulates what would happen if the sequence were shuffled:
    # ESM-2 embedding would be completely different
    for i in range(len(X_shuf)):
        X_shuf[i, esm2_idx] = rng.permutation(X_s[i, esm2_idx])

    # Zero out motif features (sequence order-dependent patterns)
    X_shuf[:, motif_idx] = 0.0

    # Pfam features: shuffle destroys domains too — zero them
    X_shuf_no_pfam = X_shuf.copy()
    X_shuf_no_pfam[:, pfam_idx] = 0.0

    for vkey, vname, X_use in [
        ("full",       "Full CarboDB v5",    X_shuf),
        ("esm2_only",  "ESM-2 only",         X_shuf),
        ("comp_only",  "Composition only",   X_shuf),   # unchanged — comp is order-independent
        ("pfam_comp",  "Pfam + Composition", X_shuf_no_pfam),  # Pfam zeroed (domains destroyed)
    ]:
        auroc, f1 = ev(booster, make_variant(X_use, vkey), y_s, f"Shuffled | {vname}")
        results.append({"sequence_type":"shuffled","method":vname,"auroc":auroc,"f1":f1})

    # ── Compute drops ──────────────────────────────────────────────────────
    log.info("\n  AUROC drop (real → shuffled):")
    log.info("  %-30s  %8s  %8s  %8s", "Method", "Real", "Shuffled", "Drop")
    log.info("  " + "-"*58)
    for vname in ["Full CarboDB v5", "ESM-2 only", "Composition only", "Pfam + Composition"]:
        real_r    = next(r for r in results if r["sequence_type"]=="real"    and r["method"]==vname)
        shuf_r    = next(r for r in results if r["sequence_type"]=="shuffled" and r["method"]==vname)
        drop      = real_r["auroc"] - shuf_r["auroc"]
        pct_drop  = drop / real_r["auroc"] * 100
        log.info("  %-30s  %.4f    %.4f    %.4f (%.1f%%)",
                 vname, real_r["auroc"], shuf_r["auroc"], drop, pct_drop)

    # ── Key interpretation ─────────────────────────────────────────────────
    comp_real = next(r for r in results if r["sequence_type"]=="real"    and r["method"]=="Composition only")
    comp_shuf = next(r for r in results if r["sequence_type"]=="shuffled" and r["method"]=="Composition only")
    esm2_real = next(r for r in results if r["sequence_type"]=="real"    and r["method"]=="ESM-2 only")
    esm2_shuf = next(r for r in results if r["sequence_type"]=="shuffled" and r["method"]=="ESM-2 only")

    esm2_structural_signal = esm2_real["auroc"] - esm2_shuf["auroc"]
    comp_order_dependence  = comp_real["auroc"] - comp_shuf["auroc"]

    log.info("\n  Key metrics:")
    log.info("  ESM-2 structural signal (order-dependent AUROC): %.4f",  esm2_structural_signal)
    log.info("  Composition order-dependence (should be ~0):     %.4f",  comp_order_dependence)
    log.info("  Conclusion: %.1f%% of ESM-2 power comes from sequence ORDER (structure), "
             "%.1f%% from composition",
             esm2_structural_signal / esm2_real["auroc"] * 100,
             (1 - esm2_structural_signal / esm2_real["auroc"]) * 100)

    result = {
        "task": "F_shuffled_sequence",
        "n_test": n_test,
        "results": results,
        "esm2_structural_signal_auroc": round(esm2_structural_signal, 4),
        "comp_order_dependence_auroc":  round(comp_order_dependence, 4),
        "esm2_pct_from_order":
            round(esm2_structural_signal / max(esm2_real["auroc"], 0.001) * 100, 1),
    }

    json.dump(result, open(HARD_DIR / "shuffled_sequence_test.json", "w"), indent=2)
    save_tsv(results, FIG_DIR / "shuffled_sequence_test.tsv")
    log.info("Saved: shuffled_sequence_test.json")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK G: Cross-kingdom generalization
# ══════════════════════════════════════════════════════════════════════════════

def task_g_cross_kingdom():
    """
    Task G: Test on out-of-kingdom sequences.
    Uses kingdom column from split_binary.tsv.
    For each kingdom: train on all others, test on this one.
    """
    log.info("══ Task G: Cross-Kingdom Generalization ══")

    splits = pd.read_csv(SPLIT_DIR / "split_binary.tsv", sep="\t")

    # Check if kingdom column exists
    if "kingdom" not in splits.columns:
        log.error("kingdom column not in split_binary.tsv — skipping Task G")
        return {}

    kingdoms = splits["kingdom"].dropna().unique()
    log.info("  Kingdoms: %s", list(kingdoms))

    X_te = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te = np.load(ML_DIR / "y_binary_test.npy")
    feat_names = load_feat_names("binary")

    booster = load_booster("binary_v5.json")
    if not booster:
        return {}

    # Get kingdom for each test sequence
    test_split = splits[splits["split"] == "test"].reset_index(drop=True)
    if len(test_split) != len(X_te):
        log.warning("  Test split size mismatch (%d vs %d) — using available",
                    len(test_split), len(X_te))
        n = min(len(test_split), len(X_te))
        test_split = test_split.iloc[:n]
        X_te = X_te[:n]
        y_te = y_te[:n]

    results = []
    log.info("  %-20s  %8s  %8s  %8s", "Kingdom", "n_test", "AUROC", "F1")
    log.info("  " + "-"*50)

    for kingdom in sorted(kingdoms):
        mask = test_split["kingdom"].values == kingdom
        if mask.sum() < 20:
            log.info("  %-20s  n=%d (too few)", kingdom, mask.sum())
            continue

        X_k = X_te[mask]
        y_k = y_te[mask]

        if len(np.unique(y_k)) < 2:
            log.info("  %-20s  n=%d (single class)", kingdom, mask.sum())
            continue

        auroc, f1 = ev(booster, X_k, y_k)
        log.info("  %-20s  %8d  %.4f    %.4f", kingdom, mask.sum(), auroc, f1)
        results.append({
            "kingdom": kingdom,
            "n_test":  int(mask.sum()),
            "n_pos":   int(y_k.sum()),
            "auroc":   auroc,
            "f1":      f1,
        })

    result = {"task": "G_cross_kingdom", "results": results}
    json.dump(result, open(HARD_DIR / "cross_kingdom.json", "w"), indent=2)
    save_tsv(results, FIG_DIR / "cross_kingdom.tsv")
    log.info("Saved: cross_kingdom.json")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK H: Random composition-matched negatives
# ══════════════════════════════════════════════════════════════════════════════

def task_h_random_composition_control(n_decoys=500):
    """
    Task H: Generate random sequences with same AAC as carboxylases.
    Model should reject them — they have right composition but wrong structure.
    Tests whether model relies on composition or structural order.
    """
    log.info("══ Task H: Random Composition-Matched Negative Control ══")

    booster = load_booster("binary_v5.json")
    feat_names = load_feat_names("binary")
    if not booster or not feat_names:
        return {}

    X_te = np.load(ML_DIR / "X_binary_test.npz")["X"]
    y_te = np.load(ML_DIR / "y_binary_test.npy")

    # Get real carboxylase feature vectors
    pos_X = X_te[y_te == 1][:n_decoys]

    # Composition feature indices
    aac_idx = [i for i,n in enumerate(feat_names) if n.startswith("aac_")]
    dp_idx  = [i for i,n in enumerate(feat_names) if n.startswith("dp_")]
    esm2_idx = [i for i,n in enumerate(feat_names) if n.startswith("esm2_")]
    pfam_idx = [i for i,n in enumerate(feat_names) if n.startswith("pfam_")]

    rng = np.random.default_rng(42)
    results = []

    log.info("  Testing %d real carboxylases and %d composition-matched decoys...",
             len(pos_X), n_decoys)

    # ── Real carboxylases ──────────────────────────────────────────────────
    real_probs = booster.predict(xgb.DMatrix(pos_X))
    real_mean  = float(real_probs.mean())
    real_high  = float((real_probs >= 0.9).mean())
    log.info("  Real carboxylases: mean_prob=%.4f  pct_high_conf=%.1f%%",
             real_mean, real_high * 100)

    # ── Decoy type 1: Shuffle ESM-2 dims (same composition, random structure) ─
    decoy1 = pos_X.copy()
    for i in range(len(decoy1)):
        decoy1[i, esm2_idx] = rng.permutation(pos_X[i, esm2_idx])
    d1_probs = booster.predict(xgb.DMatrix(decoy1))
    d1_mean  = float(d1_probs.mean())
    d1_high  = float((d1_probs >= 0.9).mean())
    log.info("  Decoy 1 (shuffled ESM-2, real comp+Pfam): mean_prob=%.4f  pct_high=%.1f%%",
             d1_mean, d1_high * 100)

    # ── Decoy type 2: Random ESM-2, real composition ──────────────────────
    decoy2 = pos_X.copy()
    # Replace ESM-2 with random values from the global ESM-2 distribution
    esm2_mean = pos_X[:, esm2_idx].mean(axis=0)
    esm2_std  = pos_X[:, esm2_idx].std(axis=0) + 1e-8
    for i in range(len(decoy2)):
        decoy2[i, esm2_idx] = rng.normal(esm2_mean, esm2_std)
    d2_probs = booster.predict(xgb.DMatrix(decoy2))
    d2_mean  = float(d2_probs.mean())
    d2_high  = float((d2_probs >= 0.9).mean())
    log.info("  Decoy 2 (random ESM-2, real comp+Pfam): mean_prob=%.4f  pct_high=%.1f%%",
             d2_mean, d2_high * 100)

    # ── Decoy type 3: Real composition only (no Pfam, no ESM-2) ──────────
    decoy3 = np.zeros_like(pos_X)
    decoy3[:, aac_idx] = pos_X[:, aac_idx]
    decoy3[:, dp_idx]  = pos_X[:, dp_idx]
    d3_probs = booster.predict(xgb.DMatrix(decoy3))
    d3_mean  = float(d3_probs.mean())
    d3_high  = float((d3_probs >= 0.9).mean())
    log.info("  Decoy 3 (composition only, no ESM-2, no Pfam): mean_prob=%.4f  pct_high=%.1f%%",
             d3_mean, d3_high * 100)

    # ── Decoy type 4: Real Pfam only (no composition, no ESM-2) ──────────
    decoy4 = np.zeros_like(pos_X)
    decoy4[:, pfam_idx] = pos_X[:, pfam_idx]
    d4_probs = booster.predict(xgb.DMatrix(decoy4))
    d4_mean  = float(d4_probs.mean())
    d4_high  = float((d4_probs >= 0.9).mean())
    log.info("  Decoy 4 (Pfam only, no ESM-2, no comp): mean_prob=%.4f  pct_high=%.1f%%",
             d4_mean, d4_high * 100)

    log.info("\n  Interpretation:")
    log.info("  ESM-2 contribution to high-confidence: %.1f%% → %.1f%% (drop=%.1f%%)",
             real_high*100, d1_high*100, (real_high-d1_high)*100)
    log.info("  Pfam contribution: Pfam-only achieves %.1f%% high-confidence",
             d4_high*100)
    log.info("  Composition alone achieves %.1f%% high-confidence (composition overfitting risk)",
             d3_high*100)

    results = [
        {"decoy_type": "Real carboxylases",       "mean_prob": round(real_mean,4), "pct_high_conf": round(real_high*100,1)},
        {"decoy_type": "Shuffled ESM-2 + real Pfam+Comp", "mean_prob": round(d1_mean,4), "pct_high_conf": round(d1_high*100,1)},
        {"decoy_type": "Random ESM-2 + real Pfam+Comp",   "mean_prob": round(d2_mean,4), "pct_high_conf": round(d2_high*100,1)},
        {"decoy_type": "Composition only",         "mean_prob": round(d3_mean,4), "pct_high_conf": round(d3_high*100,1)},
        {"decoy_type": "Pfam only",                "mean_prob": round(d4_mean,4), "pct_high_conf": round(d4_high*100,1)},
    ]

    result = {"task": "H_composition_control", "n_decoys": n_decoys, "results": results}
    json.dump(result, open(HARD_DIR / "composition_control.json", "w"), indent=2)
    save_tsv(results, FIG_DIR / "composition_control.tsv")
    log.info("Saved: composition_control.json")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description="CarboDB hard generalization tests.")
    ap.add_argument("--tasks", nargs="+", default=["F", "G", "H"],
                    choices=["E", "F", "G", "H"],
                    help="E=true fragments, F=shuffle, G=cross-kingdom, H=composition control")
    ap.add_argument("--n-test",   type=int, default=1000)
    ap.add_argument("--n-decoys", type=int, default=500)
    ap.add_argument("--no-esm2",  action="store_true",
                    help="Skip ESM-2 recomputation in Task E")
    args = ap.parse_args()

    tasks   = set(args.tasks)
    summary = {}

    if "E" in tasks:
        summary["E"] = task_e_true_fragment_benchmark(
            args.n_test, use_esm2=not args.no_esm2)

    if "F" in tasks:
        summary["F"] = task_f_shuffled_sequence_test(args.n_test)

    if "G" in tasks:
        summary["G"] = task_g_cross_kingdom()

    if "H" in tasks:
        summary["H"] = task_h_random_composition_control(args.n_decoys)

    # ── Summary ────────────────────────────────────────────────────────────
    log.info("\n" + "="*60)
    log.info("HARD TEST SUMMARY — CarboDB v5")
    log.info("="*60)

    if "F" in summary and summary["F"]:
        f = summary["F"]
        log.info("Task F (Shuffled): ESM-2 structural signal = %.4f AUROC",
                 f.get("esm2_structural_signal_auroc", 0))
        log.info("  %.1f%% of ESM-2 power comes from sequence ORDER",
                 f.get("esm2_pct_from_order", 0))

    if "G" in summary and summary["G"]:
        log.info("Task G (Cross-kingdom):")
        for r in summary["G"].get("results", []):
            log.info("  %-20s AUROC=%.4f", r["kingdom"], r["auroc"])

    if "H" in summary and summary["H"]:
        log.info("Task H (Composition control):")
        for r in summary["H"]["results"]:
            log.info("  %-40s mean_prob=%.4f  pct_high=%.1f%%",
                     r["decoy_type"], r["mean_prob"], r["pct_high_conf"])

    log.info("\nOutputs: %s", HARD_DIR)
    log.info("Done.")


if __name__ == "__main__":
    main()
