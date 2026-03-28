#!/usr/bin/env python3
"""
04a_composition.py
==================
CarboxyDB — Annotation Layer A: Sequence-derived features.

Computes all sequence features that need no external tools:
  A1  AA composition          (20 features, aac_ prefix)
  A2  Dipeptide frequency     (400 features, dp_ prefix)
  A3  Pseudo-AAC              (30 features, pse_ prefix)
  A4  Physicochemical         (~15 features, phys_ prefix)
  A5  Catalytic core          (~17 features, inv_cat_ prefix)
  A6  EC-specific motifs      (7 features, motif_ prefix)

Input:  data/primary/master.fasta  (or a chunk)
Output: data/features/composition/composition_{chunk}.tsv
        First column = cdb_id

Usage:
    # Full dataset (runs in ~30 min on single CPU)
    python scripts/04a_composition.py

    # Single chunk (for SLURM array — not needed, this script is fast enough)
    python scripts/04a_composition.py --chunk 0001

    # Test on first 100 sequences
    python scripts/04a_composition.py --test
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("04a_composition")

AA = sorted(CFG.SEQ_VALID_AA)   # 20 standard amino acids

# ── EC-specific motif patterns ────────────────────────────────────────────────
MOTIFS = {
    "motif_rubisco_kk":      re.compile(r"K[A-Z]{1,3}K"),
    "motif_rubisco_gk":      re.compile(r"G[A-Z]{0,2}K"),
    "motif_ca_hh":           re.compile(r"H[A-Z]{1,4}H"),
    "motif_ca_his_cluster":  re.compile(r"H[A-Z]{0,5}H[A-Z]{0,5}H"),
    "motif_pepc_rr":         re.compile(r"R[A-Z]{1,3}R"),
    "motif_biotin_mk":       re.compile(r"[LIVMF]X{0,3}[DE]X{3,5}[LIVMF]X{4,5}[KR]AMK",
                                         re.IGNORECASE),
    "motif_biotin_amk":      re.compile(r"AMK"),
}


# ═══════════════════════════════════════════════════════════════════════════════
# Feature extractors
# ═══════════════════════════════════════════════════════════════════════════════

def aa_composition(seq: str) -> dict:
    """A1: 20 amino acid frequencies."""
    n = len(seq)
    if n == 0:
        return {f"aac_{a}": 0.0 for a in AA}
    counts = {a: 0 for a in AA}
    for aa in seq:
        if aa in counts:
            counts[aa] += 1
    return {f"aac_{a}": counts[a] / n for a in AA}


def dipeptide_freq(seq: str) -> dict:
    """A2: 400 dipeptide frequencies."""
    n = len(seq) - 1
    if n <= 0:
        return {f"dp_{a}{b}": 0.0 for a in AA for b in AA}
    counts = {}
    for i in range(n):
        dp = seq[i:i+2]
        if all(c in CFG.SEQ_VALID_AA for c in dp):
            counts[dp] = counts.get(dp, 0) + 1
    return {f"dp_{a}{b}": counts.get(f"{a}{b}", 0) / n for a in AA for b in AA}


def pseudo_aac(seq: str, lamda: int = 10, weight: float = 0.05) -> dict:
    """
    A3: Pseudo amino acid composition (PseAAC).
    30 features: 20 standard AAC + 10 sequence-order correlation factors.
    """
    features = {}
    n = len(seq)

    # Physicochemical properties for correlation factors
    H1 = {"A":0.62,"C":0.29,"D":-0.90,"E":-0.74,"F":1.19,"G":0.48,"H":-0.40,
          "I":1.38,"K":-1.50,"L":1.06,"M":0.64,"N":-0.78,"P":0.12,"Q":-0.85,
          "R":-2.53,"S":-0.18,"T":-0.05,"V":1.08,"W":0.81,"Y":0.26}
    H2 = {"A":-0.5,"C":1.0,"D":3.0,"E":3.0,"F":-2.5,"G":0.0,"H":-0.5,
          "I":-1.8,"K":3.0,"L":-1.8,"M":-1.3,"N":0.2,"P":0.0,"Q":0.2,
          "R":3.0,"S":0.3,"T":-0.4,"V":-1.5,"W":-3.4,"Y":-2.3}

    # Standard AAC part (weighted)
    aac = aa_composition(seq)
    sum_theta = 0.0
    thetas = []
    for lam in range(1, min(lamda+1, n)):
        theta = 0.0
        for i in range(n - lam):
            a1, a2 = seq[i], seq[i+lam]
            if a1 in H1 and a2 in H1:
                theta += (H1[a1]-H1[a2])**2 + (H2[a1]-H2[a2])**2
        theta /= max(1, n - lam)
        thetas.append(theta)
        sum_theta += theta

    denom = 1 + weight * sum_theta
    for a in AA:
        features[f"pse_{a}"] = aac[f"aac_{a}"] / denom

    for j, theta in enumerate(thetas):
        features[f"pse_corr_{j+1}"] = (weight * theta) / denom

    # Pad if sequence too short for all lambda values
    for j in range(len(thetas), lamda):
        features[f"pse_corr_{j+1}"] = 0.0

    return features


def physicochemical(seq: str) -> dict:
    """A4: Physicochemical properties via BioPython ProteinAnalysis."""
    features = {}
    n = len(seq)
    features["phys_length"]      = float(n)
    features["phys_length_log"]  = float(np.log1p(n))

    try:
        pa = ProteinAnalysis(seq)
        features["phys_mw"]          = pa.molecular_weight()
        features["phys_pi"]          = pa.isoelectric_point()
        features["phys_gravy"]       = pa.gravy()
        features["phys_aromaticity"] = pa.aromaticity()
        features["phys_instability"] = pa.instability_index()
        features["phys_charge_ph7"]  = pa.charge_at_pH(7.0)
    except Exception:
        for k in ["phys_mw","phys_pi","phys_gravy","phys_aromaticity",
                  "phys_instability","phys_charge_ph7"]:
            features[k] = 0.0

    # Residue class fractions
    seq_set = list(seq)
    charged  = set("DEKR")
    aromatic = set("FWY")
    polar    = set("NQST")
    nonpolar = set("ACGILMPV")
    small    = set("AGST")

    features["phys_frac_charged"]  = sum(1 for a in seq_set if a in charged)  / max(n,1)
    features["phys_frac_aromatic"] = sum(1 for a in seq_set if a in aromatic) / max(n,1)
    features["phys_frac_polar"]    = sum(1 for a in seq_set if a in polar)    / max(n,1)
    features["phys_frac_nonpolar"] = sum(1 for a in seq_set if a in nonpolar) / max(n,1)
    features["phys_frac_small"]    = sum(1 for a in seq_set if a in small)    / max(n,1)
    features["phys_frac_glycine"]  = seq.count("G") / max(n,1)
    features["phys_frac_proline"]  = seq.count("P") / max(n,1)

    return features


def catalytic_core(seq: str) -> dict:
    """
    A5: Catalytic core features.
    Analyses the middle 50% of the sequence (core region).
    Confirmed top features in CarboxyPred v3: catalytic_K_percent (4.3%)
    """
    features = {}
    n = len(seq)
    start = n // 4
    end   = n - n // 4
    core  = seq[start:end]
    nc    = max(len(core), 1)

    cat_aa = list("DEHKCST")  # common catalytic residues
    for a in cat_aa:
        features[f"inv_cat_{a}"] = core.count(a) / nc

    # Distance statistics between catalytic residues
    positions = [i for i, aa in enumerate(core) if aa in set("DEHKCS")]
    if len(positions) >= 2:
        dists = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
        features["inv_cat_mean_dist"] = float(np.mean(dists))
        features["inv_cat_std_dist"]  = float(np.std(dists))
        features["inv_cat_min_dist"]  = float(np.min(dists))
        features["inv_cat_max_dist"]  = float(np.max(dists))
        features["inv_cat_clustering"] = len(positions) / nc
    else:
        for k in ["inv_cat_mean_dist","inv_cat_std_dist",
                  "inv_cat_min_dist","inv_cat_max_dist","inv_cat_clustering"]:
            features[k] = 0.0

    # Core composition classes
    hydrophobic = set("ACFILMPVWY")
    charged_aa  = set("DEKR")
    polar_aa    = set("NQSTH")
    aromatic_aa = set("FWY")

    features["inv_hydrophobic"] = sum(1 for a in core if a in hydrophobic) / nc
    features["inv_charged"]     = sum(1 for a in core if a in charged_aa)  / nc
    features["inv_polar"]       = sum(1 for a in core if a in polar_aa)    / nc
    features["inv_aromatic"]    = sum(1 for a in core if a in aromatic_aa) / nc

    pos_charged = sum(1 for a in core if a in "KR")
    neg_charged = sum(1 for a in core if a in "DE")
    features["inv_net_charge"]  = (pos_charged - neg_charged) / nc

    return features


def ec_motifs(seq: str) -> dict:
    """A6: EC-specific motif counts (7 features)."""
    return {name: float(len(pat.findall(seq))) for name, pat in MOTIFS.items()}


def extract_all(cdb_id: str, seq: str) -> dict:
    """Extract all A-layer features for one sequence."""
    row = {"cdb_id": cdb_id}
    row.update(aa_composition(seq))
    row.update(dipeptide_freq(seq))
    row.update(pseudo_aac(seq))
    row.update(physicochemical(seq))
    row.update(catalytic_core(seq))
    row.update(ec_motifs(seq))
    return row


# ═══════════════════════════════════════════════════════════════════════════════
# FASTA reader
# ═══════════════════════════════════════════════════════════════════════════════

def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read FASTA, return list of (cdb_id, sequence)."""
    seqs = []
    cdb_id = seq_lines = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cdb_id and seq_lines:
                    seqs.append((cdb_id, "".join(seq_lines)))
                cdb_id = line[1:].split("|")[0]  # first field = cdb_id
                seq_lines = []
            elif line:
                if seq_lines is not None:
                    seq_lines.append(line)
    if cdb_id and seq_lines:
        seqs.append((cdb_id, "".join(seq_lines)))
    return seqs


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk", default=None,
                    help="Chunk ID e.g. 0001 — process single chunk file")
    ap.add_argument("--test",  action="store_true",
                    help="Test on first 100 sequences only")
    ap.add_argument("--batch-size", type=int, default=10_000,
                    help="Sequences per output file when processing full FASTA")
    args = ap.parse_args()

    PATHS.FEAT_COMP.mkdir(parents=True, exist_ok=True)

    if args.chunk:
        # Single chunk mode (for SLURM array if needed)
        fasta_path = PATHS.INTERIM / "fasta_chunks" / f"chunk_{args.chunk}.fasta"
        if not fasta_path.exists():
            log.error("Chunk not found: %s", fasta_path)
            sys.exit(1)
        seqs = read_fasta(fasta_path)
        out_path = PATHS.FEAT_COMP / f"composition_{args.chunk}.tsv"
        process_batch(seqs, out_path)

    elif args.test:
        fasta_path = PATHS.MASTER_FASTA
        seqs = read_fasta(fasta_path)[:100]
        out_path = PATHS.FEAT_COMP / f"composition_test.tsv"
        process_batch(seqs, out_path)
        log.info("Test done — %s", out_path)

    else:
        # Full dataset — process in batches, write multiple TSV files
        log.info("Processing full master.fasta (%s)", PATHS.MASTER_FASTA)
        fasta_path = PATHS.MASTER_FASTA

        log.info("Reading FASTA...")
        seqs = read_fasta(fasta_path)
        log.info("  %d sequences loaded", len(seqs))

        batch_size = args.batch_size
        n_batches = (len(seqs) + batch_size - 1) // batch_size

        for i in range(n_batches):
            batch = seqs[i*batch_size : (i+1)*batch_size]
            out_path = PATHS.FEAT_COMP / f"composition_{i+1:04d}.tsv"

            if out_path.exists():
                log.info("  Batch %d/%d already exists — skip", i+1, n_batches)
                continue

            log.info("  Batch %d/%d (%d seqs)...", i+1, n_batches, len(batch))
            process_batch(batch, out_path)

        log.info("All batches done.")
        log.info("Next: concatenate with:")
        log.info("  python scripts/04a_composition.py --merge")


def process_batch(seqs: list, out_path: Path):
    rows = []
    for cdb_id, seq in tqdm(seqs, desc=str(out_path.name), leave=False):
        try:
            rows.append(extract_all(cdb_id, seq))
        except Exception as e:
            log.warning("  Error on %s: %s", cdb_id, e)

    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(out_path, sep="\t", index=False)
        log.info("  Written: %s (%d rows, %d features)",
                 out_path.name, len(df), len(df.columns)-1)


if __name__ == "__main__":
    main()
