#!/usr/bin/env python3
"""
04_annotate_features.py
=======================
CarboxyDB — Step 4: Extract ALL features for every sequence.

Feature layers (confirmed from working CarboxylaseDatabase runs, Dec 2025):
  A  — Sequence-derived (no external tools, always fast)
       A1: AA composition         (20 features, aa_ prefix)
       A2: Dipeptide frequency    (400 features, dp_ prefix) ← top ML features
       A3: Pseudo-AAC             (30 features, pse_ prefix)
       A4: Physicochemical        (~15 features, phys_ prefix)
       A5: Catalytic core         (~17 features, inv_cat_ prefix)
       A6: EC-specific motifs     (7 features, motif_ prefix)

  B  — Domain annotation (needs HMMER3 + Pfam-A.hmm)
       B1: Pfam binary hits       (~30 features, pfam_ prefix)
       B2: PROSITE regex          (14 features, prosite_ prefix)

  C  — Homology (needs BLAST+ and BRENDA BLAST db)
       C1: Best BLAST hit         (4 features, blast_ prefix)

  D  — MEME motif hits [PENDING]
       Binary FIMO hits           (65 features, meme_ prefix)
       Added automatically when data/features/meme/meme_hits.tsv exists

  E  — ESM-2 embeddings [HPC GPU job]
       1280 features, esm2_ prefix
       Run separately: python scripts/04b_esm2.py (submit to HPC)

IMPORTANT: Layers B and C run in BATCH mode, not per-sequence.
  B: hmmscan runs ONCE on full FASTA → parse domtblout → map to per-sequence dict
  C: blastp runs ONCE in batch mode → parse tabular output → map to per-sequence dict
  This is the critical difference from the broken approach of calling tools per-sequence.

Usage:
    # Run all layers A+B+C+D (D only if meme_hits.tsv exists)
    python scripts/04_annotate_features.py

    # Run only layer A (fast, for testing)
    python scripts/04_annotate_features.py --layers A

    # Run B only (after HMMER DB is ready)
    python scripts/04_annotate_features.py --layers B

    # Merge ESM-2 into existing feature file (after HPC job)
    python scripts/04_annotate_features.py --merge-esm2
"""

import argparse
import json
import math
import re
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("WARNING: biopython not installed — physicochemical features will be partial")

TS = datetime.now().strftime("%Y%m%d_%H%M%S")

INTERIM  = Path("data/interim")
FEAT_DIR = Path("data/processed")
FEAT_DIR.mkdir(parents=True, exist_ok=True)

PFAM_HMM  = Path("data/raw/pfam/Pfam-A.hmm")
BLAST_DB  = Path("data/features/blast/brenda_db")
MEME_HITS = Path("data/features/meme/meme_hits.tsv")

AA = "ACDEFGHIKLMNPQRSTVWY"


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER A — Sequence composition (pure Python, no external tools)
# ═══════════════════════════════════════════════════════════════════════════════

def feat_aa_composition(seq: str) -> dict:
    """20 amino acid frequencies."""
    n = max(len(seq), 1)
    return {f"aa_{aa}": seq.count(aa) / n for aa in AA}


def feat_dipeptide(seq: str) -> dict:
    """400 dipeptide frequencies (dp_XY). Confirmed top ML features."""
    n = max(len(seq) - 1, 1)
    counts: dict = {f"dp_{a}{b}": 0 for a in AA for b in AA}
    for i in range(len(seq) - 1):
        key = f"dp_{seq[i]}{seq[i+1]}"
        if key in counts:
            counts[key] += 1
    return {k: v / n for k, v in counts.items()}


def feat_pseudo_aac(seq: str, lam: int = 10, w: float = 0.05) -> dict:
    """PseAAC — pseudo amino acid composition (30 features)."""
    # Kyte-Doolittle hydrophobicity + hydrophilicity
    H  = dict(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-0.4,H=-3.2,
              I=4.5,K=-3.9,L=3.8,M=1.9,N=-3.5,P=-1.6,Q=-3.5,
              R=-4.5,S=-0.8,T=-0.7,V=4.2,W=-0.9,Y=-1.3)
    Hy = dict(A=-0.5,C=-1.0,D=3.0,E=3.0,F=-2.5,G=0.0,H=-0.5,
              I=-1.8,K=3.0,L=-1.8,M=-1.3,N=0.2,P=0.0,Q=0.2,
              R=3.0,S=0.3,T=-0.4,V=-1.5,W=-3.4,Y=-2.3)
    n = len(seq)
    feats = {f"pse_{aa}": seq.count(aa) / max(n, 1) for aa in AA}
    for lam_i in range(1, min(lam + 1, n)):
        theta = sum(
            (H.get(seq[i], 0) - H.get(seq[i + lam_i], 0)) ** 2 +
            (Hy.get(seq[i], 0) - Hy.get(seq[i + lam_i], 0)) ** 2
            for i in range(n - lam_i)
        ) / (2 * (n - lam_i))
        feats[f"pse_theta_{lam_i}"] = w * theta
    return feats


def feat_physicochemical(seq: str) -> dict:
    """MW, pI, GRAVY, aromaticity, instability, fractions."""
    n = max(len(seq), 1)
    feats = {
        "phys_length":        len(seq),
        "phys_length_log":    math.log1p(len(seq)),
        "phys_mw":            0.0,
        "phys_pi":            7.0,
        "phys_charge_ph7":    0.0,
        "phys_gravy":         0.0,
        "phys_aromaticity":   0.0,
        "phys_instability":   0.0,
        "phys_frac_glycine":  seq.count("G") / n,
        "phys_frac_proline":  seq.count("P") / n,
        "phys_frac_charged":  sum(seq.count(a) for a in "DEKR") / n,
        "phys_frac_aromatic": sum(seq.count(a) for a in "FHWY") / n,
        "phys_frac_polar":    sum(seq.count(a) for a in "NQSTY") / n,
        "phys_frac_nonpolar": sum(seq.count(a) for a in "AVILMFWP") / n,
        "phys_frac_small":    sum(seq.count(a) for a in "ACGST") / n,
    }
    if HAS_BIOPYTHON and len(seq) >= 10:
        try:
            clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
            if clean:
                pa = ProteinAnalysis(clean)
                feats["phys_mw"]          = pa.molecular_weight()
                feats["phys_pi"]          = pa.isoelectric_point()
                feats["phys_charge_ph7"]  = pa.charge_at_pH(7.0)
                feats["phys_gravy"]       = pa.gravy()
                feats["phys_aromaticity"] = pa.aromaticity()
                feats["phys_instability"] = pa.instability_index()
        except Exception:
            pass
    return feats


def feat_catalytic_core(seq: str) -> dict:
    """
    Catalytic-core features — middle 50% of protein (confirmed in working code).
    inv_cat_* prefix. Top catalytic feature: catalytic_K_percent (4.3% importance).
    """
    n = len(seq)
    start, end = n // 4, 3 * n // 4
    core = seq[start:end] if n >= 20 else seq

    feats: dict = {}
    for aa in "DEHKCST":
        feats[f"inv_cat_{aa}"] = core.count(aa) / max(len(core), 1)

    # Catalytic residue clustering (confirmed: mean/std/min/max dist)
    cat_pos = [i for i, a in enumerate(seq) if a in "DEHKC"]
    if len(cat_pos) > 1:
        dists = [cat_pos[i+1] - cat_pos[i] for i in range(len(cat_pos) - 1)]
        feats["inv_cat_mean_dist"]   = float(np.mean(dists))
        feats["inv_cat_std_dist"]    = float(np.std(dists))
        feats["inv_cat_min_dist"]    = float(np.min(dists))
        feats["inv_cat_max_dist"]    = float(np.max(dists))
        feats["inv_cat_clustering"]  = 1.0 / (np.mean(dists) + 1)
    else:
        for k in ["mean_dist","std_dist","min_dist","max_dist","clustering"]:
            feats[f"inv_cat_{k}"] = 0.0

    # Global physicochemical balance (confirmed in working code)
    feats["inv_hydrophobic"] = sum(seq.count(a) for a in "AILMFVWY") / max(n, 1)
    feats["inv_charged"]     = sum(seq.count(a) for a in "DEKR")     / max(n, 1)
    feats["inv_polar"]       = sum(seq.count(a) for a in "STNQ")     / max(n, 1)
    feats["inv_aromatic"]    = sum(seq.count(a) for a in "FWY")      / max(n, 1)
    feats["inv_net_charge"]  = (seq.count("K") + seq.count("R")
                                - seq.count("D") - seq.count("E")) / max(n, 1)
    return feats


def feat_ec_motifs(seq: str) -> dict:
    """
    7 EC-specific sequence motifs (confirmed in CarboxylaseDatabase, Dec 2025).
    Normalized per 100 residues.
    """
    n = max(len(seq), 1)
    return {
        "motif_rubisco_kk":    seq.count("KK")  / n * 100,
        "motif_rubisco_gk":    seq.count("GK")  / n * 100,
        "motif_ca_hh":         seq.count("HH")  / n * 100,
        "motif_ca_his_cluster": sum(1 for i in range(n - 3)
                                    if seq[i:i+4].count("H") >= 2) / n * 100,
        "motif_pepc_rr":       seq.count("RR")  / n * 100,
        "motif_biotin_mk":     seq.count("MK")  / n * 100,
        "motif_biotin_amk":    seq.count("AMK") / n * 100,
    }


# PROSITE patterns — confirmed from working v2 code
PROSITE_PATTERNS = {
    "PS00157_RUBISCO_LARGE":  r"[GA].[LIVM]K[GP][HY][LI]",
    "PS00158_RUBISCO_SMALL":  r"E.{5}[PS]W[KR]L.{0,1}[FYQM]",
    "PS00162_CARB_ANHYDRASE": r"C.{2}DS.{1,3}[AP]",
    "PS00188_BIOTIN_LIPOYL":  r"[LIVMF][LIMSTAC].[LIVMFSTAC].[SA]MK[MQ].{3}[LIVMFCT]",
    "PS00781_PEPC_1":         r"[LIVMF].[LIVA]K[LIVM]HG[DA]",
    "PS00393_PEPC_2":         r"G[LIVM]RCG[AP]E[PQ]",
    "PS00017_ATP_GTP_A":      r"[AG].{4}GK[ST]",
    "PS00021_CYS_ACTIVE_SITE":r"[LIVMFYWC].{1,2}[LIVMFYW].C.[SACLIVMFYW]",
    "PS00110_BIOTIN_CARBOX":  r"[LIVMFA].{1,2}[IVLM].AM.K.[LIVMFY]",
    "PS00013_ASP_PROTEASE":   r"[LIVMFGAC].{0,1}[IVFYW].{2}[LIVMFGAC].D.G.[LIVMFGAC]",
    "PS00124_PEPCK":          r"[ST][GAS].{1}[LIVMF].{2,4}[LIVMF].{2,4}[LIVM].{1,3}[DN]",
    "PS00519_RUBISCO_ACTSITE":r"D.{6}[EG].{4}[LIVMFY].{3}[LIVMFYA].{2}[KR]",
    "PS00521_LGS_SYNTHASE":   r"[LIVMF].{2}[DE].{3}[LIVMF].G.G.{2}[ST]",
    "PS01213_PEPC_FAMILY":    r"[LIVM].{2}[KR].{2}[DE].{5,7}[LIVM].{3}[LIVM]",
}


def feat_prosite(seq: str) -> dict:
    """14 PROSITE pattern features (count + binary)."""
    feats = {}
    for name, pattern in PROSITE_PATTERNS.items():
        try:
            matches = list(re.finditer(pattern, seq))
            feats[f"prosite_{name}_count"]   = len(matches)
            feats[f"prosite_{name}_present"] = 1 if matches else 0
        except re.error:
            feats[f"prosite_{name}_count"]   = 0
            feats[f"prosite_{name}_present"] = 0
    return feats


def compute_layer_A(seq: str) -> dict:
    seq = seq.upper()
    return {
        **feat_aa_composition(seq),
        **feat_dipeptide(seq),
        **feat_pseudo_aac(seq),
        **feat_physicochemical(seq),
        **feat_catalytic_core(seq),
        **feat_ec_motifs(seq),
        **feat_prosite(seq),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER B — Pfam/HMMER (BATCH mode — run once on full FASTA)
# ═══════════════════════════════════════════════════════════════════════════════

# Key Pfam domains for carboxylases (confirmed from working DB)
PFAM_BINARY = [
    "PF00016", "PF02788", "PF00101",  # RuBisCO large + small
    "PF00194", "PF03119",              # Carbonic anhydrase α + β
    "PF00311",                         # PEPC barrel
    "PF00821",                         # PEPCK
    "PF02785", "PF00364", "PF01039",  # Biotin carboxylase family
    "PF02786", "PF02787",              # CPS large chain
    "PF00289",                         # CPS large chain 2
    "PF01309",                         # VKGC
    "PF03599", "PF03590",              # CODH
    "PF00384",                         # FDH
    "PF00682",                         # HMGL-like
    "PF00101",                         # RuBisCO small
]


def run_hmmscan_batch(fasta_path: Path, n_cpu: int = 8) -> dict:
    """
    Run hmmscan ONCE on the full FASTA. Returns {uid: {pfam_acc, ...}}.
    This is the correct batch approach — NOT per-sequence subprocess calls.
    """
    if not PFAM_HMM.exists():
        print(f"  SKIP Layer B: Pfam-A.hmm not found at {PFAM_HMM}")
        print(f"  Download: wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz")
        print(f"  Then: hmmpress data/raw/pfam/Pfam-A.hmm")
        return {}

    domtbl = fasta_path.with_suffix(".domtbl")
    cmd = [
        "hmmscan", "--domtblout", str(domtbl),
        "--cpu", str(n_cpu),
        "-E", "1e-3",
        str(PFAM_HMM), str(fasta_path),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  SKIP Layer B: hmmscan failed — {e}")
        return {}

    # Parse domtblout
    hits: dict = {}
    with open(domtbl) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            target  = parts[0]           # sequence name
            pfam_ac = parts[4].split(".")[0]  # Pfam accession without version
            evalue  = float(parts[12])
            if evalue <= 1e-3:
                hits.setdefault(target, set()).add(pfam_ac)
    return hits


def pfam_features(uid: str, hits: dict) -> dict:
    uid_hits = hits.get(uid, set())
    feats = {f"pfam_{acc}": int(acc in uid_hits) for acc in PFAM_BINARY}
    feats["pfam_n_hits"] = len(uid_hits)
    feats["pfam_hits_json"] = json.dumps(sorted(uid_hits))
    return feats


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER C — BLAST (BATCH mode)
# ═══════════════════════════════════════════════════════════════════════════════

def build_blast_db(brenda_fasta: Path):
    """Build BLAST database from BRENDA positive FASTA (run once)."""
    BLAST_DB.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["makeblastdb", "-in", str(brenda_fasta),
           "-dbtype", "prot", "-out", str(BLAST_DB)]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"  BLAST DB created: {BLAST_DB}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  WARNING: makeblastdb failed — {e}")


def run_blastp_batch(fasta_path: Path, n_cpu: int = 8) -> dict:
    """
    Run blastp ONCE on full FASTA in batch mode.
    Returns {uid: {pident, evalue, best_ec, has_hit}}.
    """
    if not Path(str(BLAST_DB) + ".psq").exists() and \
       not Path(str(BLAST_DB) + ".pdb").exists():
        print(f"  SKIP Layer C: BLAST DB not found at {BLAST_DB}")
        print(f"  Build with: makeblastdb -in <brenda_positives.fasta> -dbtype prot -out {BLAST_DB}")
        return {}

    out_file = fasta_path.with_suffix(".blast.tsv")
    cmd = [
        "blastp",
        "-query",   str(fasta_path),
        "-db",      str(BLAST_DB),
        "-out",     str(out_file),
        "-outfmt",  "6 qseqid sseqid pident evalue stitle",
        "-max_target_seqs", "1",
        "-evalue",  "1e-3",
        "-num_threads", str(n_cpu),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=3600)
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired) as e:
        print(f"  SKIP Layer C: blastp failed — {e}")
        return {}

    hits: dict = {}
    with open(out_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            uid    = parts[0].split("|")[0]
            pident = float(parts[2])
            evalue = float(parts[3])
            title  = parts[4] if len(parts) > 4 else ""
            # Only keep best hit per query (first occurrence in tabular output)
            if uid not in hits:
                m = re.search(r"\|EC:([0-9.]+)", title)
                best_ec = m.group(1) if m else ""
                hits[uid] = {
                    "blast_best_pident": pident,
                    "blast_best_evalue": evalue,
                    "blast_best_ec":     best_ec,
                    "blast_has_hit":     1,
                }
    return hits


def blast_features(uid: str, hits: dict) -> dict:
    if uid in hits:
        return hits[uid]
    return {"blast_best_pident": 0.0, "blast_best_evalue": 999.0,
            "blast_best_ec": "", "blast_has_hit": 0}


# ═══════════════════════════════════════════════════════════════════════════════
# LAYER D — MEME motif hits (from pre-computed meme_hits.tsv)
# ═══════════════════════════════════════════════════════════════════════════════

def load_meme_hits() -> dict:
    """
    Load MEME FIMO hit table if it exists.
    Returns {uid: {meme_col: 0/1, ...}}.

    Expected format of data/features/meme/meme_hits.tsv:
        uniprot_id  meme_RuBisCO_001_HPGYGFL  meme_ACC_002_KPKLGL  ...
        P00875      1                           0                    ...
    """
    if not MEME_HITS.exists():
        return {}
    df = pd.read_csv(MEME_HITS, sep="\t", dtype=str)
    df = df.fillna("0")
    meme_cols = [c for c in df.columns if c != "uniprot_id"]
    result = {}
    for _, row in df.iterrows():
        uid = str(row["uniprot_id"])
        result[uid] = {col: int(str(row[col])) for col in meme_cols}
    print(f"  MEME hits loaded: {len(result):,} sequences, {len(meme_cols)} motifs")
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Main orchestrator
# ═══════════════════════════════════════════════════════════════════════════════

def annotate(input_tsv: Path, output_tsv: Path, layers: set, n_cpu: int = 8):

    df = pd.read_csv(input_tsv, sep="\t", low_memory=False)
    df["sequence"] = df["sequence"].fillna("").astype(str)
    print(f"Loaded {len(df):,} sequences from {input_tsv}")

    # ── Pre-compute batch annotations ─────────────────────────────────────────
    pfam_hits  = {}
    blast_hits = {}
    meme_hits  = {}

    if "B" in layers:
        print("\nLayer B — Writing temp FASTA for hmmscan...")
        tmp = Path(tempfile.mktemp(suffix=".fasta"))
        _write_temp_fasta(df, tmp)
        pfam_hits = run_hmmscan_batch(tmp, n_cpu=n_cpu)
        tmp.unlink(missing_ok=True)
        print(f"  Pfam hits for {len(pfam_hits):,} sequences")

    if "C" in layers:
        print("\nLayer C — Writing temp FASTA for blastp...")
        tmp = Path(tempfile.mktemp(suffix=".fasta"))
        _write_temp_fasta(df, tmp)
        blast_hits = run_blastp_batch(tmp, n_cpu=n_cpu)
        tmp.unlink(missing_ok=True)
        print(f"  BLAST hits for {len(blast_hits):,} sequences")

    if "D" in layers:
        print("\nLayer D — Loading MEME motif hits...")
        meme_hits = load_meme_hits()
        if not meme_hits:
            print(f"  SKIP: {MEME_HITS} not found (MEME subproject pending)")

    # ── Per-sequence annotation ────────────────────────────────────────────────
    records = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Annotating"):
        uid = str(row.get("uniprot_id", ""))
        seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(row.get("sequence", "")).upper())

        rec = {
            "uniprot_id":    uid,
            "ec_number":     str(row.get("ec_number",     "")),
            "label":         row.get("label",             -1),
            "evidence_tier": row.get("evidence_tier",      3),
            "km_best_mM":    row.get("km_best_mM",        None),
            "km_log10_mM":   row.get("km_log10_mM",       None),
            "sequence":      str(row.get("sequence",       "")),
            "length":        len(seq),
            "organism":      str(row.get("organism",       "")),
            "source":        str(row.get("source",         "")),
        }

        if len(seq) >= 10:
            if "A" in layers:
                rec.update(compute_layer_A(seq))
            if "B" in layers:
                rec.update(pfam_features(uid, pfam_hits))
            if "C" in layers:
                rec.update(blast_features(uid, blast_hits))
            if "D" in layers and uid in meme_hits:
                rec.update(meme_hits[uid])

        records.append(rec)

    df_out = pd.DataFrame(records)
    df_out.to_csv(output_tsv, sep="\t", index=False)

    n_feat = len(df_out.columns) - 10  # subtract metadata cols
    print(f"\n✓ Feature matrix: {len(df_out):,} sequences × {n_feat:,} features")
    print(f"✓ Saved → {output_tsv}")

    # Feature category summary
    cols = df_out.columns.tolist()
    print(f"\n  Feature categories:")
    print(f"    AA composition:     {len([c for c in cols if c.startswith('aa_')]):>5}")
    print(f"    Dipeptides:         {len([c for c in cols if c.startswith('dp_')]):>5}  ← top ML features")
    print(f"    Pseudo-AAC:         {len([c for c in cols if c.startswith('pse_')]):>5}")
    print(f"    Physicochemical:    {len([c for c in cols if c.startswith('phys_')]):>5}")
    print(f"    Catalytic core:     {len([c for c in cols if c.startswith('inv_')]):>5}")
    print(f"    EC-specific motifs: {len([c for c in cols if c.startswith('motif_')]):>5}")
    print(f"    PROSITE patterns:   {len([c for c in cols if c.startswith('prosite_')]):>5}")
    print(f"    Pfam domains:       {len([c for c in cols if c.startswith('pfam_')]):>5}")
    print(f"    BLAST homology:     {len([c for c in cols if c.startswith('blast_')]):>5}")
    print(f"    MEME motifs:        {len([c for c in cols if c.startswith('meme_')]):>5}  (pending if 0)")


def _write_temp_fasta(df: pd.DataFrame, path: Path):
    with open(path, "w") as f:
        for _, row in df.iterrows():
            seq = str(row.get("sequence", "") or "").strip()
            uid = str(row.get("uniprot_id", ""))
            if seq and uid:
                f.write(f">{uid}\n{seq}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input",  default="data/interim/master_combined.tsv")
    p.add_argument("--output", default=f"data/processed/features_{datetime.now().strftime('%Y%m%d')}.tsv")
    p.add_argument("--layers", nargs="+", default=["A", "B", "C", "D"],
                   choices=["A", "B", "C", "D"])
    p.add_argument("--n-cpu",  type=int, default=8)
    args = p.parse_args()

    annotate(
        input_tsv  = Path(args.input),
        output_tsv = Path(args.output),
        layers     = set(args.layers),
        n_cpu      = args.n_cpu,
    )


if __name__ == "__main__":
    main()
