#!/usr/bin/env python3
"""
Filter the pilot 1M FASTA down to sequences with at least one hit to the
18 carboxylase-related Pfam domains that v5 was trained on.

Reads:  data/gmgc/interim/pfam_hits/pilot_pfam_hits.tbl
        data/gmgc/interim/pilot_1M.faa
Writes: data/gmgc/interim/candidates.faa  (sequences with at least one carboxylase Pfam hit)
        data/gmgc/interim/candidates_pfam.tsv (per-sequence Pfam hit table for downstream features)
"""
import sys
from pathlib import Path
from collections import defaultdict

ROOT = Path.home() / "Projects_shared/CarboDB_v3"

# The 18 Pfam IDs v5 uses (from data/ml/feature_names_binary.json)
CARBOXY_PFAMS = {
    "PF00016",  # RuBisCO_large
    "PF00101",  # RuBisCO_small
    "PF00194",  # Carbonic_anhydrase
    "PF00289",  # Biotin_carb_N
    "PF00311",  # PEP_carboxylase
    "PF00364",  # Biotin_lipoyl
    "PF00384",  # Molybdopterin
    "PF00682",  # HMGL-like
    "PF00821",  # PEPCK_C
    "PF01039",  # Carboxyl_trans
    "PF01309",  # C5_kinase
    "PF02785",  # Biotin_carb_M
    "PF02786",  # Biotin_carb_C
    "PF02787",  # Biotin_carb_central
    "PF02788",  # RuBisCO_large_N
    "PF03119",  # PEPCK_N
    "PF03590",  # ATP_grasp_misc
    "PF03599",  # CarbamoylP_synth_lrg
}

HITS_FILE = ROOT / "data/gmgc/interim/pfam_hits/pilot_pfam_hits.tbl"
FASTA_IN  = ROOT / "data/gmgc/interim/pilot_1M.faa"
FASTA_OUT = ROOT / "data/gmgc/interim/candidates.faa"
PFAM_TSV  = ROOT / "data/gmgc/interim/candidates_pfam.tsv"

# ── 1. PARSE PFAM HITS ────────────────────────────────────────────────
# domtblout format (whitespace-separated):
# 1 target_name   2 acc          3 tlen   4 query_name  5 qacc       6 qlen   7 fE   8 fSc  ...
#                 (Pfam acc with version, e.g. PF00016.27)

print(f"Reading Pfam hits from {HITS_FILE}", file=sys.stderr)
seq_carboxy_hits = defaultdict(set)   # seq_id -> set of Pfam accs hit
seq_all_pfam = defaultdict(list)       # seq_id -> list of (acc_no_version, evalue, bitscore)
n_lines = 0
n_carboxy = 0
with open(HITS_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) < 8:
            continue
        target_name = parts[0]   # e.g., RuBisCO_large
        acc_full    = parts[1]   # e.g., PF00016.27
        query_name  = parts[3]   # e.g., GMGC10.008_717_091.UNKNOWN_1
        evalue      = float(parts[6])
        bitscore    = float(parts[7])

        # Strip version suffix to compare to CARBOXY_PFAMS
        acc = acc_full.split(".")[0] if "." in acc_full else acc_full
        seq_all_pfam[query_name].append((acc, evalue, bitscore))

        if acc in CARBOXY_PFAMS:
            seq_carboxy_hits[query_name].add(acc)
            n_carboxy += 1
        n_lines += 1

print(f"  Parsed {n_lines:,} Pfam hit lines", file=sys.stderr)
print(f"  {n_carboxy:,} hits to carboxylase Pfams", file=sys.stderr)
print(f"  {len(seq_carboxy_hits):,} sequences with at least one carboxylase Pfam hit", file=sys.stderr)

# ── 2. WRITE Pfam TSV (for downstream v5 features) ───────────────────
print(f"Writing per-sequence Pfam table to {PFAM_TSV}", file=sys.stderr)
with open(PFAM_TSV, "w") as f:
    f.write("seq_id\tpfam_acc\tevalue\tbitscore\n")
    for seq_id in seq_carboxy_hits:
        for acc, ev, bs in seq_all_pfam[seq_id]:
            f.write(f"{seq_id}\t{acc}\t{ev}\t{bs}\n")

# ── 3. STREAM-FILTER THE FASTA ────────────────────────────────────────
print(f"Streaming {FASTA_IN} -> {FASTA_OUT}", file=sys.stderr)
candidate_ids = set(seq_carboxy_hits.keys())
n_written = 0
n_seen = 0
with open(FASTA_IN) as fin, open(FASTA_OUT, "w") as fout:
    keep = False
    for line in fin:
        if line.startswith(">"):
            n_seen += 1
            sid = line[1:].strip().split()[0]
            keep = sid in candidate_ids
            if keep:
                n_written += 1
        if keep:
            fout.write(line)

print(f"  Scanned: {n_seen:,}", file=sys.stderr)
print(f"  Written: {n_written:,} candidate sequences", file=sys.stderr)
print(f"Done.", file=sys.stderr)
