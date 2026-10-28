#!/usr/bin/env python3
"""
Filter the FULL marine FASTA down to sequences hitting any of the 18
carboxylase Pfams, using hmmsearch domtblout output.

NOTE: hmmsearch output has swapped target/query columns vs hmmscan:
  hmmsearch:  col 1 = sequence_id, col 5 = Pfam acc
  hmmscan:    col 1 = Pfam name,   col 2 = Pfam acc, col 4 = sequence_id

Reads:  data/gmgc/interim_full/marine_full_pfam_hits.domtbl
        data/gmgc/interim_full/marine.95nr.faa
Writes: data/gmgc/interim_full/candidates_full.faa
        data/gmgc/interim_full/candidates_full_pfam.tsv
"""
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path.home() / "Projects_shared/CarboDB_v3"

# Same 18 v5 carboxylase Pfams as script 04
CARBOXY_PFAMS = {
    "PF00016", "PF00101", "PF00194", "PF00289", "PF00311", "PF00364",
    "PF00384", "PF00682", "PF00821", "PF01039", "PF01309", "PF02785",
    "PF02786", "PF02787", "PF02788", "PF03119", "PF03590", "PF03599",
}

HITS_FILE = ROOT / "data/gmgc/interim_full/marine_full_pfam_hits.domtbl"
FASTA_IN  = ROOT / "data/gmgc/interim_full/marine.95nr.faa"
FASTA_OUT = ROOT / "data/gmgc/interim_full/candidates_full.faa"
PFAM_TSV  = ROOT / "data/gmgc/interim_full/candidates_full_pfam.tsv"

# ── 1. PARSE hmmsearch DOMTBLOUT ─────────────────────────────────────
# hmmsearch columns (whitespace-separated):
#   1 target_name (= sequence ID)
#   2 target_acc  (= sequence acc, usually '-')
#   3 tlen        (sequence length)
#   4 query_name  (HMM name, e.g. RuBisCO_large)
#   5 query_acc   (Pfam acc, e.g. PF00016.27)
#   6 qlen        (HMM length)
#   7 fE          (full-sequence E-value)
#   8 fSc         (full-sequence bitscore)
print(f"Reading hmmsearch hits from {HITS_FILE}", file=sys.stderr)
seq_carboxy_hits = defaultdict(set)
seq_all_pfam = defaultdict(list)
n_lines = 0
n_carboxy = 0
with open(HITS_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split()
        if len(parts) < 8:
            continue
        seq_id     = parts[0]
        acc_full   = parts[4]
        evalue     = float(parts[6])
        bitscore   = float(parts[7])

        # Strip version suffix from Pfam acc
        acc = acc_full.split(".")[0] if "." in acc_full else acc_full
        seq_all_pfam[seq_id].append((acc, evalue, bitscore))
        if acc in CARBOXY_PFAMS:
            seq_carboxy_hits[seq_id].add(acc)
            n_carboxy += 1
        n_lines += 1

print(f"  Parsed {n_lines:,} domtblout lines", file=sys.stderr)
print(f"  Carboxylase Pfam hits: {n_carboxy:,}", file=sys.stderr)
print(f"  Unique sequences hitting >=1 carboxylase Pfam: {len(seq_carboxy_hits):,}",
      file=sys.stderr)

# ── 2. WRITE Pfam TSV (input for v5 feature builder) ─────────────────
print(f"\nWriting per-sequence Pfam table to {PFAM_TSV}", file=sys.stderr)
with open(PFAM_TSV, "w") as f:
    f.write("seq_id\tpfam_acc\tevalue\tbitscore\n")
    for sid in seq_carboxy_hits:
        for acc, ev, bs in seq_all_pfam[sid]:
            f.write(f"{sid}\t{acc}\t{ev}\t{bs}\n")

# ── 3. STREAM-FILTER THE FASTA ───────────────────────────────────────
print(f"\nStreaming {FASTA_IN} -> {FASTA_OUT}", file=sys.stderr)
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
            if n_seen % 5_000_000 == 0:
                print(f"  {n_seen:,} scanned, {n_written:,} kept", file=sys.stderr)
        if keep:
            fout.write(line)

print(f"\nFinal counts:", file=sys.stderr)
print(f"  Scanned: {n_seen:,}", file=sys.stderr)
print(f"  Written: {n_written:,} carboxylase candidates", file=sys.stderr)
print(f"\nDone.", file=sys.stderr)
