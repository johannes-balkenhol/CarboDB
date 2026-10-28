#!/bin/bash
# Translate the marine nucleotide FASTA to amino acids using transeq (EMBOSS).
# Frame 1 only — these are predicted ORFs from GMGC, so frame is canonical.
set -euo pipefail
cd ~/Projects_shared/CarboDB_v3/data/gmgc

NT_GZ="downloads/GMGC10.marine.95nr.fna.gz"
AA_OUT="interim/marine.95nr.faa"
AA_GZ="interim/marine.95nr.faa.gz"

echo "Started: $(date)"
echo "Input:  $NT_GZ ($(ls -lh $NT_GZ | awk '{print $5}'))"

# Decompress on the fly, translate frame 1
# -clean replaces internal stop * with X
# -trim removes terminal stops
zcat "$NT_GZ" | transeq \
    -filter \
    -frame 1 \
    -clean \
    -trim \
    > "$AA_OUT" 2> interim/transeq.log

echo "Translation done: $(date)"
echo "Output: $(ls -lh $AA_OUT | awk '{print $5}')"
echo "Sequences: $(grep -c '^>' $AA_OUT)"

# Gzip for storage
gzip "$AA_OUT"
echo "Compressed: $(ls -lh $AA_GZ | awk '{print $5}')"
echo "Finished: $(date)"
