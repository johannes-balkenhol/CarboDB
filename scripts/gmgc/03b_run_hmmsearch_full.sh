#!/bin/bash
# Run hmmsearch on the FULL translated marine FASTA (88M sequences)
# against the 18-Pfam restricted library.
#
# hmmsearch is fundamentally faster than hmmscan when queries >> HMMs:
# - hmmscan: each sequence scanned against all HMMs (slow with many seqs)
# - hmmsearch: each HMM scanned against all sequences (fast)
# Combined with reduced HMM count (18 vs 21k), expect ~50-100x speedup.
set -euo pipefail

cd ~/Projects_shared/CarboDB_v3

INPUT="data/gmgc/interim/marine.95nr.faa.gz"
PFAM_DB="data/dbs/pfam/carboxy_pfams.hmm"
WORK_DIR="data/gmgc/interim_full"
OUT_TBL="$WORK_DIR/marine_full_pfam_hits.domtbl"
OUT_LOG="$WORK_DIR/hmmsearch.log"

mkdir -p "$WORK_DIR"

N_CPUS=32   # heavy parallelization (we have 128 cores; reserve some)

echo "Started: $(date)"
echo
echo "Input:    $INPUT ($(ls -lh $INPUT | awk '{print $5}'))"
echo "HMM DB:   $PFAM_DB ($(ls -lh $PFAM_DB | awk '{print $5}'), $(grep -c '^NAME ' $PFAM_DB) HMMs)"
echo "Output:   $OUT_TBL"
echo "CPUs:     $N_CPUS"
echo

# Need uncompressed FASTA for hmmsearch (it doesn't read .gz directly via stdin reliably)
UNZIPPED="$WORK_DIR/marine.95nr.faa"
if [ ! -f "$UNZIPPED" ]; then
    echo "Decompressing FASTA..."
    pigz -dc "$INPUT" > "$UNZIPPED" || gunzip -c "$INPUT" > "$UNZIPPED"
    echo "Decompressed: $(ls -lh $UNZIPPED | awk '{print $5}')"
fi

# Run hmmsearch
# Each HMM searches all sequences; --cpu sets threads per HMM internally
# With 18 HMMs, we let hmmsearch use multiple threads per HMM
echo
echo "Running hmmsearch..."
time hmmsearch \
    --domtblout "$OUT_TBL" \
    -E 1e-3 \
    --cpu $N_CPUS \
    --noali \
    "$PFAM_DB" \
    "$UNZIPPED" \
    > "$OUT_LOG" 2>&1

echo
echo "=== Result ==="
echo "Domain table: $(ls -lh $OUT_TBL | awk '{print $5}')"
echo "Lines:        $(wc -l < $OUT_TBL)"
echo "Hit count:    $(grep -vc '^#' $OUT_TBL)"
echo
echo "Top HMMs by hit count:"
grep -v '^#' "$OUT_TBL" | awk '{print $4}' | sort | uniq -c | sort -rn | head

echo
echo "Finished: $(date)"
