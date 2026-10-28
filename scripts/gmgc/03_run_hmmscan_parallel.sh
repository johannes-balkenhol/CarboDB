#!/bin/bash
# Run hmmscan in parallel across 16 chunks of the pilot 1M FASTA.
# Uses 2 CPUs per chunk = 32 cores total. Output: 16 domtblout TSVs.
set -euo pipefail

cd ~/Projects_shared/CarboDB_v3

INPUT="data/gmgc/interim/pilot_1M.faa"
PFAM_DB="data/dbs/pfam/Pfam-A.hmm"
CHUNK_DIR="data/gmgc/interim/chunks"
HIT_DIR="data/gmgc/interim/pfam_hits"
N_CHUNKS=16
CPUS_PER_CHUNK=2

mkdir -p "$CHUNK_DIR" "$HIT_DIR"

echo "Started: $(date)"
echo "Splitting $INPUT into $N_CHUNKS chunks..."

# Split: every Nth sequence into chunk_<i>.faa
python3 -c "
import os
from pathlib import Path
n_chunks = $N_CHUNKS
chunk_dir = Path('$CHUNK_DIR')
inputs = open('$INPUT')
outputs = [open(chunk_dir / f'chunk_{i:02d}.faa', 'w') for i in range(n_chunks)]
i = -1
for line in inputs:
    if line.startswith('>'):
        i += 1
    outputs[i % n_chunks].write(line)
for o in outputs:
    o.close()
print('Chunks written.')
"

echo "Chunks:"
ls -lh "$CHUNK_DIR"/chunk_*.faa | head -5

echo
echo "Running hmmscan in parallel ($N_CHUNKS chunks x $CPUS_PER_CHUNK CPUs each)..."

# Launch all chunks in background, wait for all
pids=()
for i in $(seq -f "%02g" 0 $((N_CHUNKS-1))); do
  (
    hmmscan --domtblout "$HIT_DIR/chunk_${i}.tbl" \
            -E 1e-3 \
            --cpu $CPUS_PER_CHUNK \
            --noali \
            "$PFAM_DB" \
            "$CHUNK_DIR/chunk_${i}.faa" \
            > "$HIT_DIR/chunk_${i}.log" 2>&1
    echo "Chunk $i done at $(date)"
  ) &
  pids+=($!)
done

# Wait for all to finish
echo "Launched ${#pids[@]} chunks (PIDs: ${pids[*]})"
echo "Waiting for completion..."
for pid in "${pids[@]}"; do
  wait "$pid"
done

echo
echo "All chunks done: $(date)"
echo "Hit tables:"
ls -lh "$HIT_DIR"/chunk_*.tbl

# Concatenate all hit tables (each has its own header, the parser strips comments)
cat "$HIT_DIR"/chunk_*.tbl > "$HIT_DIR/pilot_pfam_hits.tbl"
echo
echo "Combined: $(wc -l < $HIT_DIR/pilot_pfam_hits.tbl) lines total"
echo "Finished: $(date)"
