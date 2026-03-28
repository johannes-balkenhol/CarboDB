#!/bin/bash
#SBATCH --job-name=carbodb_ankh
#SBATCH --chdir=/storage/users/job37yv/Projects/CarboDB_v3
#SBATCH --partition=hades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/ankh_%A_%a.log
#SBATCH --error=logs/ankh_%A_%a.log
#SBATCH --array=1-2380%30

# CarboxyDB — Ankh protein language model embeddings (CPU)
# Model: ankh-large (1024-dim, CPU-runnable)
#
# Submit:
#   sbatch scripts/04d_ankh.sh
#
# Monitor:
#   ls data/features/ankh/*.npy | wc -l

set -euo pipefail

PROJECT=/storage/users/job37yv/Projects/CarboDB_v3
CHUNKS=${PROJECT}/data/interim/fasta_chunks
OUT_DIR=${PROJECT}/data/features/ankh

mkdir -p ${OUT_DIR}

CHUNK_OFFSET=${CHUNK_OFFSET:-0}
CHUNK=$(printf "%04d" $((${SLURM_ARRAY_TASK_ID} + ${CHUNK_OFFSET})))
INPUT=${CHUNKS}/chunk_${CHUNK}.fasta
OUTPUT=${OUT_DIR}/ankh_${CHUNK}.tsv

if [ -f "${OUTPUT}" ]; then
    echo "Already done: ${OUTPUT}"
    exit 0
fi

if [ ! -f "${INPUT}" ]; then
    echo "Input not found: ${INPUT}"
    exit 1
fi

echo "=== Ankh chunk ${CHUNK} ==="
date

python3 << 'PYEOF'
import os, sys
import numpy as np
import pandas as pd
from pathlib import Path

# Limit CPU threads
os.environ["OMP_NUM_THREADS"]     = "8"
os.environ["MKL_NUM_THREADS"]     = "8"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

import torch
torch.set_num_threads(8)

import ankh

chunk  = "${CHUNK}"
input_path  = Path("${INPUT}")
output_path = Path("${OUTPUT}")

# Load model once
print("Loading Ankh-large model...")
model, tokenizer = ankh.load_large_model()
model.eval()
print("Model loaded.")

# Read FASTA
seqs = []
cdb_ids = []
with open(input_path) as f:
    cdb_id = seq_lines = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if cdb_id and seq_lines:
                seqs.append("".join(seq_lines))
                cdb_ids.append(cdb_id)
            cdb_id   = line[1:].split("|")[0]
            seq_lines = []
        elif line:
            if seq_lines is not None:
                seq_lines.append(line)
    if cdb_id and seq_lines:
        seqs.append("".join(seq_lines))
        cdb_ids.append(cdb_id)

print(f"Processing {len(seqs)} sequences...")

BATCH_SIZE = 8   # conservative for CPU RAM
rows = []

for i in range(0, len(seqs), BATCH_SIZE):
    batch_seqs    = seqs[i:i+BATCH_SIZE]
    batch_ids     = cdb_ids[i:i+BATCH_SIZE]

    # Truncate very long sequences (Ankh has token limit)
    batch_seqs_trunc = [s[:1022] for s in batch_seqs]

    # Tokenise
    ids = tokenizer.batch_encode_plus(
        [list(s) for s in batch_seqs_trunc],
        add_special_tokens=True,
        padding=True,
        is_split_into_words=True,
        return_tensors="pt"
    )

    with torch.no_grad():
        output = model(
            input_ids=ids["input_ids"],
            attention_mask=ids["attention_mask"]
        )

    # Mean-pool over sequence length (exclude padding)
    hidden = output.last_hidden_state  # (B, L, 1024)
    mask   = ids["attention_mask"].unsqueeze(-1).float()
    embs   = (hidden * mask).sum(dim=1) / mask.sum(dim=1)
    embs   = embs.cpu().numpy()

    for cdb_id, emb in zip(batch_ids, embs):
        row = {"cdb_id": cdb_id}
        for j, v in enumerate(emb):
            row[f"ankh_{j}"] = float(v)
        rows.append(row)

    if (i // BATCH_SIZE) % 10 == 0:
        print(f"  {i+len(batch_seqs)}/{len(seqs)}")

df = pd.DataFrame(rows)
df.to_csv(output_path, sep="\t", index=False)
print(f"Written {len(df)} rows × {len(df.columns)-1} features → {output_path}")
PYEOF

echo "Done chunk ${CHUNK}"
date
