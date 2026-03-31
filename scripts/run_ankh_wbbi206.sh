#!/bin/bash
# run_ankh_wbbi206.sh — no GNU parallel, uses bash background jobs
set -uo pipefail

WD=/storage/users/job37yv/Projects/CarboDB_v3
CHUNKS=${WD}/data/interim/fasta_chunks
OUT=${WD}/data/features/ankh
LOG=${WD}/logs/ankh_wbbi206.log
TOTAL=2381
PARALLEL=8

mkdir -p ${OUT}
cd ${WD}

echo "=== Ankh wbbi206 $(date) ===" | tee -a ${LOG}
DONE=$(ls ${OUT}/ankh_*.tsv 2>/dev/null | wc -l)
echo "Already done: ${DONE}/${TOTAL}" | tee -a ${LOG}

process_chunk() {
    local CHUNK=$(printf "%04d" $1)
    local OUTPUT=${OUT}/ankh_${CHUNK}.tsv
    local INPUT=${CHUNKS}/chunk_${CHUNK}.fasta
    [ -f "${OUTPUT}" ] && return 0
    [ ! -f "${INPUT}" ] && return 0
    OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 TOKENIZERS_PARALLELISM=false \
    CHUNK_VAR="${CHUNK}" OUTPUT_VAR="${OUTPUT}" INPUT_VAR="${INPUT}" \
    /home/job37yv/miniforge3/envs/carboxylase/bin/python3 -c "
import os, torch
os.environ['OMP_NUM_THREADS']='4'; os.environ['MKL_NUM_THREADS']='4'
torch.set_num_threads(4)
import ankh, pandas as pd
from pathlib import Path
chunk  = os.environ['CHUNK_VAR']
output = Path(os.environ['OUTPUT_VAR'])
inp    = Path(os.environ['INPUT_VAR'])
model, tokenizer = ankh.load_large_model()
model.eval()
seqs, ids = [], []
with open(inp) as f:
    cid = sl = None
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if cid and sl: seqs.append(''.join(sl)); ids.append(cid)
            cid = line[1:].split('|')[0]; sl = []
        elif line: sl.append(line)
    if cid and sl: seqs.append(''.join(sl)); ids.append(cid)
rows = []
for i in range(0, len(seqs), 4):
    bs = [s[:1022] for s in seqs[i:i+4]]; bi = ids[i:i+4]
    enc = tokenizer.batch_encode_plus([list(s) for s in bs],
        add_special_tokens=True, padding=True,
        is_split_into_words=True, return_tensors='pt')
    with torch.no_grad():
        out = model(input_ids=enc['input_ids'], attention_mask=enc['attention_mask'])
    mask = enc['attention_mask'].unsqueeze(-1).float()
    embs = ((out.last_hidden_state * mask).sum(1) / mask.sum(1)).cpu().numpy()
    for cid2, emb in zip(bi, embs):
        row = {'cdb_id': cid2}
        for j,v in enumerate(emb): row[f'ankh_{j}'] = float(v)
        rows.append(row)
pd.DataFrame(rows).to_csv(output, sep='\t', index=False)
print(f'chunk {chunk}: {len(rows)} rows')
"
}
export -f process_chunk
export OUT CHUNKS

pids=()
for i in $(seq 1 ${TOTAL}); do
    process_chunk $i >> ${LOG} 2>&1 &
    pids+=($!)
    if [ ${#pids[@]} -ge ${PARALLEL} ]; then
        wait ${pids[0]}
        pids=("${pids[@]:1}")
    fi
    if [ $((i % 50)) -eq 0 ]; then
        DONE=$(ls ${OUT}/ankh_*.tsv 2>/dev/null | wc -l)
        echo "Progress: ${DONE}/${TOTAL} $(date)" | tee -a ${LOG}
    fi
done
wait
DONE=$(ls ${OUT}/ankh_*.tsv 2>/dev/null | wc -l)
echo "=== DONE: ${DONE}/${TOTAL} $(date) ===" | tee -a ${LOG}
