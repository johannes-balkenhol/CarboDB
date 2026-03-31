#!/bin/bash
# run_hmmer_wbbi206.sh — no GNU parallel, uses bash background jobs
set -uo pipefail

WD=/storage/users/job37yv/Projects/CarboDB_v3
PFAM=${WD}/data/dbs/pfam/Pfam-A.hmm
CHUNKS=${WD}/data/interim/fasta_chunks
OUT=${WD}/data/features/domains
LOG=${WD}/logs/hmmer_wbbi206.log
TOTAL=2381
PARALLEL=20
CORES=4

mkdir -p ${OUT}
cd ${WD}

echo "=== HMMER wbbi206 $(date) ===" | tee -a ${LOG}
DONE=$(ls ${OUT}/hmmer_*.tsv 2>/dev/null | wc -l)
echo "Already done: ${DONE}/${TOTAL}" | tee -a ${LOG}

process_chunk() {
    local CHUNK=$(printf "%04d" $1)
    local OUTPUT=${OUT}/hmmer_${CHUNK}.tsv
    local INPUT=${CHUNKS}/chunk_${CHUNK}.fasta
    local DOMTBL=${OUT}/hmmer_${CHUNK}.domtblout
    [ -f "${OUTPUT}" ] && return 0
    [ ! -f "${INPUT}" ] && return 0
    hmmscan --domtblout ${DOMTBL} --noali --cpu ${CORES} -E 1e-3 \
            ${PFAM} ${INPUT} > /dev/null 2>&1
    CHUNK_VAR="${CHUNK}" DOMTBL_VAR="${DOMTBL}" \
    OUTPUT_VAR="${OUTPUT}" INPUT_VAR="${INPUT}" \
    python3 -c "
import os, json
from pathlib import Path
import pandas as pd
chunk       = os.environ['CHUNK_VAR']
domtbl      = Path(os.environ['DOMTBL_VAR'])
output      = Path(os.environ['OUTPUT_VAR'])
input_fasta = os.environ['INPUT_VAR']
CARBOXY_PFAM = ['PF00016','PF02788','PF00101','PF00194','PF03119',
    'PF00311','PF00821','PF02785','PF00364','PF01039',
    'PF02786','PF02787','PF00289','PF01309','PF03599',
    'PF03590','PF00384','PF00682']
hits = {}
if domtbl.exists() and domtbl.stat().st_size > 0:
    with open(domtbl) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.split()
            if len(parts) < 7: continue
            pfam_acc = parts[1].split('.')[0]
            cdb_id   = parts[3].split('|')[0]
            if float(parts[6]) <= 1e-3:
                if cdb_id not in hits: hits[cdb_id] = []
                if pfam_acc not in hits[cdb_id]: hits[cdb_id].append(pfam_acc)
rows = []
for cdb_id, pfam_list in hits.items():
    row = {'cdb_id': cdb_id, 'pfam_hits_json': json.dumps(sorted(pfam_list)), 'pfam_n_hits': len(pfam_list)}
    for pf in CARBOXY_PFAM: row[f'pfam_{pf}'] = 1 if pf in pfam_list else 0
    rows.append(row)
seen = set(hits.keys())
with open(input_fasta) as f:
    for line in f:
        if line.startswith('>'):
            cdb_id = line[1:].split('|')[0].strip()
            if cdb_id not in seen:
                row = {'cdb_id': cdb_id, 'pfam_hits_json': '[]', 'pfam_n_hits': 0}
                for pf in CARBOXY_PFAM: row[f'pfam_{pf}'] = 0
                rows.append(row)
pd.DataFrame(rows).to_csv(output, sep='\t', index=False)
print(f'chunk {chunk}: {len(rows)} rows')
"
    rm -f ${DOMTBL}
}
export -f process_chunk
export OUT CHUNKS PFAM CORES

pids=()
for i in $(seq 1 ${TOTAL}); do
    process_chunk $i >> ${LOG} 2>&1 &
    pids+=($!)
    if [ ${#pids[@]} -ge ${PARALLEL} ]; then
        wait ${pids[0]}
        pids=("${pids[@]:1}")
    fi
    if [ $((i % 100)) -eq 0 ]; then
        DONE=$(ls ${OUT}/hmmer_*.tsv 2>/dev/null | wc -l)
        echo "Progress: ${DONE}/${TOTAL} (launched $i) $(date)" | tee -a ${LOG}
    fi
done
wait
DONE=$(ls ${OUT}/hmmer_*.tsv 2>/dev/null | wc -l)
echo "=== DONE: ${DONE}/${TOTAL} $(date) ===" | tee -a ${LOG}
