#!/bin/bash
#SBATCH --job-name=carbodb_hmmer
#SBATCH --partition=hades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=200G
#SBATCH --output=/storage/users/job37yv/Projects_shared/CarboDB_v3/logs/hmmer_%A_%a.log
#SBATCH --error=/storage/users/job37yv/Projects_shared/CarboDB_v3/logs/hmmer_%A_%a.log
#SBATCH --array=1-2380%50
# %50 = max 50 jobs running simultaneously to avoid overwhelming the node

# CarboxyDB — HMMER Pfam domain annotation
# Runs hmmscan on one FASTA chunk per array job
#
# Submit:
#   sbatch scripts/04b_hmmer.sh
#
# Check progress:
#   squeue -u $USER
#   ls data/features/domains/ | wc -l   # should reach 2380
#
# After all jobs done, merge:
#   python scripts/04b_hmmer_merge.py

set -euo pipefail

PROJECT=/storage/users/job37yv/Projects_shared/CarboDB_v3
PFAM_HMM=${PROJECT}/data/dbs/pfam/Pfam-A.hmm
CHUNKS=${PROJECT}/data/interim/fasta_chunks
OUT_DIR=${PROJECT}/data/features/domains
LOGS=${PROJECT}/logs

mkdir -p ${OUT_DIR}

# Pad array task ID to 4 digits
CHUNK=$(printf "%04d" ${SLURM_ARRAY_TASK_ID})
INPUT=${CHUNKS}/chunk_${CHUNK}.fasta
OUTPUT=${OUT_DIR}/hmmer_${CHUNK}.tsv
DOMTBL=${OUT_DIR}/hmmer_${CHUNK}.domtblout

# Skip if already done
if [ -f "${OUTPUT}" ]; then
    echo "Already done: ${OUTPUT}"
    exit 0
fi

if [ ! -f "${INPUT}" ]; then
    echo "Input not found: ${INPUT}"
    exit 1
fi

echo "=== HMMER chunk ${CHUNK} ==="
echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
date

# Run hmmscan
hmmscan \
    --domtblout ${DOMTBL} \
    --noali \
    --cpu ${SLURM_CPUS_PER_TASK} \
    -E 1e-3 \
    --domE 1e-3 \
    ${PFAM_HMM} \
    ${INPUT} \
    > /dev/null

echo "hmmscan done, parsing results..."

# Parse domtblout to TSV — pass variables via environment (heredoc with 'PYEOF' blocks expansion)
CHUNK_VAR="${CHUNK}" DOMTBL_VAR="${DOMTBL}" OUTPUT_VAR="${OUTPUT}" INPUT_VAR="${INPUT}" \
python3 << PYEOF
import os, sys, json
from pathlib import Path
import pandas as pd

chunk  = os.environ["CHUNK_VAR"]
domtbl = Path(os.environ["DOMTBL_VAR"])
output = Path(os.environ["OUTPUT_VAR"])
input_fasta = os.environ["INPUT_VAR"]

if not domtbl.exists() or domtbl.stat().st_size == 0:
    output.write_text("cdb_id\tpfam_hits_json\tpfam_n_hits\n")
    print(f"No hits for chunk {chunk}")
    sys.exit(0)

hits = {}
with open(domtbl) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        pfam_acc = parts[1].split(".")[0]
        query    = parts[3]
        cdb_id   = query.split("|")[0]
        evalue   = float(parts[6])
        if evalue <= 1e-3:
            if cdb_id not in hits:
                hits[cdb_id] = []
            if pfam_acc not in hits[cdb_id]:
                hits[cdb_id].append(pfam_acc)

CARBOXY_PFAM = [
    "PF00016","PF02788","PF00101","PF00194","PF03119",
    "PF00311","PF00821","PF02785","PF00364","PF01039",
    "PF02786","PF02787","PF00289","PF01309","PF03599",
    "PF03590","PF00384","PF00682",
]

rows = []
for cdb_id, pfam_list in hits.items():
    row = {"cdb_id": cdb_id, "pfam_hits_json": json.dumps(sorted(pfam_list)), "pfam_n_hits": len(pfam_list)}
    for pf in CARBOXY_PFAM:
        row[f"pfam_{pf}"] = 1 if pf in pfam_list else 0
    rows.append(row)

seen = set(hits.keys())
with open(input_fasta) as f:
    for line in f:
        if line.startswith(">"):
            cdb_id = line[1:].split("|")[0].strip()
            if cdb_id not in seen:
                row = {"cdb_id": cdb_id, "pfam_hits_json": "[]", "pfam_n_hits": 0}
                for pf in CARBOXY_PFAM:
                    row[f"pfam_{pf}"] = 0
                rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(output, sep="\t", index=False)
print(f"Written {len(df)} rows to {output}")
PYEOF

echo "Done chunk ${CHUNK}"
date
