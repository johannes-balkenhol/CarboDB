#!/bin/bash
#SBATCH --job-name=carbodb_ipr
#SBATCH --partition=hades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/ipr_%A_%a.log
#SBATCH --error=logs/ipr_%A_%a.log
#SBATCH --array=1-2380%20
# %20 = max 20 jobs at once (InterProScan is memory-heavy)

# CarboxyDB — InterProScan annotation
# Applications: Pfam, ProSiteProfiles, ProSitePatterns,
#               PANTHER, Gene3D, TIGRFAM, SUPERFAMILY, CDD, HAMAP
#
# Submit:
#   sbatch scripts/04c_interproscan.sh
#
# Monitor:
#   squeue -u $USER -n carbodb_ipr
#   ls data/features/interpro/*.tsv | wc -l

set -euo pipefail

PROJECT=/storage/users/job37yv/Projects/CarboDB_v3
IPR=${PROJECT}/data/dbs/interpro/interproscan-5.72-103.0/interproscan.sh
CHUNKS=${PROJECT}/data/interim/fasta_chunks
OUT_DIR=${PROJECT}/data/features/interpro
TMPDIR=${PROJECT}/data/interim/ipr_tmp

mkdir -p ${OUT_DIR} ${TMPDIR}

CHUNK_OFFSET=${CHUNK_OFFSET:-0}
CHUNK=$(printf "%04d" $((${SLURM_ARRAY_TASK_ID} + ${CHUNK_OFFSET})))
INPUT=${CHUNKS}/chunk_${CHUNK}.fasta
OUTPUT=${OUT_DIR}/ipr_${CHUNK}.tsv

if [ -f "${OUTPUT}" ]; then
    echo "Already done: ${OUTPUT}"
    exit 0
fi

if [ ! -f "${INPUT}" ]; then
    echo "Input not found: ${INPUT}"
    exit 1
fi

echo "=== InterProScan chunk ${CHUNK} ==="
date

# Strip stop codons (*) from sequences — InterProScan rejects them
CLEAN_INPUT=${TMPDIR}/clean_${CHUNK}.fasta
python3 -c "
with open('${INPUT}') as fin, open('${CLEAN_INPUT}', 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            fout.write(line)
        else:
            fout.write(line.strip().replace('*','') + '\n')
"

# Run InterProScan
${IPR} \
    -i  ${CLEAN_INPUT} \
    -o  ${OUTPUT} \
    -f  tsv \
    -appl Pfam,ProSiteProfiles,ProSitePatterns,PANTHER,Gene3D,TIGRFAM,SUPERFAMILY,CDD,HAMAP \
    -dp \
    --cpu ${SLURM_CPUS_PER_TASK} \
    -T  ${TMPDIR}/ipr_${CHUNK}_tmp \
    --disable-precalc

# Clean up temp files
rm -f ${CLEAN_INPUT}
rm -rf ${TMPDIR}/ipr_${CHUNK}_tmp

# Parse TSV to structured format
python3 << 'PYEOF'
import json
import sys
from pathlib import Path
import pandas as pd

chunk  = "${CHUNK}"
output = Path("${OUTPUT}")
parsed = output.with_suffix(".parsed.tsv")

if not output.exists() or output.stat().st_size == 0:
    parsed.write_text("cdb_id\tpanther_family\tpanther_subfamily\tcath_superfamily\tn_panther\tn_gene3d\tn_tigrfam\tn_prosite_prof\tn_prosite_pat\traw_ipr_json\n")
    print(f"Empty output for chunk {chunk}")
    sys.exit(0)

# InterProScan TSV columns:
# 0=protein_id 1=md5 2=length 3=database 4=accession 5=description
# 6=start 7=stop 8=evalue 9=status 10=date 11=interpro_acc 12=interpro_desc
# 13=GO 14=pathway
cols = ["protein_id","md5","length","database","accession","description",
        "start","stop","evalue","status","date","interpro_acc","interpro_desc",
        "go_terms","pathways"]

try:
    df = pd.read_csv(output, sep="\t", header=None, names=cols[:15],
                     on_bad_lines="skip", dtype=str).fillna("")
except Exception as e:
    print(f"Parse error: {e}")
    sys.exit(0)

# cdb_id is first field of protein_id (format: CDB000001|uniprot|ec|label...)
df["cdb_id"] = df["protein_id"].str.split("|").str[0]

# Aggregate per cdb_id
rows = []
for cdb_id, grp in df.groupby("cdb_id"):
    row = {"cdb_id": cdb_id}

    # PANTHER
    panther = grp[grp["database"] == "PANTHER"]
    if not panther.empty:
        acc = panther["accession"].iloc[0]
        parts = acc.split(":")
        row["panther_family"]    = parts[0] if parts else ""
        row["panther_subfamily"] = parts[1] if len(parts) > 1 else ""
        row["n_panther"] = len(panther)
    else:
        row["panther_family"] = row["panther_subfamily"] = ""
        row["n_panther"] = 0

    # Gene3D / CATH
    gene3d = grp[grp["database"] == "Gene3D"]
    row["gene3d_domains_json"] = json.dumps(gene3d["accession"].tolist())
    row["cath_superfamily"]    = gene3d["accession"].iloc[0] if not gene3d.empty else ""
    row["n_gene3d"] = len(gene3d)

    # TIGRFAM
    tigr = grp[grp["database"] == "TIGRFAM"]
    row["tigrfam_hits_json"] = json.dumps(tigr["accession"].tolist())
    row["n_tigrfam"] = len(tigr)

    # SUPERFAMILY
    sfam = grp[grp["database"] == "SUPERFAMILY"]
    row["superfamily_json"] = json.dumps(sfam["accession"].tolist())

    # CDD
    cdd = grp[grp["database"] == "CDD"]
    row["cdd_hits_json"] = json.dumps(cdd["accession"].tolist())

    # HAMAP
    hamap = grp[grp["database"] == "HAMAP"]
    row["hamap_hits_json"] = json.dumps(hamap["accession"].tolist())

    # ProSite Profiles
    psp = grp[grp["database"] == "ProSiteProfiles"]
    row["prosite_profiles_json"] = json.dumps(psp["accession"].tolist())
    row["n_prosite_prof"] = len(psp)

    # ProSite Patterns
    pspa = grp[grp["database"] == "ProSitePatterns"]
    row["prosite_patterns_json"] = json.dumps(pspa["accession"].tolist())
    row["n_prosite_pat"] = len(pspa)

    # Full raw JSON
    row["raw_ipr_json"] = grp[["database","accession","description","interpro_acc"]].to_json(orient="records")

    rows.append(row)

# Add zero rows for sequences with no hits
seen = set(r["cdb_id"] for r in rows)
with open("${INPUT}") as f:
    for line in f:
        if line.startswith(">"):
            cdb_id = line[1:].split("|")[0].strip()
            if cdb_id not in seen:
                rows.append({
                    "cdb_id": cdb_id,
                    "panther_family": "", "panther_subfamily": "",
                    "cath_superfamily": "",
                    "gene3d_domains_json": "[]", "tigrfam_hits_json": "[]",
                    "superfamily_json": "[]", "cdd_hits_json": "[]",
                    "hamap_hits_json": "[]",
                    "prosite_profiles_json": "[]", "prosite_patterns_json": "[]",
                    "n_panther": 0, "n_gene3d": 0, "n_tigrfam": 0,
                    "n_prosite_prof": 0, "n_prosite_pat": 0,
                    "raw_ipr_json": "[]",
                })

out_df = pd.DataFrame(rows)
out_df.to_csv(parsed, sep="\t", index=False)

# Replace raw IPR output with parsed version
output.unlink()
parsed.rename(output)
print(f"Parsed {len(out_df)} sequences for chunk {chunk}")
PYEOF

echo "Done chunk ${CHUNK}"
date
