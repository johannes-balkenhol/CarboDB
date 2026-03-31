#!/bin/bash
# run_interproscan_wbbi206.sh — no GNU parallel
set -uo pipefail

WD=/storage/users/job37yv/Projects/CarboDB_v3
IPR=${WD}/data/dbs/interpro/interproscan-5.72-103.0/interproscan.sh
CHUNKS=${WD}/data/interim/fasta_chunks
OUT=${WD}/data/features/interpro
TMP=${WD}/data/interim/ipr_tmp
LOG=${WD}/logs/ipr_wbbi206.log
TOTAL=2381
PARALLEL=4
CORES=8

mkdir -p ${OUT} ${TMP}
cd ${WD}

echo "=== InterProScan wbbi206 $(date) ===" | tee -a ${LOG}
DONE=$(ls ${OUT}/ipr_*.tsv 2>/dev/null | wc -l)
echo "Already done: ${DONE}/${TOTAL}" | tee -a ${LOG}

process_chunk() {
    local CHUNK=$(printf "%04d" $1)
    local OUTPUT=${OUT}/ipr_${CHUNK}.tsv
    local INPUT=${CHUNKS}/chunk_${CHUNK}.fasta
    local CLEAN=${TMP}/clean_${CHUNK}.fasta
    local ITMP=${TMP}/ipr_${CHUNK}_tmp
    [ -f "${OUTPUT}" ] && return 0
    [ ! -f "${INPUT}" ] && return 0
    python3 -c "
with open('${INPUT}') as fin, open('${CLEAN}','w') as fout:
    for line in fin:
        fout.write(line) if line.startswith('>') else fout.write(line.strip().replace('*','')+'\n')
"
    ${IPR} -i ${CLEAN} -o ${OUTPUT} -f tsv \
        -appl Pfam,ProSiteProfiles,ProSitePatterns,PANTHER,Gene3D,TIGRFAM,SUPERFAMILY,CDD,HAMAP \
        -dp --cpu ${CORES} -T ${ITMP} --disable-precalc > /dev/null 2>&1
    rm -f ${CLEAN}; rm -rf ${ITMP}
    CHUNK_VAR="${CHUNK}" OUTPUT_VAR="${OUTPUT}" INPUT_VAR="${INPUT}" \
    python3 -c "
import os, json, sys
from pathlib import Path
import pandas as pd
chunk=os.environ['CHUNK_VAR']; output=Path(os.environ['OUTPUT_VAR']); inp=os.environ['INPUT_VAR']
parsed=output.with_suffix('.parsed.tsv')
cols=['protein_id','md5','length','database','accession','description','start','stop','evalue','status','date','interpro_acc','interpro_desc','go_terms','pathways']
rows=[]
if output.exists() and output.stat().st_size>0:
    try:
        df=pd.read_csv(output,sep='\t',header=None,names=cols[:15],on_bad_lines='skip',dtype=str).fillna('')
        df['cdb_id']=df['protein_id'].str.split('|').str[0]
        for cdb_id,grp in df.groupby('cdb_id'):
            row={'cdb_id':cdb_id}
            for db,key,n in [('PANTHER','panther_family',1),('Gene3D','cath_superfamily',1),('TIGRFAM','tigrfam_hits_json',0),('SUPERFAMILY','superfamily_json',0),('CDD','cdd_hits_json',0),('HAMAP','hamap_hits_json',0),('ProSiteProfiles','prosite_profiles_json',0),('ProSitePatterns','prosite_patterns_json',0)]:
                sub=grp[grp['database']==db]
                if db=='PANTHER': row['panther_family']=sub['accession'].iloc[0].split(':')[0] if not sub.empty else ''; row['panther_subfamily']=sub['accession'].iloc[0].split(':')[1] if not sub.empty and ':' in sub['accession'].iloc[0] else ''; row['n_panther']=len(sub)
                elif db=='Gene3D': row['cath_superfamily']=sub['accession'].iloc[0] if not sub.empty else ''; row['gene3d_domains_json']=json.dumps(sub['accession'].tolist()); row['n_gene3d']=len(sub)
                elif db=='TIGRFAM': row['tigrfam_hits_json']=json.dumps(sub['accession'].tolist()); row['n_tigrfam']=len(sub)
                elif db=='ProSiteProfiles': row['prosite_profiles_json']=json.dumps(sub['accession'].tolist()); row['n_prosite_prof']=len(sub)
                elif db=='ProSitePatterns': row['prosite_patterns_json']=json.dumps(sub['accession'].tolist()); row['n_prosite_pat']=len(sub)
                else: row[key]=json.dumps(sub['accession'].tolist())
            row['raw_ipr_json']=grp[['database','accession','interpro_acc']].to_json(orient='records')
            rows.append(row)
    except Exception as e: print(f'Parse error: {e}',file=sys.stderr)
seen=set(r['cdb_id'] for r in rows)
with open(inp) as f:
    for line in f:
        if line.startswith('>'):
            cdb_id=line[1:].split('|')[0].strip()
            if cdb_id not in seen:
                rows.append({'cdb_id':cdb_id,'panther_family':'','panther_subfamily':'','cath_superfamily':'','gene3d_domains_json':'[]','tigrfam_hits_json':'[]','superfamily_json':'[]','cdd_hits_json':'[]','hamap_hits_json':'[]','prosite_profiles_json':'[]','prosite_patterns_json':'[]','n_panther':0,'n_gene3d':0,'n_tigrfam':0,'n_prosite_prof':0,'n_prosite_pat':0,'raw_ipr_json':'[]'})
pd.DataFrame(rows).to_csv(parsed,sep='\t',index=False)
output.unlink(missing_ok=True); parsed.rename(output)
print(f'chunk {chunk}: {len(rows)} rows')
"
}
export -f process_chunk
export OUT CHUNKS IPR TMP CORES

pids=()
for i in $(seq 1 ${TOTAL}); do
    process_chunk $i >> ${LOG} 2>&1 &
    pids+=($!)
    if [ ${#pids[@]} -ge ${PARALLEL} ]; then
        wait ${pids[0]}
        pids=("${pids[@]:1}")
    fi
    if [ $((i % 50)) -eq 0 ]; then
        DONE=$(ls ${OUT}/ipr_*.tsv 2>/dev/null | wc -l)
        echo "Progress: ${DONE}/${TOTAL} $(date)" | tee -a ${LOG}
    fi
done
wait
DONE=$(ls ${OUT}/ipr_*.tsv 2>/dev/null | wc -l)
echo "=== DONE: ${DONE}/${TOTAL} $(date) ===" | tee -a ${LOG}
