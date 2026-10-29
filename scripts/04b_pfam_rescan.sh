#!/usr/bin/env bash
# 04b_pfam_rescan.sh -- CarboDB v5.2: rerun hmmscan vs Pfam-A on all sequences, KEEP domtblout.
# Same parameters as v5 (04b_hmmer.sh): hmmscan --noali -E 1e-3 --domE 1e-3
# Runs on local SSD, resumable per chunk, each finished chunk copied to /storage immediately.
# Usage: bash scripts/04b_pfam_rescan.sh <LOCAL_SSD_DIR> [NJOBS=16] [CPU_PER_JOB=4]
set -euo pipefail
SCR=${1:?give a local SSD directory}; NJOBS=${2:-16}; CPJ=${3:-4}; CHUNK=2000
DB=/storage/users/job37yv/Projects/CarboDB_v3/data/primary/carbodb.sqlite
PFAM_SRC=/storage/users/job37yv/Projects/CarboDB_v3/data/dbs/pfam/Pfam-A.hmm
DATA=/storage/users/job37yv/Projects/CarboDB_data
OUT=$DATA/primary/annotations/pfam/domtblout
RAWDB=$DATA/raw/dbs/pfam
W=$SCR/carbodb_pfam_rescan; mkdir -p $W/chunks $W/work $OUT $RAWDB
command -v hmmscan >/dev/null || { echo "hmmscan not in PATH (conda activate carboxylase)"; exit 1; }
echo "free on SSD: $(df -h $SCR | tail -1 | awk '{print $4}')"
if [ ! -s $RAWDB/Pfam-A.hmm ]; then
  cp -L ${PFAM_SRC}* $RAWDB/
  (cd $RAWDB && sha256sum Pfam-A.hmm > Pfam-A.hmm.sha256)
  ls -la "$(dirname "$(readlink -f $PFAM_SRC)")" > $RAWDB/SOURCE_listing.txt
  echo "copied from $(readlink -f $PFAM_SRC) on $(date +%F)" > $RAWDB/SOURCE.md
fi
[ -s $W/Pfam-A.hmm ] || cp $RAWDB/Pfam-A.hmm* $W/
NMODELS=$(grep -c '^NAME' $W/Pfam-A.hmm)
if [ ! -f $W/chunks/.done ]; then
  for pass in 1 2; do
    if [ $pass = 1 ]; then COND="s.label IN (1,2) OR p.is_co2_pred=1 OR d.sequence_id IS NULL"
    else COND="NOT (s.label IN (1,2) OR p.is_co2_pred=1 OR d.sequence_id IS NULL)"; fi
    sqlite3 -readonly $DB "SELECT '>'||s.cdb_id||'|'||s.uniprot_id||'|'||s.ec_number||'|label='||s.label, s.sequence
      FROM sequences s
      LEFT JOIN (SELECT sequence_id, MAX(is_co2_pred) AS is_co2_pred FROM predictions GROUP BY sequence_id) p ON p.sequence_id=s.id
      LEFT JOIN features_domains d ON d.sequence_id=s.id
      WHERE $COND ORDER BY s.id;" -separator $'\n' \
    | awk -v n=$CHUNK -v pre="$W/chunks/p${pass}_" '/^>/{if(c%n==0){if(f)close(f); f=sprintf("%s%05d.fasta",pre,c/n)} c++} {print > f}'
    echo "pass $pass: $(ls $W/chunks/p${pass}_*.fasta | wc -l) chunks"
  done
  touch $W/chunks/.done
fi
cat > $OUT/run.json <<J
{"stage":"04b_pfam_rescan","tool":"hmmscan","tool_version":"$(hmmscan -h | sed -n 2p | sed 's/^# //')",
 "params":"--noali -E 1e-3 --domE 1e-3 --cpu $CPJ","pfam_sha256":"$(cut -d' ' -f1 $RAWDB/Pfam-A.hmm.sha256)",
 "pfam_models":$NMODELS,"host":"$(hostname)","started":"$(date -Is)","chunk_size":$CHUNK,"parallel_jobs":$NJOBS}
J
export W OUT CPJ
run_chunk() {
  f=$1; b=$(basename $f .fasta)
  [ -s $OUT/$b.domtblout.gz ] && return 0
  hmmscan --noali --cpu $CPJ -E 1e-3 --domE 1e-3 --domtblout $W/work/$b.domtblout $W/Pfam-A.hmm $f > /dev/null \
    && gzip -f $W/work/$b.domtblout && cp $W/work/$b.domtblout.gz $OUT/$b.domtblout.gz.part \
    && mv $OUT/$b.domtblout.gz.part $OUT/$b.domtblout.gz && rm $W/work/$b.domtblout.gz \
    && echo "$(date +%T) done $b ($(ls $OUT/*.domtblout.gz | wc -l) total)"
}
export -f run_chunk
ls $W/chunks/p1_*.fasta $W/chunks/p2_*.fasta | xargs -P $NJOBS -I{} bash -c 'run_chunk {}'
echo "ALL DONE $(date -Is): $(ls $OUT/*.domtblout.gz | wc -l) chunks in $OUT"
