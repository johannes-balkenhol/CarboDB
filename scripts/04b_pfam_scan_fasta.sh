#!/usr/bin/env bash
# 04b_pfam_scan_fasta.sh -- hmmscan vs Pfam-A for any FASTA (e.g. GMGC candidates), KEEP domtblout.
# Same settings as v5 / 04b_pfam_rescan.sh: --noali -E 1e-3 --domE 1e-3. Local SSD, resumable per chunk.
# Usage: bash scripts/04b_pfam_scan_fasta.sh <FASTA> <NAME> <LOCAL_SSD_DIR> [NJOBS=15] [CPU_PER_JOB=4]
set -euo pipefail
FA=${1:?FASTA}; NAME=${2:?name}; SCR=${3:?local SSD dir}; NJOBS=${4:-15}; CPJ=${5:-4}; CHUNK=2000
DATA=/storage/users/job37yv/Projects/CarboDB_data
PFAM=$DATA/raw/dbs/pfam/Pfam-A.hmm
OUT=$DATA/primary/annotations/pfam_$NAME/domtblout
W=$SCR/pfam_$NAME; mkdir -p $W/chunks $W/work $OUT
command -v hmmscan >/dev/null || { echo "hmmscan not in PATH (conda activate carboxylase)"; exit 1; }
[ -s $W/Pfam-A.hmm ] || cp $PFAM* $W/
if [ ! -f $W/chunks/.done ]; then
  awk -v n=$CHUNK -v pre="$W/chunks/c_" '/^>/{if(c%n==0){if(f)close(f); f=sprintf("%s%05d.fasta",pre,c/n)} c++} {print > f}' $FA
  touch $W/chunks/.done
fi
echo "$(ls $W/chunks/c_*.fasta | wc -l) chunks, $(grep -c '>' $FA) sequences from $FA"
cat > $OUT/run.json <<J
{"stage":"04b_pfam_scan_fasta","name":"$NAME","input":"$FA","input_sha256":"$(sha256sum $FA | cut -d' ' -f1)",
 "tool":"hmmscan","tool_version":"$(hmmscan -h | sed -n 2p | sed 's/^# //')","params":"--noali -E 1e-3 --domE 1e-3 --cpu $CPJ",
 "pfam_sha256":"$(cut -d' ' -f1 $DATA/raw/dbs/pfam/Pfam-A.hmm.sha256)","host":"$(hostname)","started":"$(date -Is)","parallel_jobs":$NJOBS}
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
ls $W/chunks/c_*.fasta | xargs -P $NJOBS -I{} bash -c 'run_chunk {}'
echo "ALL DONE $(date -Is): $(ls $OUT/*.domtblout.gz | wc -l) chunks in $OUT"
