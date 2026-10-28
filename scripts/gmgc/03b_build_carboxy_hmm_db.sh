#!/bin/bash
# Build a restricted HMM library containing only the 18 carboxylase Pfams.
# Uses HMM NAMES (not accessions) because that's what hmmfetch SSI indexes.
set -euo pipefail

cd ~/Projects_shared/CarboDB_v3

PFAM_FULL=data/dbs/pfam/Pfam-A.hmm
PFAM_OUT=data/dbs/pfam/carboxy_pfams.hmm
PFAM_NAMES=data/dbs/pfam/carboxy_pfam_names.txt

if [ ! -f "$PFAM_NAMES" ]; then
    echo "ERROR: $PFAM_NAMES missing. Run the acc->name mapping block first."
    exit 1
fi

echo "Names to extract:"
cat "$PFAM_NAMES"
echo

# hmmfetch already built .h3m.ssi index, but it indexes by NAME, so we're good
# Extract the 18 HMMs to a new file (use -f for "fetch list from file")
echo "Extracting 18 HMMs..."
hmmfetch -f "$PFAM_FULL" "$PFAM_NAMES" > "$PFAM_OUT"

# Press the new HMM file for hmmsearch
echo "Pressing carboxy HMM database..."
rm -f "${PFAM_OUT}.h3"?
hmmpress "$PFAM_OUT"

echo
echo "=== Result ==="
ls -lh "$PFAM_OUT"*
echo
echo "HMMs in restricted library:"
grep "^NAME " "$PFAM_OUT"
