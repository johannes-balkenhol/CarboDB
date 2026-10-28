#!/usr/bin/env python3
"""
Random 1M-sequence subsample from the translated marine FASTA for the pilot scan.
Reservoir sampling — single streaming pass, modest peak RAM.
"""
import gzip
import random
import sys
from pathlib import Path

ROOT      = Path.home() / "Projects_shared/CarboDB_v3"
INPUT     = ROOT / "data/gmgc/interim/marine.95nr.faa.gz"
OUTPUT    = ROOT / "data/gmgc/interim/pilot_1M.faa"
N_SAMPLES = 1_000_000
SEED      = 42

random.seed(SEED)


def stream_fasta(path):
    """Yield (header, seq) tuples from a gzipped FASTA."""
    with gzip.open(path, "rt") as f:
        header = None
        buf = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header = line
                buf = []
            else:
                buf.append(line)
        if header is not None:
            yield header, "".join(buf)


reservoir = []
n_seen = 0
for entry in stream_fasta(INPUT):
    n_seen += 1
    if len(reservoir) < N_SAMPLES:
        reservoir.append(entry)
    else:
        j = random.randint(0, n_seen - 1)
        if j < N_SAMPLES:
            reservoir[j] = entry
    if n_seen % 1_000_000 == 0:
        print(f"  Seen {n_seen:,} / reservoir {len(reservoir):,}", file=sys.stderr)

print(f"Total scanned: {n_seen:,}", file=sys.stderr)
print(f"Sampled:       {len(reservoir):,}", file=sys.stderr)
print(f"Writing to {OUTPUT}", file=sys.stderr)

with open(OUTPUT, "w") as f:
    for header, seq in reservoir:
        f.write(header + "\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

print("Done.", file=sys.stderr)
