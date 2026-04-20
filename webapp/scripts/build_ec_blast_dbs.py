#!/usr/bin/env python3
"""
One-time setup: build per-EC BLAST databases with CD-HIT 70% reduction.

Output:
  data/blast_ec_dbs/raw/ec_{EC}.fa        — raw extracted sequences
  data/blast_ec_dbs/clustered/ec_{EC}.fa  — CD-HIT clustered reps + .phr/.pin/.psq indices
  data/blast_ec_dbs/manifest.json         — per-EC stats

FASTA header format: >uniprot_id|EC|source|reviewed|organism
  lets us recover everything without extra DB hits when parsing BLAST output.

Runtime: ~20-40 min for all 27 EC classes, dominated by 4.1.1.39 (154k seqs).
"""

import sqlite3
import subprocess
import json
import time
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
os.chdir(ROOT)

DB_PATH = Path("data/primary/carbodb.sqlite")
OUT_DIR = Path("data/blast_ec_dbs")
RAW_DIR = OUT_DIR / "raw"
CLUST_DIR = OUT_DIR / "clustered"
IDENTITY = 0.70
THREADS = 4
MEM_MB = 4000


def get_ec_classes(conn):
    cur = conn.execute("""
        SELECT ec_number, COUNT(*) as n
        FROM sequences WHERE label=1 AND seq_valid=1
        GROUP BY ec_number
        ORDER BY n DESC
    """)
    return [(row[0], row[1]) for row in cur]


def dump_fasta(conn, ec, fa_path):
    cur = conn.execute("""
        SELECT uniprot_id, organism, source, reviewed, sequence
        FROM sequences
        WHERE label=1 AND seq_valid=1 AND ec_number=?
    """, (ec,))
    n = 0
    with open(fa_path, "w") as f:
        for uid, org, src, rev, seq in cur:
            if not seq or not uid:
                continue
            org_clean = (org or "unknown").replace("|", "_").replace(" ", "_")[:80]
            header = f">{uid}|{ec}|{src}|{rev}|{org_clean}"
            f.write(f"{header}\n{seq}\n")
            n += 1
    return n


def cluster_fasta(input_fa, output_fa, pident=IDENTITY):
    cmd = [
        "cd-hit",
        "-i", str(input_fa),
        "-o", str(output_fa),
        "-c", str(pident),
        "-n", "5",
        "-M", str(MEM_MB),
        "-T", str(THREADS),
        "-d", "0",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    CD-HIT FAILED: {result.stderr[:500]}", file=sys.stderr)
        return False
    return True


def makeblastdb_fa(fa_path, title):
    cmd = [
        "makeblastdb",
        "-in", str(fa_path),
        "-dbtype", "prot",
        "-title", title,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    makeblastdb FAILED: {result.stderr[:500]}", file=sys.stderr)
        return False
    return True


def count_seqs(fa_path):
    n = 0
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def main():
    if not DB_PATH.exists():
        print(f"DB not found: {DB_PATH}", file=sys.stderr)
        sys.exit(1)

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    CLUST_DIR.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(DB_PATH)
    ec_list = get_ec_classes(conn)

    print(f"Building BLAST DBs for {len(ec_list)} EC classes")
    print(f"Output: {OUT_DIR}/")
    print(f"CD-HIT identity: {IDENTITY}, threads: {THREADS}, mem: {MEM_MB} MB\n")

    manifest = {}
    total_start = time.time()

    for i, (ec, n_total) in enumerate(ec_list, 1):
        safe_ec = ec.replace(".", "_")
        raw_fa = RAW_DIR / f"ec_{safe_ec}.fa"
        clust_fa = CLUST_DIR / f"ec_{safe_ec}.fa"

        print(f"[{i:>2}/{len(ec_list)}] EC {ec} ({n_total} seqs) ... ", end="", flush=True)
        t0 = time.time()

        n_dumped = dump_fasta(conn, ec, raw_fa)
        if n_dumped == 0:
            print("skip (empty after dump)")
            continue

        # For tiny classes, skip CD-HIT (not worth overhead, nothing to reduce)
        if n_dumped < 50:
            import shutil
            shutil.copy(raw_fa, clust_fa)
            n_clust = n_dumped
        else:
            ok = cluster_fasta(raw_fa, clust_fa)
            if not ok:
                print("FAILED")
                continue
            n_clust = count_seqs(clust_fa)

        ok = makeblastdb_fa(clust_fa, f"CarboDB_{safe_ec}")
        if not ok:
            print("makeblastdb failed")
            continue

        dt = time.time() - t0
        reduction = round(100 * (1 - n_clust / max(n_dumped, 1)), 1)
        manifest[ec] = {
            "n_original": n_dumped,
            "n_clustered": n_clust,
            "reduction_pct": reduction,
            "db_path": str(clust_fa),
            "build_seconds": round(dt, 1),
        }
        print(f"-> {n_clust} ({reduction}% reduction, {dt:.1f}s)")

    # Save manifest
    manifest_path = OUT_DIR / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    conn.close()

    total_dt = time.time() - total_start
    total_clustered = sum(v["n_clustered"] for v in manifest.values())
    total_original = sum(v["n_original"] for v in manifest.values())
    print(f"\n=== DONE ===")
    print(f"Total: {total_original} -> {total_clustered} "
          f"({100*(1-total_clustered/max(total_original,1)):.1f}% reduction)")
    print(f"Walltime: {total_dt/60:.1f} min")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
