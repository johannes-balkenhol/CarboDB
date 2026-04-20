#!/usr/bin/env python3
"""
Build per-EC BLAST databases containing ONLY sequences with experimental Km.

Output:
  data/blast_ec_dbs_exp/ec_{EC}.fa       — raw extracted sequences + BLAST indices
  data/blast_ec_dbs_exp/manifest.json    — per-EC stats

No CD-HIT reduction needed — experimental pools are small (max ~1000).

Runtime: ~30 seconds for all ECs combined.
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
OUT_DIR = Path("data/blast_ec_dbs_exp")


def get_ec_classes_with_exp_km(conn):
    cur = conn.execute("""
        SELECT s.ec_number, COUNT(DISTINCT s.uniprot_id) as n
        FROM sequences s
        JOIN km_evidence k ON k.uniprot_id = s.uniprot_id
        WHERE s.label=1 AND s.seq_valid=1 AND k.evidence_tier=1
        GROUP BY s.ec_number
        ORDER BY n DESC
    """)
    return [(row[0], row[1]) for row in cur]


def dump_fasta(conn, ec, fa_path):
    cur = conn.execute("""
        SELECT DISTINCT s.uniprot_id, s.organism, s.source, s.reviewed, s.sequence
        FROM sequences s
        JOIN km_evidence k ON k.uniprot_id = s.uniprot_id
        WHERE s.label=1 AND s.seq_valid=1 AND s.ec_number=? AND k.evidence_tier=1
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


def makeblastdb_fa(fa_path, title):
    result = subprocess.run(
        ["makeblastdb", "-in", str(fa_path), "-dbtype", "prot", "-title", title],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"    makeblastdb FAILED: {result.stderr[:300]}", file=sys.stderr)
        return False
    return True


def main():
    if not DB_PATH.exists():
        print(f"DB not found: {DB_PATH}", file=sys.stderr)
        sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(DB_PATH)
    ec_list = get_ec_classes_with_exp_km(conn)

    print(f"Building experimental-Km BLAST DBs for {len(ec_list)} EC classes\n")

    manifest = {}
    total_start = time.time()

    for i, (ec, n_total) in enumerate(ec_list, 1):
        safe_ec = ec.replace(".", "_")
        fa = OUT_DIR / f"ec_{safe_ec}.fa"
        print(f"[{i:>2}/{len(ec_list)}] EC {ec} ({n_total} seqs w/ exp Km) ... ",
              end="", flush=True)
        t0 = time.time()
        n_dumped = dump_fasta(conn, ec, fa)
        if n_dumped == 0:
            print("skip (empty)")
            continue
        if not makeblastdb_fa(fa, f"CarboDB_exp_{safe_ec}"):
            print("makeblastdb failed")
            continue
        dt = time.time() - t0
        manifest[ec] = {
            "n_sequences": n_dumped,
            "db_path": str(fa),
            "build_seconds": round(dt, 2),
        }
        print(f"-> {n_dumped} seqs, {dt:.2f}s")

    with open(OUT_DIR / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)

    conn.close()
    total_dt = time.time() - total_start
    total = sum(v["n_sequences"] for v in manifest.values())
    print(f"\n=== DONE ===")
    print(f"Total: {total} experimental-Km sequences across {len(manifest)} ECs")
    print(f"Walltime: {total_dt:.1f}s")


if __name__ == "__main__":
    main()
