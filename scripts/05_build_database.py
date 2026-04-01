#!/usr/bin/env python3
"""
05_build_database.py
====================
CarboxyDB — Step 05: Build carbodb.sqlite from master.tsv + feature TSVs.

Creates the SQLite database using schema v2.2. Injects all available
feature tables. Partial feature tables (InterPro, Ankh) are injected
as far as available — run this script again later to fill remaining rows.

Usage:
    # Full build (all available features)
    python scripts/05_build_database.py

    # Dry run — check counts without writing
    python scripts/05_build_database.py --dry-run

    # Update only specific tables (for incremental updates)
    python scripts/05_build_database.py --update interpro
    python scripts/05_build_database.py --update ankh

Output:
    data/primary/carbodb.sqlite   (~20-50GB depending on embeddings)

Runtime:
    ~2-4h for full build on 2.38M sequences
    ~30min for core tables only (no embeddings)
"""

import argparse
import json
import logging
import sqlite3
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("05_build_database")

DB_PATH    = PATHS.PRIMARY / "carbodb.sqlite"
SCHEMA_SQL = Path(__file__).resolve().parent / "schema.sql"

# ── Chunk sizes for memory-efficient loading ───────────────────────────────
CHUNK_SEQUENCES   = 50_000
CHUNK_FEATURES    = 10_000


# ═══════════════════════════════════════════════════════════════════════════
# Database helpers
# ═══════════════════════════════════════════════════════════════════════════

def get_connection(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA cache_size = -2000000")   # 2GB cache
    conn.execute("PRAGMA temp_store = MEMORY")
    return conn


def apply_schema(conn: sqlite3.Connection):
    log.info("Applying schema from %s", SCHEMA_SQL)
    with open(SCHEMA_SQL) as f:
        sql = f.read()
    conn.executescript(sql)
    conn.commit()
    log.info("Schema applied")


def table_count(conn: sqlite3.Connection, table: str) -> int:
    try:
        return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    except Exception:
        return 0


# ═══════════════════════════════════════════════════════════════════════════
# Step 1: sequences + id_map
# ═══════════════════════════════════════════════════════════════════════════

def build_sequences(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building sequences + id_map ──")

    existing = table_count(conn, "sequences")
    if existing > 0:
        log.info("  sequences already has %d rows — skipping", existing)
        return

    master = PATHS.MASTER_TSV
    log.info("  Reading %s ...", master)

    total = 0
    for chunk in tqdm(
        pd.read_csv(master, sep="\t", dtype=str, chunksize=CHUNK_SEQUENCES),
        desc="sequences"
    ):
        chunk = chunk.where(chunk.notna(), None)

        if dry_run:
            total += len(chunk)
            continue

        # sequences table
        seq_rows = []
        idmap_rows = []
        ec_rows = []
        km_rows = []

        for _, row in chunk.iterrows():
            seq_rows.append({
                "cdb_id":      row.get("cdb_id"),
                "uniprot_id":  row.get("uniprot_id"),
                "ec_number":   row.get("ec_number"),
                "label":       int(row.get("label", 0)),
                "source":      row.get("source"),
                "sequence":    row.get("sequence"),
                "length":      int(row.get("length", 0)) if row.get("length") else None,
                "organism":    row.get("organism"),
                "reviewed":    int(row.get("reviewed", 0)) if row.get("reviewed") else 0,
                "km_best_mM":  float(row["km_best_mM"]) if row.get("km_best_mM") else None,
                "km_log10_mM": float(row["km_log10_mM"]) if row.get("km_log10_mM") else None,
                "seq_valid":   1,
            })

            idmap_rows.append({
                "cdb_id":     row.get("cdb_id"),
                "uniprot_id": row.get("uniprot_id"),
            })

            # ec_evidence tier 1/2 from source
            ec_tier = 1 if row.get("source") in ("brenda",) else 2
            ec_rows.append({
                "uniprot_id":   row.get("uniprot_id"),
                "ec_number":    row.get("ec_number"),
                "source":       row.get("source", ""),
                "evidence_tier": ec_tier,
                "confidence":   None,
                "model_version": None,
            })

            # km_evidence tier 1 for experimental Km
            if row.get("km_best_mM"):
                km_rows.append({
                    "uniprot_id":   row.get("uniprot_id"),
                    "ec_number":    row.get("ec_number"),
                    "km_value_mM":  float(row["km_best_mM"]),
                    "km_log10_mM":  float(row["km_log10_mM"]) if row.get("km_log10_mM") else None,
                    "km_unit":      "mM",
                    "substrate":    "CO2",
                    "source":       row.get("source", "brenda"),
                    "evidence_tier": 1,
                    "commentary":   None,
                    "model_version": None,
                })

        conn.executemany("""
            INSERT OR IGNORE INTO sequences
            (cdb_id, uniprot_id, ec_number, label, source, sequence, length,
             organism, reviewed, km_best_mM, km_log10_mM, seq_valid)
            VALUES
            (:cdb_id, :uniprot_id, :ec_number, :label, :source, :sequence,
             :length, :organism, :reviewed, :km_best_mM, :km_log10_mM, :seq_valid)
        """, seq_rows)

        conn.executemany("""
            INSERT OR IGNORE INTO id_map (cdb_id, uniprot_id)
            VALUES (:cdb_id, :uniprot_id)
        """, idmap_rows)

        # Insert ec_evidence with sequence_id lookup
        conn.executemany("""
            INSERT OR IGNORE INTO ec_evidence
            (sequence_id, uniprot_id, ec_number, source, evidence_tier, confidence, model_version)
            SELECT s.id, :uniprot_id, :ec_number, :source, :evidence_tier, :confidence, :model_version
            FROM sequences s WHERE s.uniprot_id = :uniprot_id
        """, ec_rows)

        if km_rows:
            conn.executemany("""
                INSERT OR IGNORE INTO km_evidence
                (sequence_id, uniprot_id, ec_number, km_value_mM, km_log10_mM,
                 km_unit, substrate, source, evidence_tier, commentary, model_version)
                SELECT s.id, :uniprot_id, :ec_number, :km_value_mM, :km_log10_mM,
                       :km_unit, :substrate, :source, :evidence_tier, :commentary, :model_version
                FROM sequences s WHERE s.uniprot_id = :uniprot_id
            """, km_rows)

        conn.commit()
        total += len(chunk)

    log.info("  sequences: %d rows inserted", total)
    log.info("  id_map:    %d rows", table_count(conn, "id_map"))
    log.info("  ec_evidence: %d rows", table_count(conn, "ec_evidence"))
    log.info("  km_evidence: %d rows", table_count(conn, "km_evidence"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 2: Composition features
# ═══════════════════════════════════════════════════════════════════════════

def build_composition(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building features_composition ──")

    existing = table_count(conn, "features_composition")
    if existing > 0:
        log.info("  features_composition already has %d rows — skipping", existing)
        return

    feat_dir = PATHS.FEAT_COMP
    tsv_files = sorted(feat_dir.glob("composition_[0-9]*.tsv"))
    log.info("  Found %d composition TSV files", len(tsv_files))

    total = 0
    for tsv in tqdm(tsv_files, desc="composition"):
        df = pd.read_csv(tsv, sep="\t", dtype=str)
        df = df.where(df.notna(), None)

        if dry_run:
            total += len(df)
            continue

        # Get all column names from first file to build INSERT
        cols = [c for c in df.columns if c != "cdb_id"]

        for chunk_df in [df[i:i+CHUNK_FEATURES] for i in range(0, len(df), CHUNK_FEATURES)]:
            rows = []
            for _, row in chunk_df.iterrows():
                r = {"cdb_id": row["cdb_id"]}
                for c in cols:
                    try:
                        r[c] = float(row[c]) if row[c] is not None else None
                    except (ValueError, TypeError):
                        r[c] = row[c]  # keep as text for JSON cols
                rows.append(r)

            placeholders = ", ".join([f":{c}" for c in ["cdb_id"] + cols])
            col_names = ", ".join(["uniprot_id", "sequence_id"] + cols)

            conn.executemany(f"""
                INSERT OR IGNORE INTO features_composition
                (sequence_id, uniprot_id, {", ".join(cols)})
                SELECT s.id, s.uniprot_id, {", ".join([f':{{c}}' for c in cols]).replace('{c}', '').replace(':', '')}
                FROM sequences s WHERE s.cdb_id = :cdb_id
            """.replace(
                f", {', '.join([f':{{c}}' for c in cols]).replace('{c}', '').replace(':', '')}",
                ", " + ", ".join([f":{c}" for c in cols])
            ), rows)

            conn.commit()
            total += len(chunk_df)

    log.info("  features_composition: %d rows", table_count(conn, "features_composition"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 2b: Composition — cleaner version using direct INSERT
# ═══════════════════════════════════════════════════════════════════════════

def build_composition_v2(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building features_composition ──")

    existing = table_count(conn, "features_composition")
    if existing > 0:
        log.info("  features_composition already has %d rows — skipping", existing)
        return

    feat_dir = PATHS.FEAT_COMP
    tsv_files = sorted(feat_dir.glob("composition_[0-9]*.tsv"))
    log.info("  Found %d composition TSV files", len(tsv_files))

    if not tsv_files:
        log.warning("  No composition TSV files found")
        return

    # Build uid→id lookup
    log.info("  Building cdb_id→sequence_id lookup...")
    uid_map = dict(conn.execute("SELECT cdb_id, id FROM sequences").fetchall())
    log.info("  Lookup: %d entries", len(uid_map))

    total = 0
    for tsv in tqdm(tsv_files, desc="composition"):
        df = pd.read_csv(tsv, sep="\t")
        if dry_run:
            total += len(df)
            continue

        df["sequence_id"] = df["cdb_id"].map(uid_map)
        df["uniprot_id"]  = df["cdb_id"].map(
            dict(conn.execute("SELECT cdb_id, uniprot_id FROM sequences").fetchall())
        )
        df = df.dropna(subset=["sequence_id"])

        cols = [c for c in df.columns if c not in ("cdb_id",)]
        col_str = ", ".join(cols)
        ph_str  = ", ".join([f":{c}" for c in cols])

        conn.executemany(
            f"INSERT OR IGNORE INTO features_composition ({col_str}) VALUES ({ph_str})",
            df[cols].where(df[cols].notna(), None).to_dict("records")
        )
        conn.commit()
        total += len(df)

    log.info("  features_composition: %d rows inserted", table_count(conn, "features_composition"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 3: HMMER domain features
# ═══════════════════════════════════════════════════════════════════════════

def build_domains(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building features_domains ──")

    existing = table_count(conn, "features_domains")
    if existing > 0:
        log.info("  features_domains already has %d rows — skipping", existing)
        return

    feat_dir = PATHS.FEAT_DOMAINS
    tsv_files = sorted(feat_dir.glob("hmmer_*.tsv"))
    log.info("  Found %d HMMER TSV files", len(tsv_files))

    if not tsv_files:
        log.warning("  No HMMER TSV files found")
        return

    uid_map = dict(conn.execute(
        "SELECT cdb_id, id FROM sequences"
    ).fetchall())
    uniprot_map = dict(conn.execute(
        "SELECT cdb_id, uniprot_id FROM sequences"
    ).fetchall())

    total = 0
    for tsv in tqdm(tsv_files, desc="domains"):
        df = pd.read_csv(tsv, sep="\t")
        if dry_run:
            total += len(df)
            continue

        df["sequence_id"] = df["cdb_id"].map(uid_map)
        df["uniprot_id"]  = df["cdb_id"].map(uniprot_map)
        df = df.dropna(subset=["sequence_id"])
        df["sequence_id"] = df["sequence_id"].astype(int)

        cols = [c for c in df.columns if c != "cdb_id"]
        col_str = ", ".join(cols)
        ph_str  = ", ".join([f":{c}" for c in cols])

        conn.executemany(
            f"INSERT OR IGNORE INTO features_domains ({col_str}) VALUES ({ph_str})",
            df[cols].where(df[cols].notna(), None).to_dict("records")
        )
        conn.commit()
        total += len(df)

    log.info("  features_domains: %d rows", table_count(conn, "features_domains"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 4: InterProScan features
# ═══════════════════════════════════════════════════════════════════════════

def build_interpro(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building features_interpro ──")

    feat_dir = PATHS.FEAT_INTERPRO
    tsv_files = sorted(feat_dir.glob("ipr_*.tsv"))
    existing  = table_count(conn, "features_interpro")
    log.info("  Found %d InterPro TSV files (existing rows: %d)", len(tsv_files), existing)

    if not tsv_files:
        log.warning("  No InterPro TSV files found — skipping")
        return

    uid_map     = dict(conn.execute("SELECT cdb_id, id FROM sequences").fetchall())
    uniprot_map = dict(conn.execute("SELECT cdb_id, uniprot_id FROM sequences").fetchall())
    done_uids   = set(r[0] for r in conn.execute(
        "SELECT uniprot_id FROM features_interpro"
    ).fetchall())

    total = 0
    for tsv in tqdm(tsv_files, desc="interpro"):
        df = pd.read_csv(tsv, sep="\t")
        if dry_run:
            total += len(df)
            continue

        df["sequence_id"] = df["cdb_id"].map(uid_map)
        df["uniprot_id"]  = df["cdb_id"].map(uniprot_map)
        df = df.dropna(subset=["sequence_id"])
        df = df[~df["uniprot_id"].isin(done_uids)]

        if df.empty:
            continue

        df["sequence_id"] = df["sequence_id"].astype(int)
        cols = [c for c in df.columns if c != "cdb_id"]
        col_str = ", ".join(cols)
        ph_str  = ", ".join([f":{c}" for c in cols])

        conn.executemany(
            f"INSERT OR IGNORE INTO features_interpro ({col_str}) VALUES ({ph_str})",
            df[cols].where(df[cols].notna(), None).to_dict("records")
        )
        conn.commit()
        total += len(df)

    log.info("  features_interpro: %d rows", table_count(conn, "features_interpro"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 5: ESM-2 embeddings (store as float32 blob)
# ═══════════════════════════════════════════════════════════════════════════

def build_esm2(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building features_esm2 ──")

    existing = table_count(conn, "features_esm2")
    feat_dir = PATHS.FEAT_ESM2
    tsv_files = sorted(feat_dir.glob("esm2_*.tsv"))
    log.info("  Found %d ESM-2 TSV files (existing rows: %d)", len(tsv_files), existing)

    if not tsv_files:
        log.warning("  No ESM-2 TSV files found — skipping")
        return

    uid_map     = dict(conn.execute("SELECT cdb_id, id FROM sequences").fetchall())
    uniprot_map = dict(conn.execute("SELECT cdb_id, uniprot_id FROM sequences").fetchall())
    done_uids   = set(r[0] for r in conn.execute(
        "SELECT uniprot_id FROM features_esm2"
    ).fetchall())

    ESM2_DIM = 1280
    emb_cols = [f"esm2_{i}" for i in range(ESM2_DIM)]

    total = 0
    for tsv in tqdm(tsv_files, desc="ESM-2"):
        df = pd.read_csv(tsv, sep="\t")
        if dry_run:
            total += len(df)
            continue

        df["sequence_id"] = df["cdb_id"].map(uid_map)
        df["uniprot_id"]  = df["cdb_id"].map(uniprot_map)
        df = df.dropna(subset=["sequence_id"])
        df = df[~df["uniprot_id"].isin(done_uids)]

        if df.empty:
            continue

        df["sequence_id"] = df["sequence_id"].astype(int)

        rows = []
        for _, row in df.iterrows():
            emb = np.array([row[c] for c in emb_cols
                           if c in row], dtype=np.float32)
            rows.append({
                "sequence_id":   int(row["sequence_id"]),
                "uniprot_id":    row["uniprot_id"],
                "embedding_blob": emb.tobytes(),
                "model_version":  "esm2_t33_650M_UR50D",
                "computed_at":    TS,
            })

        conn.executemany("""
            INSERT OR IGNORE INTO features_esm2
            (sequence_id, uniprot_id, embedding_blob, model_version, computed_at)
            VALUES (:sequence_id, :uniprot_id, :embedding_blob, :model_version, :computed_at)
        """, rows)
        conn.commit()
        total += len(df)

    log.info("  features_esm2: %d rows", table_count(conn, "features_esm2"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 6: Ankh embeddings
# ═══════════════════════════════════════════════════════════════════════════

def build_ankh(conn: sqlite3.Connection, dry_run: bool = False):
    log.info("── Building features_ankh ──")

    existing  = table_count(conn, "features_ankh")
    feat_dir  = PATHS.FEAT_ANKH
    tsv_files = sorted(feat_dir.glob("ankh_*.tsv"))
    log.info("  Found %d Ankh TSV files (existing rows: %d)", len(tsv_files), existing)

    if not tsv_files:
        log.warning("  No Ankh TSV files found — skipping")
        return

    uid_map     = dict(conn.execute("SELECT cdb_id, id FROM sequences").fetchall())
    uniprot_map = dict(conn.execute("SELECT cdb_id, uniprot_id FROM sequences").fetchall())
    done_uids   = set(r[0] for r in conn.execute(
        "SELECT uniprot_id FROM features_ankh"
    ).fetchall())

    ANKH_DIM = 1536
    emb_cols = [f"ankh_{i}" for i in range(ANKH_DIM)]

    total = 0
    for tsv in tqdm(tsv_files, desc="Ankh"):
        df = pd.read_csv(tsv, sep="\t")
        if dry_run:
            total += len(df)
            continue

        df["sequence_id"] = df["cdb_id"].map(uid_map)
        df["uniprot_id"]  = df["cdb_id"].map(uniprot_map)
        df = df.dropna(subset=["sequence_id"])
        df = df[~df["uniprot_id"].isin(done_uids)]

        if df.empty:
            continue

        df["sequence_id"] = df["sequence_id"].astype(int)

        rows = []
        for _, row in df.iterrows():
            emb = np.array([row[c] for c in emb_cols
                           if c in row], dtype=np.float32)
            rows.append({
                "sequence_id":    int(row["sequence_id"]),
                "uniprot_id":     row["uniprot_id"],
                "embedding_blob": emb.tobytes(),
                "model_version":  "ankh-large",
                "computed_at":    TS,
            })

        conn.executemany("""
            INSERT OR IGNORE INTO features_ankh
            (sequence_id, uniprot_id, embedding_blob, model_version, computed_at)
            VALUES (:sequence_id, :uniprot_id, :embedding_blob, :model_version, :computed_at)
        """, rows)
        conn.commit()
        total += len(df)

    log.info("  features_ankh: %d rows", table_count(conn, "features_ankh"))


# ═══════════════════════════════════════════════════════════════════════════
# Step 7: Update db_metadata
# ═══════════════════════════════════════════════════════════════════════════

def update_metadata(conn: sqlite3.Connection):
    log.info("── Updating db_metadata ──")
    stats = {
        "n_sequences":          table_count(conn, "sequences"),
        "n_label1":             conn.execute("SELECT COUNT(*) FROM sequences WHERE label=1").fetchone()[0],
        "n_label2":             conn.execute("SELECT COUNT(*) FROM sequences WHERE label=2").fetchone()[0],
        "n_label0":             conn.execute("SELECT COUNT(*) FROM sequences WHERE label=0").fetchone()[0],
        "n_km_gold":            conn.execute("SELECT COUNT(*) FROM sequences WHERE km_best_mM IS NOT NULL").fetchone()[0],
        "n_features_comp":      table_count(conn, "features_composition"),
        "n_features_domains":   table_count(conn, "features_domains"),
        "n_features_interpro":  table_count(conn, "features_interpro"),
        "n_features_esm2":      table_count(conn, "features_esm2"),
        "n_features_ankh":      table_count(conn, "features_ankh"),
        "built_at":             TS,
    }
    for k, v in stats.items():
        conn.execute(
            "INSERT OR REPLACE INTO db_metadata (key, value) VALUES (?, ?)",
            (k, str(v))
        )
    conn.commit()
    log.info("  Metadata updated")
    for k, v in stats.items():
        log.info("    %-30s %s", k, v)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def print_summary(conn: sqlite3.Connection):
    log.info("\n" + "=" * 70)
    log.info("DATABASE SUMMARY")
    log.info("=" * 70)
    tables = [
        "sequences", "id_map", "ec_evidence", "km_evidence",
        "features_composition", "features_domains", "features_interpro",
        "features_esm2", "features_ankh", "features_blast",
        "features_expert_motifs", "predictions", "confidence_scores",
    ]
    for t in tables:
        n = table_count(conn, t)
        log.info("  %-30s %12d rows", t, n)

    size_mb = DB_PATH.stat().st_size / 1e6
    log.info("\n  DB size: %.1f MB (%.1f GB)", size_mb, size_mb / 1000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--update",  default=None,
                    help="Update only one table: interpro|ankh|esm2|domains|composition")
    args = ap.parse_args()

    if args.dry_run:
        log.info("DRY RUN — no database writes")

    PATHS.PRIMARY.mkdir(parents=True, exist_ok=True)

    log.info("Database path: %s", DB_PATH)
    log.info("Schema:        %s", SCHEMA_SQL)

    conn = get_connection(DB_PATH)

    t0 = time.time()

    if args.update:
        # Incremental update mode
        log.info("Incremental update: %s", args.update)
        apply_schema(conn)  # safe — uses CREATE IF NOT EXISTS
        if args.update == "interpro":
            build_interpro(conn, args.dry_run)
        elif args.update == "ankh":
            build_ankh(conn, args.dry_run)
        elif args.update == "esm2":
            build_esm2(conn, args.dry_run)
        elif args.update == "domains":
            build_domains(conn, args.dry_run)
        elif args.update == "composition":
            build_composition_v2(conn, args.dry_run)
        update_metadata(conn)
    else:
        # Full build
        apply_schema(conn)
        build_sequences(conn, args.dry_run)
        build_composition_v2(conn, args.dry_run)
        build_domains(conn, args.dry_run)
        build_interpro(conn, args.dry_run)
        build_esm2(conn, args.dry_run)
        build_ankh(conn, args.dry_run)
        update_metadata(conn)

    elapsed = time.time() - t0
    log.info("\nTotal time: %.1f min", elapsed / 60)

    if not args.dry_run:
        print_summary(conn)

    conn.close()
    log.info("Done. Next: python scripts/06_cluster_sequences.py")


if __name__ == "__main__":
    main()
