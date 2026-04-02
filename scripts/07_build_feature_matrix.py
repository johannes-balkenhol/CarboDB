#!/usr/bin/env python3
"""
07_build_feature_matrix.py
==========================
CarboDB — Step 07: Assemble ML feature matrices from carbodb.sqlite.

Builds three feature matrices (one per ML task):
  X_binary:  composition + domains + ESM-2  → binary classification
  X_ec:      composition + domains + ESM-2  → EC class prediction
  X_km:      composition + domains + ESM-2  → Km regression

Features used (v5 model):
  - Composition:  489 features (aac + physico + catalytic + motifs + dp + pse)
  - Pfam domains: 19 features (18 binary + n_hits)
  - ESM-2:       1280 features (mean-pooled embeddings, float32 blob)
  Total:         1,788 features

InterPro + Ankh added automatically when available (v6 model, ~3,350 features).

Output (numpy .npz files for fast loading):
  data/ml/X_binary_train.npz   X_binary_val.npz   X_binary_test.npz
  data/ml/y_binary_train.npy   y_binary_val.npy   y_binary_test.npy
  data/ml/X_ec_train.npz       X_ec_val.npz       X_ec_test.npz
  data/ml/y_ec_train.npy       y_ec_val.npy       y_ec_test.npy
  data/ml/X_km_train.npz       X_km_val.npz       X_km_test.npz
  data/ml/y_km_train.npy       y_km_val.npy       y_km_test.npy
  data/ml/feature_names.json   — ordered list of feature names
  data/ml/ec_label_map.json    — EC string → integer label
  data/ml/matrix_summary.json  — shapes, feature counts, stats

Usage:
  python scripts/07_build_feature_matrix.py
  python scripts/07_build_feature_matrix.py --tasks binary ec km
  python scripts/07_build_feature_matrix.py --no-embeddings  (composition+domains only)
"""

import argparse
import json
import logging
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("07_build_feature_matrix")

DB_PATH    = PATHS.PRIMARY / "carbodb.sqlite"
SPLIT_DIR  = PATHS.PRIMARY.parent / "splits"
ML_DIR     = PATHS.PRIMARY.parent / "ml"

# Composition columns to extract (matches schema exactly)
AAC_COLS  = [f"aac_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"]
PHYS_COLS = [
    "phys_length", "phys_length_log", "phys_mw", "phys_pi",
    "phys_gravy", "phys_aromaticity", "phys_instability",
    "phys_charge_ph7", "phys_frac_charged", "phys_frac_aromatic",
    "phys_frac_polar", "phys_frac_nonpolar", "phys_frac_small",
    "phys_frac_glycine", "phys_frac_proline",
]
INV_COLS  = [
    "inv_cat_D", "inv_cat_E", "inv_cat_H", "inv_cat_K",
    "inv_cat_C", "inv_cat_S", "inv_cat_T",
    "inv_cat_mean_dist", "inv_cat_std_dist",
    "inv_cat_min_dist", "inv_cat_max_dist", "inv_cat_clustering",
    "inv_hydrophobic", "inv_charged", "inv_polar",
    "inv_aromatic", "inv_net_charge",
]
MOTIF_COLS = [
    "motif_rubisco_kk", "motif_rubisco_gk",
    "motif_ca_hh", "motif_ca_his_cluster",
    "motif_pepc_rr", "motif_biotin_mk", "motif_biotin_amk",
]
DP_COLS   = [f"dp_{a1}{a2}"
             for a1 in "ACDEFGHIKLMNPQRSTVWY"
             for a2 in "ACDEFGHIKLMNPQRSTVWY"]
PSE_COLS  = ([f"pse_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"] +
             [f"pse_corr_{i}" for i in range(1, 11)])

COMP_COLS = AAC_COLS + PHYS_COLS + INV_COLS + MOTIF_COLS + DP_COLS + PSE_COLS

# Pfam domain columns
PFAM_COLS = [
    "pfam_PF00016", "pfam_PF02788", "pfam_PF00101", "pfam_PF00194",
    "pfam_PF03119", "pfam_PF00311", "pfam_PF00821", "pfam_PF02785",
    "pfam_PF00364", "pfam_PF01039", "pfam_PF02786", "pfam_PF02787",
    "pfam_PF00289", "pfam_PF01309", "pfam_PF03599", "pfam_PF03590",
    "pfam_PF00384", "pfam_PF00682", "pfam_n_hits",
]

# InterPro scalar columns (optional, v6)
IPR_COLS = [
    "n_panther", "n_gene3d", "n_tigrfam", "n_prosite_prof", "n_prosite_pat",
]

ESM2_DIM = 1280
ANKH_DIM = 1536


# ═══════════════════════════════════════════════════════════════════════════
# Database helpers
# ═══════════════════════════════════════════════════════════════════════════

def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA cache_size = -2000000")
    conn.execute("PRAGMA temp_store = MEMORY")
    return conn


def check_table_coverage(conn: sqlite3.Connection) -> dict:
    """Report how many sequences have each feature table filled."""
    tables = [
        "features_composition", "features_domains",
        "features_interpro", "features_esm2", "features_ankh",
    ]
    total = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
    coverage = {"total_sequences": total}
    for t in tables:
        try:
            n = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
            coverage[t] = {"n": n, "pct": round(100 * n / total, 1)}
        except Exception:
            coverage[t] = {"n": 0, "pct": 0.0}
    return coverage


# ═══════════════════════════════════════════════════════════════════════════
# Feature loading
# ═══════════════════════════════════════════════════════════════════════════

def load_composition(conn: sqlite3.Connection, cdb_ids: list) -> pd.DataFrame:
    """Load composition features for given cdb_ids."""
    log.info("  Loading composition features (%d cols)...", len(COMP_COLS))

    cols_str = ", ".join(["c.uniprot_id"] + [f"c.{col}" for col in COMP_COLS])

    # Load in chunks of 50k
    results = []
    chunk_size = 900
    for i in tqdm(range(0, len(cdb_ids), chunk_size), desc="comp chunks"):
        chunk = cdb_ids[i:i+chunk_size]
        placeholders = ",".join(["?" * len(chunk)])
        query = f"""
            SELECT {cols_str}
            FROM features_composition c
            JOIN sequences s ON s.id = c.sequence_id
            WHERE s.cdb_id IN ({placeholders})
        """
        df = pd.read_sql_query(query, conn, params=chunk)
        results.append(df)

    if not results:
        return pd.DataFrame(columns=["uniprot_id"] + COMP_COLS)

    df = pd.concat(results, ignore_index=True)
    log.info("  Composition: %d rows loaded", len(df))
    return df


def load_domains(conn: sqlite3.Connection, cdb_ids: list) -> pd.DataFrame:
    """Load Pfam domain features."""
    log.info("  Loading domain features (%d cols)...", len(PFAM_COLS))

    cols_str = ", ".join(["s.uniprot_id"] + [f"d.{col}" for col in PFAM_COLS])
    results = []
    chunk_size = 900

    for i in tqdm(range(0, len(cdb_ids), chunk_size), desc="domain chunks"):
        chunk = cdb_ids[i:i+chunk_size]
        placeholders = ",".join(["?"] * len(chunk))
        query = f"""
            SELECT {cols_str}
            FROM features_domains d
            JOIN sequences s ON s.id = d.sequence_id
            WHERE s.cdb_id IN ({placeholders})
        """
        df = pd.read_sql_query(query, conn, params=chunk)
        results.append(df)

    if not results:
        return pd.DataFrame(columns=["uniprot_id"] + PFAM_COLS)

    df = pd.concat(results, ignore_index=True)
    log.info("  Domains: %d rows loaded", len(df))
    return df


def load_interpro(conn: sqlite3.Connection, cdb_ids: list) -> pd.DataFrame:
    """Load InterPro scalar features (optional)."""
    try:
        n = conn.execute("SELECT COUNT(*) FROM features_interpro").fetchone()[0]
        if n == 0:
            log.info("  InterPro: no data yet — skipping")
            return pd.DataFrame(columns=["uniprot_id"] + IPR_COLS)
    except Exception:
        return pd.DataFrame(columns=["uniprot_id"] + IPR_COLS)

    log.info("  Loading InterPro features (%d cols, %d rows)...", len(IPR_COLS), n)
    cols_str = ", ".join(["s.uniprot_id"] + [f"i.{col}" for col in IPR_COLS])
    results = []
    chunk_size = 900

    for i in tqdm(range(0, len(cdb_ids), chunk_size), desc="ipr chunks"):
        chunk = cdb_ids[i:i+chunk_size]
        placeholders = ",".join(["?"] * len(chunk))
        query = f"""
            SELECT {cols_str}
            FROM features_interpro i
            JOIN sequences s ON s.id = i.sequence_id
            WHERE s.cdb_id IN ({placeholders})
        """
        df = pd.read_sql_query(query, conn, params=chunk)
        results.append(df)

    df = pd.concat(results, ignore_index=True) if results else pd.DataFrame(
        columns=["uniprot_id"] + IPR_COLS)
    log.info("  InterPro: %d rows loaded", len(df))
    return df


def load_esm2_embeddings(conn: sqlite3.Connection, cdb_ids: list) -> pd.DataFrame:
    """Load ESM-2 embeddings from blob storage → float32 matrix."""
    log.info("  Loading ESM-2 embeddings (1280-dim)...")

    results = []
    chunk_size = 900  # smaller chunks — each row is 5KB

    for i in tqdm(range(0, len(cdb_ids), chunk_size), desc="ESM-2 chunks"):
        chunk = cdb_ids[i:i+chunk_size]
        placeholders = ",".join(["?"] * len(chunk))
        query = f"""
            SELECT s.uniprot_id, e.embedding_blob
            FROM features_esm2 e
            JOIN sequences s ON s.id = e.sequence_id
            WHERE s.cdb_id IN ({placeholders})
        """
        rows = conn.execute(query, chunk).fetchall()
        for uid, blob in rows:
            if blob:
                emb = np.frombuffer(blob, dtype=np.float32)
                if len(emb) == ESM2_DIM:
                    results.append({"uniprot_id": uid, **{
                        f"esm2_{j}": float(emb[j]) for j in range(ESM2_DIM)
                    }})

    if not results:
        log.warning("  ESM-2: no embeddings found")
        return pd.DataFrame(columns=["uniprot_id"] +
                           [f"esm2_{j}" for j in range(ESM2_DIM)])

    df = pd.DataFrame(results)
    log.info("  ESM-2: %d embeddings loaded", len(df))
    return df


def load_ankh_embeddings(conn: sqlite3.Connection, cdb_ids: list) -> pd.DataFrame:
    """Load Ankh embeddings (optional, partial)."""
    try:
        n = conn.execute("SELECT COUNT(*) FROM features_ankh").fetchone()[0]
        if n == 0:
            log.info("  Ankh: no data yet — skipping")
            return pd.DataFrame(columns=["uniprot_id"])
    except Exception:
        return pd.DataFrame(columns=["uniprot_id"])

    log.info("  Loading Ankh embeddings (1536-dim, %d rows)...", n)
    results = []
    chunk_size = 900

    for i in tqdm(range(0, len(cdb_ids), chunk_size), desc="Ankh chunks"):
        chunk = cdb_ids[i:i+chunk_size]
        placeholders = ",".join(["?"] * len(chunk))
        query = f"""
            SELECT s.uniprot_id, a.embedding_blob
            FROM features_ankh a
            JOIN sequences s ON s.id = a.sequence_id
            WHERE s.cdb_id IN ({placeholders})
        """
        rows = conn.execute(query, chunk).fetchall()
        for uid, blob in rows:
            if blob:
                emb = np.frombuffer(blob, dtype=np.float32)
                if len(emb) == ANKH_DIM:
                    results.append({"uniprot_id": uid, **{
                        f"ankh_{j}": float(emb[j]) for j in range(ANKH_DIM)
                    }})

    df = pd.DataFrame(results) if results else pd.DataFrame(
        columns=["uniprot_id"] + [f"ankh_{j}" for j in range(ANKH_DIM)])
    log.info("  Ankh: %d embeddings loaded", len(df))
    return df


# ═══════════════════════════════════════════════════════════════════════════
# Matrix assembly
# ═══════════════════════════════════════════════════════════════════════════

def assemble_matrix(
    split_df: pd.DataFrame,
    conn: sqlite3.Connection,
    use_embeddings: bool = True,
    use_interpro: bool = True,
    use_ankh: bool = False,
) -> tuple:
    """
    Assemble feature matrix X and target vector y for a given split DataFrame.
    Returns (X, y, feature_names, uid_index).
    """
    cdb_ids = split_df["cdb_id"].tolist()
    log.info("  Assembling matrix for %d sequences...", len(cdb_ids))

    # Load all feature tables
    comp   = load_composition(conn, cdb_ids)
    domain = load_domains(conn, cdb_ids)
    ipr    = load_interpro(conn, cdb_ids) if use_interpro else None
    esm2   = load_esm2_embeddings(conn, cdb_ids) if use_embeddings else None
    ankh   = load_ankh_embeddings(conn, cdb_ids) if use_ankh else None

    # Merge all on uniprot_id
    # Start from split_df to preserve order and include all sequences
    merged = split_df[["cdb_id", "uniprot_id"]].copy()

    # Get uniprot_id from sequences table for cdb_id mapping
    uid_map = dict(conn.execute(
        "SELECT cdb_id, uniprot_id FROM sequences"
    ).fetchall())
    merged["uniprot_id"] = merged["cdb_id"].map(uid_map)

    merged = merged.merge(comp,   on="uniprot_id", how="left")
    merged = merged.merge(domain, on="uniprot_id", how="left")

    if ipr is not None and len(ipr) > 0 and "n_panther" in ipr.columns:
        merged = merged.merge(ipr, on="uniprot_id", how="left")

    if esm2 is not None and len(esm2) > 0:
        merged = merged.merge(esm2, on="uniprot_id", how="left")

    if ankh is not None and len(ankh) > 0:
        merged = merged.merge(ankh, on="uniprot_id", how="left")

    # Determine feature columns (everything except id cols)
    id_cols = {"cdb_id", "uniprot_id"}
    feat_cols = [c for c in merged.columns if c not in id_cols]

    log.info("  Feature columns: %d", len(feat_cols))
    log.info("  Missing values: %.1f%%",
             100 * merged[feat_cols].isna().sum().sum() /
             (len(merged) * len(feat_cols)))

    # Fill missing with 0 (mean imputation alternative, but 0 is safe for binary/count features)
    X = merged[feat_cols].fillna(0.0).values.astype(np.float32)
    uid_index = merged["uniprot_id"].values

    return X, feat_cols, uid_index


# ═══════════════════════════════════════════════════════════════════════════
# Main build functions
# ═══════════════════════════════════════════════════════════════════════════

def build_binary_matrices(conn: sqlite3.Connection,
                          use_embeddings: bool, dry_run: bool):
    log.info("══ Building BINARY classification matrices ══")
    split_file = SPLIT_DIR / "split_binary.tsv"
    if not split_file.exists():
        log.error("Split file not found: %s — run script 06 first", split_file)
        return

    splits = pd.read_csv(split_file, sep="\t")
    log.info("Total: %d  positives: %d  negatives: %d",
             len(splits), (splits["label"]==1).sum(), (splits["label"]==0).sum())

    for split_name in ["train", "val", "test"]:
        sub = splits[splits["split"] == split_name].copy()
        log.info("── %s: %d sequences ──", split_name, len(sub))

        X, feat_cols, uids = assemble_matrix(sub, conn, use_embeddings)
        y = sub["label"].values.astype(np.int32)

        if not dry_run:
            np.savez_compressed(ML_DIR / f"X_binary_{split_name}.npz", X=X)
            np.save(ML_DIR / f"y_binary_{split_name}.npy", y)
            np.save(ML_DIR / f"uid_binary_{split_name}.npy", uids)
            log.info("  Saved X_%s: %s  y_%s: %s",
                     split_name, X.shape, split_name, y.shape)

    # Save feature names once
    if not dry_run:
        with open(ML_DIR / "feature_names_binary.json", "w") as f:
            json.dump(feat_cols, f)
        log.info("  Feature names saved: %d features", len(feat_cols))


def build_ec_matrices(conn: sqlite3.Connection,
                      use_embeddings: bool, dry_run: bool):
    log.info("══ Building EC CLASS PREDICTION matrices ══")
    split_file = SPLIT_DIR / "split_ec.tsv"
    if not split_file.exists():
        log.error("Split file not found: %s — run script 06 first", split_file)
        return

    splits = pd.read_csv(split_file, sep="\t")

    # Build EC → integer label map
    ec_classes = sorted(splits["ec_number"].unique())
    ec_label_map = {ec: i for i, ec in enumerate(ec_classes)}
    log.info("EC classes: %d", len(ec_classes))

    if not dry_run:
        with open(ML_DIR / "ec_label_map.json", "w") as f:
            json.dump(ec_label_map, f, indent=2)

    for split_name in ["train", "val", "test"]:
        sub = splits[splits["split"] == split_name].copy()
        log.info("── %s: %d sequences ──", split_name, len(sub))

        X, feat_cols, uids = assemble_matrix(sub, conn, use_embeddings)
        y = sub["ec_number"].map(ec_label_map).values.astype(np.int32)

        if not dry_run:
            np.savez_compressed(ML_DIR / f"X_ec_{split_name}.npz", X=X)
            np.save(ML_DIR / f"y_ec_{split_name}.npy", y)
            np.save(ML_DIR / f"uid_ec_{split_name}.npy", uids)
            log.info("  Saved X_%s: %s  y_%s: %s",
                     split_name, X.shape, split_name, y.shape)

    if not dry_run:
        with open(ML_DIR / "feature_names_ec.json", "w") as f:
            json.dump(feat_cols, f)


def build_km_matrices(conn: sqlite3.Connection,
                      use_embeddings: bool, dry_run: bool):
    log.info("══ Building Km REGRESSION matrices ══")
    split_file = SPLIT_DIR / "split_km.tsv"
    if not split_file.exists():
        log.error("Split file not found: %s — run script 06 first", split_file)
        return

    splits = pd.read_csv(split_file, sep="\t")
    log.info("Total Km sequences: %d  log10 range: %.2f to %.2f",
             len(splits), splits["km_log10_mM"].min(), splits["km_log10_mM"].max())

    for split_name in ["train", "val", "test"]:
        sub = splits[splits["split"] == split_name].copy()
        log.info("── %s: %d sequences ──", split_name, len(sub))

        X, feat_cols, uids = assemble_matrix(sub, conn, use_embeddings)
        y = sub["km_log10_mM"].values.astype(np.float32)

        if not dry_run:
            np.savez_compressed(ML_DIR / f"X_km_{split_name}.npz", X=X)
            np.save(ML_DIR / f"y_km_{split_name}.npy", y)
            np.save(ML_DIR / f"uid_km_{split_name}.npy", uids)
            log.info("  Saved X_%s: %s  y_%s: %s  (log10 Km)",
                     split_name, X.shape, split_name, y.shape)

    if not dry_run:
        with open(ML_DIR / "feature_names_km.json", "w") as f:
            json.dump(feat_cols, f)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", nargs="+", default=["binary", "ec", "km"],
                    choices=["binary", "ec", "km"])
    ap.add_argument("--no-embeddings", action="store_true",
                    help="Skip ESM-2/Ankh embeddings (composition+domains only)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not DB_PATH.exists():
        log.error("Database not found: %s — run script 05 first", DB_PATH)
        sys.exit(1)

    ML_DIR.mkdir(parents=True, exist_ok=True)

    conn = get_conn()

    # Report feature coverage
    log.info("Feature coverage in database:")
    coverage = check_table_coverage(conn)
    for k, v in coverage.items():
        if isinstance(v, dict):
            log.info("  %-30s %d rows (%.1f%%)", k, v["n"], v["pct"])
        else:
            log.info("  %-30s %d", k, v)

    use_embeddings = not args.no_embeddings

    if "binary" in args.tasks:
        build_binary_matrices(conn, use_embeddings, args.dry_run)

    if "ec" in args.tasks:
        build_ec_matrices(conn, use_embeddings, args.dry_run)

    if "km" in args.tasks:
        build_km_matrices(conn, use_embeddings, args.dry_run)

    # Save summary
    if not args.dry_run:
        summary = {
            "created_at":     TS,
            "model_version":  "v5",
            "use_embeddings": use_embeddings,
            "features": {
                "composition": len(COMP_COLS),
                "domains":     len(PFAM_COLS),
                "esm2":        ESM2_DIM if use_embeddings else 0,
                "interpro":    len(IPR_COLS),
                "ankh":        0,  # added in v6
            },
            "tasks": args.tasks,
        }
        with open(ML_DIR / "matrix_summary.json", "w") as f:
            json.dump(summary, f, indent=2)
        log.info("Saved matrix summary: %s", ML_DIR / "matrix_summary.json")

    conn.close()
    log.info("Done. Next: python scripts/08_train_models.py")


if __name__ == "__main__":
    main()
