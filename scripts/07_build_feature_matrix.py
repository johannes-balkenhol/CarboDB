#!/usr/bin/env python3
"""
07_build_feature_matrix.py  v2
================================
CarboDB — Step 07: Assemble ML feature matrices from carbodb.sqlite.
Uses temp table approach to avoid SQLite 999 parameter limit.

Features (v5 model):
  Composition:  489 (aac + physico + catalytic + motifs + dipeptide + pse)
  Pfam domains:  19 (18 binary + n_hits)
  ESM-2:       1280 (mean-pooled float32 blob)
  Total:       1,788

Output: data/ml/  X_*.npz  y_*.npy  uid_*.npy  feature_names_*.json
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

DB_PATH   = PATHS.PRIMARY / "carbodb.sqlite"
SPLIT_DIR = PATHS.PRIMARY.parent / "splits"
ML_DIR    = PATHS.PRIMARY.parent / "ml"

# ── Feature column definitions ────────────────────────────────────────────
AAC_COLS   = [f"aac_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"]
PHYS_COLS  = ["phys_length","phys_length_log","phys_mw","phys_pi",
               "phys_gravy","phys_aromaticity","phys_instability",
               "phys_charge_ph7","phys_frac_charged","phys_frac_aromatic",
               "phys_frac_polar","phys_frac_nonpolar","phys_frac_small",
               "phys_frac_glycine","phys_frac_proline"]
INV_COLS   = ["inv_cat_D","inv_cat_E","inv_cat_H","inv_cat_K",
               "inv_cat_C","inv_cat_S","inv_cat_T",
               "inv_cat_mean_dist","inv_cat_std_dist",
               "inv_cat_min_dist","inv_cat_max_dist","inv_cat_clustering",
               "inv_hydrophobic","inv_charged","inv_polar",
               "inv_aromatic","inv_net_charge"]
MOTIF_COLS = ["motif_rubisco_kk","motif_rubisco_gk",
               "motif_ca_hh","motif_ca_his_cluster",
               "motif_pepc_rr","motif_biotin_mk","motif_biotin_amk"]
DP_COLS    = [f"dp_{a1}{a2}"
              for a1 in "ACDEFGHIKLMNPQRSTVWY"
              for a2 in "ACDEFGHIKLMNPQRSTVWY"]
PSE_COLS   = ([f"pse_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"] +
              [f"pse_corr_{i}" for i in range(1, 11)])
COMP_COLS  = AAC_COLS + PHYS_COLS + INV_COLS + MOTIF_COLS + DP_COLS + PSE_COLS

PFAM_COLS  = ["pfam_PF00016","pfam_PF02788","pfam_PF00101","pfam_PF00194",
               "pfam_PF03119","pfam_PF00311","pfam_PF00821","pfam_PF02785",
               "pfam_PF00364","pfam_PF01039","pfam_PF02786","pfam_PF02787",
               "pfam_PF00289","pfam_PF01309","pfam_PF03599","pfam_PF03590",
               "pfam_PF00384","pfam_PF00682","pfam_n_hits"]

IPR_COLS   = ["n_panther","n_gene3d","n_tigrfam","n_prosite_prof","n_prosite_pat"]

ESM2_DIM   = 1280
ANKH_DIM   = 1536


# ── DB helpers ────────────────────────────────────────────────────────────

def get_conn():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA cache_size = -2000000")
    conn.execute("PRAGMA temp_store = MEMORY")
    return conn


def make_tmp(conn, cdb_ids):
    """Load cdb_ids into a temp table — avoids 999 param limit."""
    conn.execute("DROP TABLE IF EXISTS _tmp_ids")
    conn.execute("CREATE TEMP TABLE _tmp_ids (cdb_id TEXT)")
    conn.executemany("INSERT INTO _tmp_ids VALUES (?)", [(x,) for x in cdb_ids])


def check_coverage(conn):
    total = conn.execute("SELECT COUNT(*) FROM sequences").fetchone()[0]
    for t in ["features_composition","features_domains",
               "features_interpro","features_esm2","features_ankh"]:
        try:
            n = conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
            log.info("  %-30s %d / %d (%.0f%%)", t, n, total, 100*n/total)
        except Exception:
            log.info("  %-30s not available", t)


# ── Feature loaders ───────────────────────────────────────────────────────

def load_composition(conn, cdb_ids):
    log.info("  Loading composition (%d features)...", len(COMP_COLS))
    make_tmp(conn, cdb_ids)
    cols = ", ".join([f"c.{c}" for c in COMP_COLS])
    sql  = f"""
        SELECT s.uniprot_id, {cols}
        FROM features_composition c
        JOIN sequences s ON s.id = c.sequence_id
        JOIN _tmp_ids t  ON t.cdb_id = s.cdb_id
    """
    df = pd.read_sql_query(sql, conn)
    log.info("  Composition: %d rows", len(df))
    return df


def load_domains(conn, cdb_ids):
    log.info("  Loading domains (%d features)...", len(PFAM_COLS))
    make_tmp(conn, cdb_ids)
    cols = ", ".join([f"d.{c}" for c in PFAM_COLS])
    sql  = f"""
        SELECT s.uniprot_id, {cols}
        FROM features_domains d
        JOIN sequences s ON s.id = d.sequence_id
        JOIN _tmp_ids t  ON t.cdb_id = s.cdb_id
    """
    df = pd.read_sql_query(sql, conn)
    log.info("  Domains: %d rows", len(df))
    return df


def load_interpro(conn, cdb_ids):
    try:
        n = conn.execute("SELECT COUNT(*) FROM features_interpro").fetchone()[0]
    except Exception:
        n = 0
    if n == 0:
        log.info("  InterPro: no data yet — skipping")
        return pd.DataFrame(columns=["uniprot_id"] + IPR_COLS)

    log.info("  Loading InterPro (%d scalar features, %d rows in DB)...",
             len(IPR_COLS), n)
    make_tmp(conn, cdb_ids)
    cols = ", ".join([f"i.{c}" for c in IPR_COLS])
    sql  = f"""
        SELECT s.uniprot_id, {cols}
        FROM features_interpro i
        JOIN sequences s ON s.id = i.sequence_id
        JOIN _tmp_ids t  ON t.cdb_id = s.cdb_id
    """
    df = pd.read_sql_query(sql, conn)
    log.info("  InterPro: %d rows loaded", len(df))
    return df


def load_esm2(conn, cdb_ids):
    try:
        n = conn.execute("SELECT COUNT(*) FROM features_esm2").fetchone()[0]
    except Exception:
        n = 0
    if n == 0:
        log.warning("  ESM-2: no data — skipping")
        return pd.DataFrame(columns=["uniprot_id"])

    log.info("  Loading ESM-2 blobs (%d rows in DB)...", n)
    make_tmp(conn, cdb_ids)
    rows = conn.execute("""
        SELECT s.uniprot_id, e.embedding_blob
        FROM features_esm2 e
        JOIN sequences s ON s.id = e.sequence_id
        JOIN _tmp_ids t  ON t.cdb_id = s.cdb_id
    """).fetchall()

    records = []
    emb_cols = [f"esm2_{j}" for j in range(ESM2_DIM)]
    for uid, blob in tqdm(rows, desc="ESM-2 decode"):
        if blob:
            emb = np.frombuffer(blob, dtype=np.float32)
            if len(emb) == ESM2_DIM:
                r = {"uniprot_id": uid}
                for j, v in enumerate(emb):
                    r[f"esm2_{j}"] = float(v)
                records.append(r)

    df = pd.DataFrame(records) if records else pd.DataFrame(
        columns=["uniprot_id"] + emb_cols)
    log.info("  ESM-2: %d embeddings loaded", len(df))
    return df


def load_ankh(conn, cdb_ids):
    try:
        n = conn.execute("SELECT COUNT(*) FROM features_ankh").fetchone()[0]
    except Exception:
        n = 0
    if n == 0:
        log.info("  Ankh: no data yet — skipping")
        return pd.DataFrame(columns=["uniprot_id"])

    log.info("  Loading Ankh blobs (%d rows in DB)...", n)
    make_tmp(conn, cdb_ids)
    rows = conn.execute("""
        SELECT s.uniprot_id, a.embedding_blob
        FROM features_ankh a
        JOIN sequences s ON s.id = a.sequence_id
        JOIN _tmp_ids t  ON t.cdb_id = s.cdb_id
    """).fetchall()

    records = []
    for uid, blob in tqdm(rows, desc="Ankh decode"):
        if blob:
            emb = np.frombuffer(blob, dtype=np.float32)
            if len(emb) == ANKH_DIM:
                r = {"uniprot_id": uid}
                for j, v in enumerate(emb):
                    r[f"ankh_{j}"] = float(v)
                records.append(r)

    df = pd.DataFrame(records) if records else pd.DataFrame(columns=["uniprot_id"])
    log.info("  Ankh: %d embeddings loaded", len(df))
    return df


# ── Matrix assembly ───────────────────────────────────────────────────────

def assemble(split_df, conn, use_esm2=True, use_interpro=True, use_ankh=False):
    cdb_ids = split_df["cdb_id"].tolist()
    log.info("  Assembling matrix for %d sequences...", len(cdb_ids))

    uid_map = dict(conn.execute(
        "SELECT cdb_id, uniprot_id FROM sequences").fetchall())

    base = split_df[["cdb_id"]].copy()
    base["uniprot_id"] = base["cdb_id"].map(uid_map)

    comp   = load_composition(conn, cdb_ids)
    domain = load_domains(conn, cdb_ids)
    ipr    = load_interpro(conn, cdb_ids) if use_interpro else None
    esm2   = load_esm2(conn, cdb_ids)    if use_esm2    else None
    ankh   = load_ankh(conn, cdb_ids)    if use_ankh    else None

    merged = base.merge(comp,   on="uniprot_id", how="left")
    merged = merged.merge(domain, on="uniprot_id", how="left")

    if ipr is not None and len(ipr) > 1:
        merged = merged.merge(ipr, on="uniprot_id", how="left")

    if esm2 is not None and len(esm2) > 1:
        merged = merged.merge(esm2, on="uniprot_id", how="left")

    if ankh is not None and len(ankh) > 1:
        merged = merged.merge(ankh, on="uniprot_id", how="left")

    feat_cols = [c for c in merged.columns if c not in ("cdb_id","uniprot_id")]
    missing_pct = 100 * merged[feat_cols].isna().sum().sum() / (len(merged)*len(feat_cols))
    log.info("  Features: %d  Missing: %.1f%%", len(feat_cols), missing_pct)

    X    = merged[feat_cols].fillna(0.0).values.astype(np.float32)
    uids = merged["uniprot_id"].values
    return X, feat_cols, uids


# ── Build tasks ───────────────────────────────────────────────────────────

def build_task(task, conn, use_esm2, dry_run):
    split_file = SPLIT_DIR / f"split_{task}.tsv"
    if not split_file.exists():
        log.error("Missing split file: %s — run script 06 first", split_file)
        return None

    splits = pd.read_csv(split_file, sep="\t")
    log.info("══ Task: %s — %d sequences ══", task, len(splits))

    feat_cols_saved = None
    results = {}

    for split_name in ["train", "val", "test"]:
        sub = splits[splits["split"] == split_name].copy()
        log.info("── %s: %d sequences ──", split_name, len(sub))

        X, feat_cols, uids = assemble(sub, conn, use_esm2=use_esm2)

        if task == "binary":
            y = sub["label"].values.astype(np.int32)
        elif task == "ec":
            if feat_cols_saved is None:
                ec_classes = sorted(splits["ec_number"].unique())
                ec_map = {ec: i for i, ec in enumerate(ec_classes)}
                if not dry_run:
                    with open(ML_DIR / "ec_label_map.json", "w") as f:
                        json.dump(ec_map, f, indent=2)
            y = sub["ec_number"].map(ec_map).values.astype(np.int32)
        elif task == "km":
            y = sub["km_log10_mM"].values.astype(np.float32)

        results[split_name] = (X, y, uids)
        feat_cols_saved = feat_cols

        if not dry_run:
            np.savez_compressed(ML_DIR / f"X_{task}_{split_name}.npz", X=X)
            np.save(ML_DIR / f"y_{task}_{split_name}.npy", y)
            np.save(ML_DIR / f"uid_{task}_{split_name}.npy", uids)
            log.info("  Saved X_%s_%s: %s  y: %s",
                     task, split_name, X.shape, y.shape)

    if not dry_run and feat_cols_saved:
        with open(ML_DIR / f"feature_names_{task}.json", "w") as f:
            json.dump(feat_cols_saved, f)
        log.info("  Feature names saved: %d", len(feat_cols_saved))

    return feat_cols_saved


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", nargs="+", default=["binary","ec","km"],
                    choices=["binary","ec","km"])
    ap.add_argument("--no-embeddings", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not DB_PATH.exists():
        log.error("DB not found: %s", DB_PATH)
        sys.exit(1)

    ML_DIR.mkdir(parents=True, exist_ok=True)
    conn = get_conn()

    log.info("Feature coverage:")
    check_coverage(conn)

    use_esm2 = not args.no_embeddings

    feat_cols = None
    for task in args.tasks:
        feat_cols = build_task(task, conn, use_esm2, args.dry_run)

    if not args.dry_run and feat_cols:
        summary = {
            "created_at":   TS,
            "model_version": "v5",
            "n_features":   len(feat_cols),
            "use_esm2":     use_esm2,
            "tasks":        args.tasks,
        }
        with open(ML_DIR / "matrix_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

    conn.close()
    log.info("Done. Next: python scripts/08_train_models.py")


if __name__ == "__main__":
    main()
