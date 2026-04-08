#!/usr/bin/env python3
"""
10_predict_all.py
=================
CarboDB — Step 10: Run all three v5 models on all 2,380,446 sequences.
Writes predictions to the carbodb.sqlite predictions table.
Also computes confidence scores and fills confidence_scores table.

Pipeline per sequence:
  1. Binary classifier  → is_co2_pred + co2_prob
  2. EC classifier      → ec_pred + ec_prob  (all sequences)
  3. Km regressor       → km_pred_mM         (predicted carboxylases only,
                                               trainable EC classes only)
  4. Confidence score   → method_agreement + confidence_label

Feature loading strategy:
  - Reads features directly from carbodb.sqlite (no re-loading TSV files)
  - Processes in chunks of 10,000 sequences
  - Uses same feature vector as training (1,793 features)
  - EC one-hot + kingdom added for Km prediction

Output:
  carbodb.sqlite — predictions table filled (one row per sequence per model)
  carbodb.sqlite — confidence_scores table filled

Usage:
  python scripts/10_predict_all.py
  python scripts/10_predict_all.py --dry-run      (first 1000 seqs only)
  python scripts/10_predict_all.py --tasks binary ec km
  python scripts/10_predict_all.py --update-km    (re-run Km only)
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
import xgboost as xgb
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, setup_logging

log = setup_logging("10_predict_all")

DB_PATH   = PATHS.PRIMARY / "carbodb.sqlite"
ML_DIR    = PATHS.PRIMARY.parent / "ml"
MODEL_DIR = PATHS.PRIMARY.parent / "models"

MODEL_VERSION  = "v5"
CHUNK_SIZE     = 10_000
ESM2_DIM       = 1280

# EC one-hot columns (from training)
EC_OH_COLS = [
    "ec_oh_1.1.1.39","ec_oh_1.2.7.7","ec_oh_4.1.1.31","ec_oh_4.1.1.32",
    "ec_oh_4.1.1.39","ec_oh_4.1.1.49","ec_oh_4.1.1.90","ec_oh_4.2.1.1",
    "ec_oh_6.3.3.3","ec_oh_6.3.4.14","ec_oh_6.3.4.16","ec_oh_6.3.4.18",
    "ec_oh_6.3.5.5","ec_oh_6.4.1.1","ec_oh_6.4.1.2","ec_oh_6.4.1.3",
    "ec_oh_6.4.1.4",
]
KINGDOM_COLS = [
    "kingdom_archaea","kingdom_bacteria","kingdom_fungi","kingdom_plant"
]

TRAINABLE_KM_EC = {
    "4.2.1.1","4.1.1.39","4.1.1.31","4.1.1.49",
    "6.3.4.14","4.1.1.32","6.4.1.1","6.4.1.4","6.4.1.2","6.4.1.3"
}


# ── DB helpers ────────────────────────────────────────────────────────────

def get_conn():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA cache_size = -2000000")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("PRAGMA journal_mode = WAL")
    return conn


def table_count(conn, table):
    try:
        return conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    except Exception:
        return 0


# ── Feature loading ───────────────────────────────────────────────────────

COMP_COLS = (
    [f"aac_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"] +
    ["phys_length","phys_length_log","phys_mw","phys_pi","phys_gravy",
     "phys_aromaticity","phys_instability","phys_charge_ph7","phys_frac_charged",
     "phys_frac_aromatic","phys_frac_polar","phys_frac_nonpolar","phys_frac_small",
     "phys_frac_glycine","phys_frac_proline"] +
    ["inv_cat_D","inv_cat_E","inv_cat_H","inv_cat_K","inv_cat_C","inv_cat_S",
     "inv_cat_T","inv_cat_mean_dist","inv_cat_std_dist","inv_cat_min_dist",
     "inv_cat_max_dist","inv_cat_clustering","inv_hydrophobic","inv_charged",
     "inv_polar","inv_aromatic","inv_net_charge"] +
    ["motif_rubisco_kk","motif_rubisco_gk","motif_ca_hh","motif_ca_his_cluster",
     "motif_pepc_rr","motif_biotin_mk","motif_biotin_amk"] +
    [f"dp_{a1}{a2}" for a1 in "ACDEFGHIKLMNPQRSTVWY"
                    for a2 in "ACDEFGHIKLMNPQRSTVWY"] +
    [f"pse_{a}" for a in "ACDEFGHIKLMNPQRSTVWY"] +
    [f"pse_corr_{i}" for i in range(1, 11)]
)

PFAM_COLS = [
    "pfam_PF00016","pfam_PF02788","pfam_PF00101","pfam_PF00194","pfam_PF03119",
    "pfam_PF00311","pfam_PF00821","pfam_PF02785","pfam_PF00364","pfam_PF01039",
    "pfam_PF02786","pfam_PF02787","pfam_PF00289","pfam_PF01309","pfam_PF03599",
    "pfam_PF03590","pfam_PF00384","pfam_PF00682","pfam_n_hits",
]

IPR_COLS = [
    "n_panther","n_gene3d","n_tigrfam","n_prosite_prof","n_prosite_pat",
]


def infer_kingdom(organism):
    if not isinstance(organism, str):
        return "bacteria"
    o = organism.lower()
    if any(x in o for x in ["arabidopsis","oryza","zea","solanum","nicotiana",
                              "spinacia","triticum","glycine","brassica",
                              "chlamydomonas","chlorella","cyanophora"]):
        return "plant"
    if any(x in o for x in ["methan","sulfolobus","pyrococcus","halobacter",
                              "archaeo","thermococcus","archaeoglobus"]):
        return "archaea"
    if any(x in o for x in ["saccharomyces","aspergillus","candida"]):
        return "fungi"
    return "bacteria"


def load_features_chunk(conn, sequence_ids):
    """Load all features for a chunk of sequence IDs from DB."""
    conn.execute("DROP TABLE IF EXISTS _pred_ids")
    conn.execute("CREATE TEMP TABLE _pred_ids (id INTEGER)")
    conn.executemany("INSERT INTO _pred_ids VALUES (?)",
                     [(int(i),) for i in sequence_ids])

    # Sequences (for EC, kingdom, organism)
    seqs_df = pd.read_sql_query("""
        SELECT s.id, s.cdb_id, s.uniprot_id, s.ec_number, s.organism
        FROM sequences s
        JOIN _pred_ids p ON p.id = s.id
    """, conn)

    # Composition
    comp_cols_str = ", ".join([f"c.{col}" for col in COMP_COLS])
    comp_df = pd.read_sql_query(f"""
        SELECT s.id, {comp_cols_str}
        FROM features_composition c
        JOIN sequences s ON s.id = c.sequence_id
        JOIN _pred_ids p ON p.id = s.id
    """, conn)

    # Domains
    pfam_cols_str = ", ".join([f"d.{col}" for col in PFAM_COLS])
    dom_df = pd.read_sql_query(f"""
        SELECT s.id, {pfam_cols_str}
        FROM features_domains d
        JOIN sequences s ON s.id = d.sequence_id
        JOIN _pred_ids p ON p.id = s.id
    """, conn)

    # InterPro
    ipr_cols_str = ", ".join([f"i.{col}" for col in IPR_COLS])
    ipr_df = pd.read_sql_query(f"""
        SELECT s.id, {ipr_cols_str}
        FROM features_interpro i
        JOIN sequences s ON s.id = i.sequence_id
        JOIN _pred_ids p ON p.id = s.id
    """, conn)

    # ESM-2 blobs
    rows = conn.execute("""
        SELECT s.id, e.embedding_blob
        FROM features_esm2 e
        JOIN sequences s ON s.id = e.sequence_id
        JOIN _pred_ids p ON p.id = s.id
    """).fetchall()
    esm2_data = {}
    for sid, blob in rows:
        if blob:
            emb = np.frombuffer(blob, dtype=np.float32)
            if len(emb) == ESM2_DIM:
                esm2_data[sid] = emb

    # Merge all features
    merged = seqs_df.merge(comp_df,  on="id", how="left")
    merged = merged.merge(dom_df,   on="id", how="left")
    merged = merged.merge(ipr_df,   on="id", how="left")

    # Add ESM-2
    esm2_mat = np.zeros((len(merged), ESM2_DIM), dtype=np.float32)
    for i, sid in enumerate(merged["id"].values):
        if sid in esm2_data:
            esm2_mat[i] = esm2_data[sid]

    esm2_df = pd.DataFrame(
        esm2_mat,
        columns=[f"esm2_{j}" for j in range(ESM2_DIM)],
        index=merged.index
    )
    merged = pd.concat([merged, esm2_df], axis=1)

    return merged


def build_X_base(merged_df):
    """Build base feature matrix (1793 features, no EC/kingdom)."""
    feat_cols = COMP_COLS + PFAM_COLS + IPR_COLS + \
                [f"esm2_{j}" for j in range(ESM2_DIM)]
    X = merged_df[feat_cols].fillna(0.0).values.astype(np.float32)
    return X


def build_X_km(merged_df):
    """Build Km feature matrix (1793 + 17 EC one-hot + 4 kingdom = 1814)."""
    X_base = build_X_base(merged_df)

    # EC one-hot
    ec_oh = np.zeros((len(merged_df), len(EC_OH_COLS)), dtype=np.float32)
    ec_col_map = {col.replace("ec_oh_", ""): i for i, col in enumerate(EC_OH_COLS)}
    for i, ec in enumerate(merged_df["ec_number"].values):
        if ec in ec_col_map:
            ec_oh[i, ec_col_map[ec]] = 1.0

    # Kingdom one-hot
    kingdom_oh = np.zeros((len(merged_df), len(KINGDOM_COLS)), dtype=np.float32)
    kingdom_col_map = {col.replace("kingdom_", ""): i
                       for i, col in enumerate(KINGDOM_COLS)}
    for i, org in enumerate(merged_df["organism"].values):
        k = infer_kingdom(org)
        if k in kingdom_col_map:
            kingdom_oh[i, kingdom_col_map[k]] = 1.0

    return np.hstack([X_base, ec_oh, kingdom_oh])


# ── Prediction helpers ────────────────────────────────────────────────────

def predict_binary(booster, X):
    return booster.predict(xgb.DMatrix(X))  # probabilities


def predict_ec(booster, X, ec_map):
    probs = booster.predict(xgb.DMatrix(X)).reshape(len(X), -1)
    inv_map = {v: k for k, v in ec_map.items()}
    pred_int = probs.argmax(axis=1)
    pred_ec  = [inv_map.get(int(i), "unknown") for i in pred_int]
    pred_prob = probs.max(axis=1)
    return pred_ec, pred_prob


def predict_km(booster, X, ec_series):
    """Only predict Km for sequences in trainable EC classes."""
    mask = ec_series.isin(TRAINABLE_KM_EC).values
    km_pred = np.full(len(X), np.nan, dtype=np.float32)
    if mask.sum() > 0:
        log10_pred = booster.predict(xgb.DMatrix(X[mask]))
        km_pred[mask] = 10 ** log10_pred  # convert log10 → mM
    return km_pred


def confidence_label(co2_prob, ec_prob):
    """Assign confidence label based on model agreement."""
    if co2_prob >= 0.95 and ec_prob >= 0.90:
        return "high"
    elif co2_prob >= 0.80 and ec_prob >= 0.70:
        return "medium"
    elif co2_prob >= 0.50:
        return "low"
    else:
        return "review"


# ── Main prediction loop ──────────────────────────────────────────────────

def predict_all(tasks, dry_run=False):
    log.info("Loading models...")
    binary_booster = xgb.Booster()
    binary_booster.load_model(MODEL_DIR / "binary_v5.json")
    log.info("  binary_v5.json loaded")

    ec_booster = xgb.Booster()
    ec_booster.load_model(MODEL_DIR / "ec_v5.json")
    ec_map = json.load(open(ML_DIR / "ec_label_map_fixed.json"))
    log.info("  ec_v5.json loaded (%d classes)", len(ec_map))

    km_booster = xgb.Booster()
    km_booster.load_model(MODEL_DIR / "km_v5_weighted.json")
    log.info("  km_v5_weighted.json loaded")

    conn = get_conn()
    total_seqs = table_count(conn, "sequences")
    log.info("Total sequences: %d", total_seqs)

    if dry_run:
        total_seqs = min(1000, total_seqs)
        log.info("DRY RUN — processing first %d sequences", total_seqs)

    # Get all sequence IDs
    all_ids = [r[0] for r in conn.execute(
        f"SELECT id FROM sequences ORDER BY id LIMIT {total_seqs}"
    ).fetchall()]

    # Clear existing predictions for this model version
    if not dry_run:
        conn.execute("DELETE FROM predictions WHERE model_version = ?",
                     (MODEL_VERSION,))
        conn.execute("DELETE FROM confidence_scores")
        conn.commit()
        log.info("Cleared existing v5 predictions")

    pred_rows   = []
    conf_rows   = []
    n_processed = 0
    t0 = time.time()

    for chunk_start in tqdm(range(0, len(all_ids), CHUNK_SIZE),
                             desc="predict chunks"):
        chunk_ids = all_ids[chunk_start:chunk_start + CHUNK_SIZE]

        try:
            merged = load_features_chunk(conn, chunk_ids)
        except Exception as e:
            log.warning("  Chunk %d failed: %s", chunk_start, e)
            continue

        if len(merged) == 0:
            continue

        X_base = build_X_base(merged)
        X_km   = build_X_km(merged)

        # Binary predictions
        co2_probs = predict_binary(binary_booster, X_base)
        is_co2    = (co2_probs >= 0.5).astype(int)

        # EC predictions
        ec_preds, ec_probs = predict_ec(ec_booster, X_base, ec_map)

        # Km predictions
        km_preds = predict_km(km_booster, X_km, merged["ec_number"])

        # Build rows
        for i in range(len(merged)):
            row = merged.iloc[i]
            seq_id    = int(row["id"])
            uid       = row["uniprot_id"]
            km_val    = None if np.isnan(km_preds[i]) else float(km_preds[i])
            km_log    = None if km_val is None else float(np.log10(km_val))

            pred_rows.append({
                "sequence_id":   seq_id,
                "uniprot_id":    uid,
                "model_version": MODEL_VERSION,
                "is_co2_pred":   int(is_co2[i]),
                "co2_prob":      float(co2_probs[i]),
                "ec_pred":       ec_preds[i],
                "ec_prob":       float(ec_probs[i]),
                "km_pred_mM":    km_val,
                "km_pred_log10": km_log,
            })

            conf_rows.append({
                "sequence_id":      seq_id,
                "uniprot_id":       uid,
                "method_agreement": 1,  # single model, updated in v6
                "ec_confidence":    float(ec_probs[i]),
                "km_confidence":    float(ec_probs[i]) if km_val else None,
                "overall_score":    float((co2_probs[i] + ec_probs[i]) / 2),
                "confidence_label": confidence_label(
                    float(co2_probs[i]), float(ec_probs[i])),
            })

        n_processed += len(merged)

        # Write to DB every 5 chunks
        if len(pred_rows) >= CHUNK_SIZE * 5 and not dry_run:
            _flush_predictions(conn, pred_rows, conf_rows)
            elapsed = time.time() - t0
            rate    = n_processed / elapsed
            eta_min = (total_seqs - n_processed) / rate / 60
            log.info("  %d / %d  (%.0f seq/s  ETA %.0f min)",
                     n_processed, total_seqs, rate, eta_min)
            pred_rows = []
            conf_rows = []

    # Final flush
    if pred_rows and not dry_run:
        _flush_predictions(conn, pred_rows, conf_rows)

    # Update metadata
    if not dry_run:
        conn.execute("INSERT OR REPLACE INTO db_metadata VALUES (?,?)",
                     ("predictions_v5_count", str(table_count(conn, "predictions"))))
        conn.execute("INSERT OR REPLACE INTO db_metadata VALUES (?,?)",
                     ("predictions_v5_date", TS))
        conn.commit()

    elapsed = time.time() - t0
    log.info("Processed %d sequences in %.1f min", n_processed, elapsed / 60)
    log.info("Predictions table: %d rows", table_count(conn, "predictions"))
    log.info("Confidence table:  %d rows", table_count(conn, "confidence_scores"))
    conn.close()


def _flush_predictions(conn, pred_rows, conf_rows):
    conn.executemany("""
        INSERT OR REPLACE INTO predictions
        (sequence_id, uniprot_id, model_version, is_co2_pred, co2_prob,
         ec_pred, ec_prob, km_pred_mM, km_pred_log10)
        VALUES
        (:sequence_id, :uniprot_id, :model_version, :is_co2_pred, :co2_prob,
         :ec_pred, :ec_prob, :km_pred_mM, :km_pred_log10)
    """, pred_rows)

    conn.executemany("""
        INSERT OR REPLACE INTO confidence_scores
        (sequence_id, uniprot_id, method_agreement, ec_confidence,
         km_confidence, overall_score, confidence_label)
        VALUES
        (:sequence_id, :uniprot_id, :method_agreement, :ec_confidence,
         :km_confidence, :overall_score, :confidence_label)
    """, conf_rows)

    conn.commit()


# ── Summary stats after prediction ───────────────────────────────────────

def print_summary(conn):
    log.info("\n" + "=" * 60)
    log.info("PREDICTION SUMMARY")
    log.info("=" * 60)

    total = table_count(conn, "predictions")
    log.info("Total predictions: %d", total)

    # Predicted carboxylases
    n_pos = conn.execute(
        "SELECT COUNT(*) FROM predictions WHERE is_co2_pred=1 AND model_version=?",
        (MODEL_VERSION,)).fetchone()[0]
    log.info("Predicted carboxylases: %d (%.1f%%)", n_pos, 100*n_pos/total)

    # EC distribution of predicted carboxylases
    ec_dist = conn.execute("""
        SELECT ec_pred, COUNT(*) as n
        FROM predictions
        WHERE is_co2_pred=1 AND model_version=?
        GROUP BY ec_pred ORDER BY n DESC LIMIT 10
    """, (MODEL_VERSION,)).fetchall()
    log.info("Top predicted EC classes:")
    for ec, n in ec_dist:
        log.info("  %-15s  %d", ec, n)

    # Confidence distribution
    conf_dist = conn.execute("""
        SELECT confidence_label, COUNT(*) as n
        FROM confidence_scores GROUP BY confidence_label
    """).fetchall()
    log.info("Confidence distribution:")
    for label, n in conf_dist:
        log.info("  %-10s  %d", label, n)

    # Km predictions
    n_km = conn.execute(
        "SELECT COUNT(*) FROM predictions WHERE km_pred_mM IS NOT NULL AND model_version=?",
        (MODEL_VERSION,)).fetchone()[0]
    log.info("Sequences with Km prediction: %d", n_km)

    if n_km > 0:
        km_stats = conn.execute("""
            SELECT AVG(km_pred_log10), MIN(km_pred_log10), MAX(km_pred_log10)
            FROM predictions WHERE km_pred_mM IS NOT NULL AND model_version=?
        """, (MODEL_VERSION,)).fetchone()
        log.info("Predicted Km range: 10^%.2f to 10^%.2f mM (mean=10^%.2f)",
                 km_stats[1], km_stats[2], km_stats[0])


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run",   action="store_true",
                    help="Process first 1000 sequences only")
    ap.add_argument("--tasks",     nargs="+", default=["binary", "ec", "km"])
    ap.add_argument("--update-km", action="store_true",
                    help="Re-run Km predictions only")
    args = ap.parse_args()

    if not DB_PATH.exists():
        log.error("Database not found: %s", DB_PATH)
        sys.exit(1)

    log.info("CarboDB v5 — predict all sequences")
    log.info("Database: %s", DB_PATH)
    log.info("Dry run: %s", args.dry_run)

    predict_all(args.tasks, dry_run=args.dry_run)

    if not args.dry_run:
        conn = get_conn()
        print_summary(conn)
        conn.close()

    log.info("Done. Database updated with v5 predictions.")
    log.info("Next: CarboDB-App webapp development")


if __name__ == "__main__":
    main()
