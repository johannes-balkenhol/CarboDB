#!/usr/bin/env python3
"""
14_biological_analysis.py
=========================
CarboDB — Step 14: Biological top-hits analysis.

Tasks:
  A. Validation controls — known carboxylases predicted correctly
  B. Novel high-Km candidates — extreme Km predictions not in BRENDA training
  C. Novel low-Km candidates — ultra-high-affinity predictions
  D. Taxonomic enrichment — which kingdoms/phyla have most carboxylases
  E. EC class distribution across organisms
  F. Top predicted carboxylases by confidence (novel discoveries)

Output: data/biological/
  validation_controls.json      known sequences + predictions
  novel_high_km.tsv             top 500 novel high-Km candidates
  novel_low_km.tsv              top 500 novel low-Km candidates  
  taxonomic_summary.json        kingdom/organism breakdowns
  ec_distribution.json          EC class distribution across taxa
  top_novel_carboxylases.tsv    high-confidence novel predictions

Usage:
  python scripts/14_biological_analysis.py
  python scripts/14_biological_analysis.py --tasks A B C
  python scripts/14_biological_analysis.py --top-n 200
"""

import argparse
import json
import sqlite3
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        if isinstance(obj, (np.bool_,)): return bool(obj)
        return super().default(obj)
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, ROOT, TS, setup_logging

log = setup_logging("14_biological_analysis")

PRIMARY   = ROOT / "data" / "primary"
DB_PATH   = PRIMARY / "carbodb.sqlite"
BIO_DIR   = ROOT / "data" / "biological"
FIG_DIR   = BIO_DIR / "figures"

BIO_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Known control sequences for validation
KNOWN_CONTROLS = [
    # (uniprot_id, common_name, ec_expected, km_expected_mM, notes)
    ("P00875", "Spinach RuBisCO (Form I)",          "4.1.1.39", 0.01,  "C3 plant, high-affinity Form I"),
    ("P00880", "R. rubrum RuBisCO (Form II)",        "4.1.1.39", 1.0,   "Photosynthetic bacterium, low-affinity Form II"),
    ("P00918", "Human carbonic anhydrase II",         "4.2.1.1",  8.0,   "Fast CA isoform, high Km"),
    ("P19819", "Human carbonic anhydrase I",          "4.2.1.1",  4.7,   "Slow CA isoform"),
    ("P00864", "E. coli PEP carboxylase",             "4.1.1.31", 0.9,   "C4 acid cycle, bacterial"),
    ("P11498", "Human pyruvate carboxylase",          "6.4.1.1",  0.4,   "Mitochondrial, biotin-dependent"),
    ("P05165", "Human propionyl-CoA carboxylase",     "6.4.1.3",  0.5,   "Biotin-dependent, rare EC"),
    ("P0ABD5", "E. coli acetyl-CoA carboxylase",     "6.4.1.2",  0.2,   "Fatty acid synthesis, ACC"),
    ("P15977", "Maize PEPC (C4)",                    "4.1.1.31", 0.07,  "C4 plant, high-affinity PEPC"),
    ("Q8TYR5", "Archaeal RuBisCO (Form III)",         "4.1.1.39", 20.0,  "Thermophile, very high Km"),
]

EC_NAMES = {
    "4.1.1.39": "RuBisCO",
    "4.2.1.1":  "Carbonic anhydrase",
    "6.3.4.16": "ACC biotin carboxylase",
    "6.3.4.14": "Pyruvate carboxylase",
    "6.3.5.5":  "Carbamoyl-P synthase",
    "6.3.4.18": "3-MCC",
    "4.1.1.49": "PEPC",
    "6.3.3.3":  "Dethiobiotin synthase",
    "4.1.1.31": "PEPCK",
    "4.1.1.112":"2-OG carboxylase",
    "4.1.1.32": "PEPCK-GTP",
    "6.4.1.1":  "Pyruvate carboxylase",
    "6.4.1.2":  "ACC",
    "6.4.1.3":  "Propionyl-CoA carboxylase",
    "6.4.1.4":  "3-MCC",
}


# ══════════════════════════════════════════════════════════════════════════════
# DB connection
# ══════════════════════════════════════════════════════════════════════════════

def get_conn():
    if not DB_PATH.exists():
        log.error("Database not found: %s", DB_PATH)
        sys.exit(1)
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def df_query(conn, sql, params=()):
    return pd.read_sql_query(sql, conn, params=params)


# ══════════════════════════════════════════════════════════════════════════════
# TASK A: Validation controls
# ══════════════════════════════════════════════════════════════════════════════

def task_a_validation_controls(conn):
    log.info("══ Task A: Validation Controls ══")

    results = []
    for uniprot_id, name, ec_expected, km_expected, notes in KNOWN_CONTROLS:
        row = df_query(conn, """
            SELECT
                s.cdb_id, s.uniprot_id, s.organism, s.length,
                s.ec_number AS ec_experimental,
                s.km_best_mM AS km_experimental,
                s.reviewed, s.source,
                p.is_co2_pred, p.co2_prob,
                p.ec_pred, p.ec_prob,
                p.km_pred_mM, p.km_pred_log10,
                c.confidence_label
            FROM sequences s
            JOIN predictions p ON p.sequence_id = s.id
            JOIN confidence_scores c ON c.sequence_id = s.id
            WHERE s.uniprot_id = ?
        """, (uniprot_id,))

        if len(row) == 0:
            log.warning("  %s (%s): NOT IN DATABASE", uniprot_id, name)
            results.append({
                "uniprot_id":     uniprot_id,
                "name":           name,
                "ec_expected":    ec_expected,
                "km_expected_mM": km_expected,
                "notes":          notes,
                "in_database":    False,
                "status":         "MISSING",
            })
            continue

        r = dict(row.iloc[0])
        ec_correct  = r["ec_pred"] == ec_expected
        is_carb     = bool(r["is_co2_pred"])
        km_pred     = r["km_pred_mM"]
        km_fold_err = abs(np.log10(km_pred) - np.log10(km_expected)) if km_pred and km_expected else None

        status_parts = []
        if is_carb:             status_parts.append("carboxylase ✓")
        else:                   status_parts.append("carboxylase ✗")
        if ec_correct:          status_parts.append(f"EC ✓ ({ec_expected})")
        else:                   status_parts.append(f"EC ✗ (pred={r['ec_pred']} expected={ec_expected})")
        if km_fold_err is not None:
            if km_fold_err < 0.5:   status_parts.append(f"Km ✓ ({km_pred:.4f} vs {km_expected} mM, {km_fold_err:.2f} log10)")
            elif km_fold_err < 1.0: status_parts.append(f"Km ~ ({km_pred:.4f} vs {km_expected} mM, {km_fold_err:.2f} log10)")
            else:                   status_parts.append(f"Km ✗ ({km_pred:.4f} vs {km_expected} mM, {km_fold_err:.2f} log10)")

        log.info("  %s (%s): %s", uniprot_id, name, " | ".join(status_parts))

        results.append({
            "uniprot_id":          uniprot_id,
            "name":                name,
            "organism":            r["organism"],
            "length":              r["length"],
            "notes":               notes,
            "in_database":         True,
            "source":              r["source"],
            "reviewed":            bool(r["reviewed"]),

            "ec_expected":         ec_expected,
            "ec_predicted":        r["ec_pred"],
            "ec_prob":             round(float(r["ec_prob"]), 4) if r["ec_prob"] else None,
            "ec_correct":          ec_correct,

            "is_carboxylase_pred": is_carb,
            "co2_prob":            round(float(r["co2_prob"]), 4),
            "confidence":          r["confidence_label"],

            "km_expected_mM":      km_expected,
            "km_predicted_mM":     round(float(km_pred), 4) if km_pred else None,
            "km_fold_error_log10": round(km_fold_err, 3) if km_fold_err is not None else None,
            "km_within_2fold":     km_fold_err is not None and km_fold_err < np.log10(2),
            "km_within_5fold":     km_fold_err is not None and km_fold_err < np.log10(5),

            "status":              " | ".join(status_parts),
        })

    # Summary
    n_in_db   = sum(1 for r in results if r["in_database"])
    n_carb_ok = sum(1 for r in results if r.get("is_carboxylase_pred"))
    n_ec_ok   = sum(1 for r in results if r.get("ec_correct"))
    n_km_2x   = sum(1 for r in results if r.get("km_within_2fold"))
    n_km_5x   = sum(1 for r in results if r.get("km_within_5fold"))

    log.info("  Summary: %d/%d in DB | %d/%d carboxylase ✓ | %d/%d EC ✓ | %d/%d Km within 2-fold | %d/%d Km within 5-fold",
             n_in_db, len(results), n_carb_ok, n_in_db, n_ec_ok, n_in_db,
             n_km_2x, n_in_db, n_km_5x, n_in_db)

    result = {
        "task": "A_validation_controls",
        "n_controls": len(results),
        "n_in_database": n_in_db,
        "n_carboxylase_correct": n_carb_ok,
        "n_ec_correct": n_ec_ok,
        "n_km_within_2fold": n_km_2x,
        "n_km_within_5fold": n_km_5x,
        "controls": results,
    }

    json.dump(result, open(BIO_DIR / "validation_controls.json", "w"), indent=2, cls=NumpyEncoder)
    pd.DataFrame(results).to_csv(FIG_DIR / "validation_controls.tsv", sep="\t", index=False)
    log.info("  Saved: validation_controls.json")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK B+C: Novel Km candidates
# ══════════════════════════════════════════════════════════════════════════════

def task_bc_novel_km_candidates(conn, top_n=500):
    log.info("══ Tasks B+C: Novel Km Candidates ══")

    # Novel = not in BRENDA (km_experimental IS NULL) + high confidence + has Km prediction
    base_sql = """
        SELECT
            s.cdb_id, s.uniprot_id, s.organism, s.length,
            s.ec_number AS ec_experimental,
            s.km_best_mM AS km_experimental,
            s.reviewed, s.source,
            p.ec_pred, p.ec_prob,
            p.km_pred_mM, p.km_pred_log10,
            p.co2_prob,
            c.confidence_label
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        JOIN confidence_scores c ON c.sequence_id = s.id
        WHERE p.is_co2_pred = 1
          AND p.km_pred_mM IS NOT NULL
          AND s.km_best_mM IS NULL
          AND c.confidence_label IN ('high', 'medium')
    """

    # Task B: High Km (potential extremophiles / low-affinity variants)
    log.info("  Task B: Novel high-Km candidates")
    high_km = df_query(conn, base_sql + " ORDER BY p.km_pred_mM DESC LIMIT ?", (top_n,))
    high_km["rank"] = range(1, len(high_km) + 1)
    high_km["ec_name"] = high_km["ec_pred"].map(EC_NAMES).fillna("")
    high_km.to_csv(FIG_DIR / "novel_high_km.tsv", sep="\t", index=False)
    log.info("  Top high-Km: %.3f mM (%s, %s)",
             high_km.iloc[0]["km_pred_mM"],
             high_km.iloc[0]["uniprot_id"],
             high_km.iloc[0]["organism"])
    log.info("  Top-5 high-Km:")
    for _, row in high_km.head(5).iterrows():
        log.info("    %s | %s | %s | %.3f mM | conf=%s",
                 row["uniprot_id"], row["organism"][:40],
                 row["ec_pred"], row["km_pred_mM"], row["confidence_label"])

    # Task C: Low Km (ultra-high-affinity — relevant for carbon capture)
    log.info("  Task C: Novel low-Km candidates")
    low_km = df_query(conn, base_sql + " ORDER BY p.km_pred_mM ASC LIMIT ?", (top_n,))
    low_km["rank"] = range(1, len(low_km) + 1)
    low_km["ec_name"] = low_km["ec_pred"].map(EC_NAMES).fillna("")
    low_km.to_csv(FIG_DIR / "novel_low_km.tsv", sep="\t", index=False)
    log.info("  Top low-Km: %.6f mM (%s, %s)",
             low_km.iloc[0]["km_pred_mM"],
             low_km.iloc[0]["uniprot_id"],
             low_km.iloc[0]["organism"])
    log.info("  Top-5 low-Km:")
    for _, row in low_km.head(5).iterrows():
        log.info("    %s | %s | %s | %.6f mM | conf=%s",
                 row["uniprot_id"], row["organism"][:40],
                 row["ec_pred"], row["km_pred_mM"], row["confidence_label"])

    # EC class breakdown of top candidates
    log.info("  EC class breakdown of top-%d high-Km:", top_n)
    for ec, grp in high_km.groupby("ec_pred"):
        log.info("    %s (%s): n=%d, mean_km=%.3f mM, range=%.4f-%.3f mM",
                 ec, EC_NAMES.get(ec, ""), len(grp),
                 grp["km_pred_mM"].mean(),
                 grp["km_pred_mM"].min(),
                 grp["km_pred_mM"].max())

    result = {
        "task": "BC_novel_km_candidates",
        "n_high_km": len(high_km),
        "n_low_km":  len(low_km),
        "high_km_max_mM":  round(float(high_km["km_pred_mM"].max()), 3),
        "high_km_top5":    high_km.head(5)[["uniprot_id","organism","ec_pred","km_pred_mM","confidence_label"]].to_dict("records"),
        "low_km_min_mM":   round(float(low_km["km_pred_mM"].min()), 6),
        "low_km_top5":     low_km.head(5)[["uniprot_id","organism","ec_pred","km_pred_mM","confidence_label"]].to_dict("records"),
    }

    json.dump(result, open(BIO_DIR / "novel_km_candidates.json", "w"), indent=2, cls=NumpyEncoder)
    log.info("  Saved: novel_high_km.tsv, novel_low_km.tsv, novel_km_candidates.json")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK D: Taxonomic enrichment
# ══════════════════════════════════════════════════════════════════════════════

def task_d_taxonomic_analysis(conn):
    log.info("══ Task D: Taxonomic Analysis ══")

    # Kingdom-level breakdown
    log.info("  Kingdom breakdown of predicted carboxylases...")
    kingdom_sql = """
        SELECT
            CASE
                WHEN s.organism LIKE '%virus%' OR s.organism LIKE '%phage%' THEN 'Viruses'
                WHEN s.organism LIKE '%archaea%' OR s.organism LIKE '%archaeon%'
                     OR s.organism LIKE '%thermophil%' AND s.source != 'swissprot' THEN 'Archaea'
                WHEN s.organism LIKE '% sp.%' OR s.organism LIKE '%uncultured%'
                     OR s.organism LIKE '%metagenom%' THEN 'Unclassified/Metagenome'
                WHEN s.reviewed = 1 THEN 'SwissProt (reviewed)'
                ELSE 'TrEMBL (unreviewed)'
            END AS kingdom_group,
            COUNT(*) AS n_total,
            SUM(p.is_co2_pred) AS n_carboxylase,
            ROUND(AVG(CASE WHEN p.is_co2_pred=1 THEN p.km_pred_mM END), 4) AS mean_km_pred,
            ROUND(AVG(CASE WHEN p.is_co2_pred=1 THEN p.co2_prob END), 4) AS mean_co2_prob
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        GROUP BY kingdom_group
        ORDER BY n_carboxylase DESC
    """
    kingdom_df = df_query(conn, kingdom_sql)
    kingdom_df["pct_carboxylase"] = (kingdom_df["n_carboxylase"] / kingdom_df["n_total"] * 100).round(2)
    log.info("  Kingdom breakdown:")
    for _, row in kingdom_df.iterrows():
        log.info("    %-35s n_carb=%6d  %%=%.1f%%  mean_Km=%.4f mM",
                 row["kingdom_group"], row["n_carboxylase"],
                 row["pct_carboxylase"],
                 row["mean_km_pred"] if row["mean_km_pred"] else 0)

    # Top organisms by carboxylase count
    log.info("  Top organisms by predicted carboxylase count...")
    top_organisms_sql = """
        SELECT
            s.organism,
            COUNT(*) AS n_carboxylase,
            GROUP_CONCAT(DISTINCT p.ec_pred) AS ec_classes,
            ROUND(AVG(p.km_pred_mM), 4) AS mean_km_pred,
            ROUND(AVG(p.co2_prob), 4) AS mean_confidence,
            SUM(CASE WHEN s.reviewed=1 THEN 1 ELSE 0 END) AS n_reviewed
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        WHERE p.is_co2_pred = 1
          AND s.organism IS NOT NULL
          AND s.organism != ''
        GROUP BY s.organism
        HAVING n_carboxylase >= 5
        ORDER BY n_carboxylase DESC
        LIMIT 50
    """
    top_orgs = df_query(conn, top_organisms_sql)
    log.info("  Top-10 organisms by carboxylase count:")
    for _, row in top_orgs.head(10).iterrows():
        log.info("    %-45s n=%5d  mean_Km=%.4f mM  ECs=%s",
                 row["organism"][:45], row["n_carboxylase"],
                 row["mean_km_pred"] if row["mean_km_pred"] else 0,
                 row["ec_classes"][:50] if row["ec_classes"] else "")

    # Top organisms by mean predicted Km (high-Km organisms)
    log.info("  Top organisms by mean predicted Km (high-affinity = low Km)...")
    high_km_orgs_sql = """
        SELECT
            s.organism,
            COUNT(*) AS n_carboxylase,
            ROUND(AVG(p.km_pred_mM), 4) AS mean_km_pred,
            ROUND(MIN(p.km_pred_mM), 6) AS min_km_pred,
            ROUND(MAX(p.km_pred_mM), 3) AS max_km_pred,
            GROUP_CONCAT(DISTINCT p.ec_pred) AS ec_classes
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        WHERE p.is_co2_pred = 1
          AND p.km_pred_mM IS NOT NULL
          AND s.organism IS NOT NULL
        GROUP BY s.organism
        HAVING n_carboxylase >= 3
        ORDER BY mean_km_pred DESC
        LIMIT 30
    """
    high_km_orgs = df_query(conn, high_km_orgs_sql)
    log.info("  Top-10 organisms by highest mean Km (potential extremophiles):")
    for _, row in high_km_orgs.head(10).iterrows():
        log.info("    %-45s mean_Km=%.3f mM  n=%d  ECs=%s",
                 row["organism"][:45], row["mean_km_pred"],
                 row["n_carboxylase"],
                 row["ec_classes"][:40] if row["ec_classes"] else "")

    # Low-Km organisms (high affinity — carbon capture candidates)
    low_km_orgs_sql = """
        SELECT
            s.organism,
            COUNT(*) AS n_carboxylase,
            ROUND(AVG(p.km_pred_mM), 6) AS mean_km_pred,
            GROUP_CONCAT(DISTINCT p.ec_pred) AS ec_classes
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        WHERE p.is_co2_pred = 1
          AND p.km_pred_mM IS NOT NULL
          AND s.organism IS NOT NULL
        GROUP BY s.organism
        HAVING n_carboxylase >= 3
        ORDER BY mean_km_pred ASC
        LIMIT 30
    """
    low_km_orgs = df_query(conn, low_km_orgs_sql)
    log.info("  Top-10 organisms by lowest mean Km (high-affinity carboxylases):")
    for _, row in low_km_orgs.head(10).iterrows():
        log.info("    %-45s mean_Km=%.6f mM  n=%d",
                 row["organism"][:45], row["mean_km_pred"], row["n_carboxylase"])

    # Save outputs
    kingdom_df.to_csv(FIG_DIR / "kingdom_breakdown.tsv", sep="\t", index=False)
    top_orgs.to_csv(FIG_DIR / "top_organisms_by_count.tsv", sep="\t", index=False)
    high_km_orgs.to_csv(FIG_DIR / "top_organisms_high_km.tsv", sep="\t", index=False)
    low_km_orgs.to_csv(FIG_DIR / "top_organisms_low_km.tsv", sep="\t", index=False)

    result = {
        "task": "D_taxonomic_analysis",
        "kingdom_breakdown": kingdom_df.to_dict("records"),
        "top_organisms_by_count": top_orgs.head(20).to_dict("records"),
        "top_organisms_high_km": high_km_orgs.head(20).to_dict("records"),
        "top_organisms_low_km": low_km_orgs.head(20).to_dict("records"),
    }
    json.dump(result, open(BIO_DIR / "taxonomic_summary.json", "w"), indent=2, cls=NumpyEncoder)
    log.info("  Saved: taxonomic_summary.json + 4 TSV files")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK E: EC class distribution
# ══════════════════════════════════════════════════════════════════════════════

def task_e_ec_distribution(conn):
    log.info("══ Task E: EC Class Distribution ══")

    # Overall EC distribution
    ec_dist_sql = """
        SELECT
            p.ec_pred,
            COUNT(*) AS n_predicted,
            ROUND(AVG(p.ec_prob), 4) AS mean_ec_prob,
            ROUND(AVG(p.km_pred_mM), 4) AS mean_km_pred,
            ROUND(MIN(p.km_pred_mM), 6) AS min_km_pred,
            ROUND(MAX(p.km_pred_mM), 3) AS max_km_pred,
            SUM(CASE WHEN s.reviewed=1 THEN 1 ELSE 0 END) AS n_reviewed,
            SUM(CASE WHEN s.km_best_mM IS NOT NULL THEN 1 ELSE 0 END) AS n_with_exp_km
        FROM predictions p
        JOIN sequences s ON s.id = p.sequence_id
        WHERE p.is_co2_pred = 1
        GROUP BY p.ec_pred
        ORDER BY n_predicted DESC
    """
    ec_dist = df_query(conn, ec_dist_sql)
    ec_dist["ec_name"] = ec_dist["ec_pred"].map(EC_NAMES).fillna("unknown")

    log.info("  EC class distribution (all 510K predicted carboxylases):")
    for _, row in ec_dist.iterrows():
        log.info("    %s (%s): n=%6d  mean_Km=%-8s  n_exp_Km=%d",
                 row["ec_pred"], row["ec_name"][:25],
                 row["n_predicted"],
                 f"{row['mean_km_pred']:.4f} mM" if row["mean_km_pred"] else "N/A",
                 row["n_with_exp_km"])

    # EC x organism source breakdown
    ec_source_sql = """
        SELECT
            p.ec_pred,
            s.source,
            COUNT(*) AS n
        FROM predictions p
        JOIN sequences s ON s.id = p.sequence_id
        WHERE p.is_co2_pred = 1
        GROUP BY p.ec_pred, s.source
        ORDER BY p.ec_pred, n DESC
    """
    ec_source = df_query(conn, ec_source_sql)

    # Km distribution per EC class
    km_per_ec_sql = """
        SELECT
            p.ec_pred,
            COUNT(*) AS n_with_km,
            ROUND(AVG(p.km_pred_mM), 4) AS mean_km,
            ROUND(AVG(p.km_pred_log10), 4) AS mean_log10_km,
            ROUND(MIN(p.km_pred_log10), 4) AS min_log10_km,
            ROUND(MAX(p.km_pred_log10), 4) AS max_log10_km
        FROM predictions p
        WHERE p.is_co2_pred = 1
          AND p.km_pred_mM IS NOT NULL
        GROUP BY p.ec_pred
        ORDER BY mean_km DESC
    """
    km_per_ec = df_query(conn, km_per_ec_sql)
    km_per_ec["ec_name"] = km_per_ec["ec_pred"].map(EC_NAMES).fillna("")

    log.info("  Km distribution per EC class:")
    for _, row in km_per_ec.iterrows():
        log.info("    %s (%s): n=%d  mean_Km=%.4f mM  range=10^%.2f to 10^%.2f",
                 row["ec_pred"], row["ec_name"][:20],
                 row["n_with_km"], row["mean_km"],
                 row["min_log10_km"], row["max_log10_km"])

    ec_dist.to_csv(FIG_DIR / "ec_distribution.tsv", sep="\t", index=False)
    km_per_ec.to_csv(FIG_DIR / "km_per_ec_distribution.tsv", sep="\t", index=False)
    ec_source.to_csv(FIG_DIR / "ec_by_source.tsv", sep="\t", index=False)

    result = {
        "task": "E_ec_distribution",
        "ec_distribution": ec_dist.to_dict("records"),
        "km_per_ec": km_per_ec.to_dict("records"),
    }
    json.dump(result, open(BIO_DIR / "ec_distribution.json", "w"), indent=2, cls=NumpyEncoder)
    log.info("  Saved: ec_distribution.json + 3 TSV files")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# TASK F: Top novel carboxylases
# ══════════════════════════════════════════════════════════════════════════════

def task_f_top_novel_carboxylases(conn, top_n=500):
    log.info("══ Task F: Top Novel Carboxylase Candidates ══")

    # Novel = unreviewed (TrEMBL) + high confidence + no experimental EC annotation
    novel_sql = """
        SELECT
            s.cdb_id, s.uniprot_id, s.organism, s.length,
            s.source, s.reviewed,
            s.ec_number AS ec_experimental,
            s.km_best_mM AS km_experimental,
            p.ec_pred, p.ec_prob, p.co2_prob,
            p.km_pred_mM, p.km_pred_log10,
            c.confidence_label, c.overall_score
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        JOIN confidence_scores c ON c.sequence_id = s.id
        WHERE p.is_co2_pred = 1
          AND s.reviewed = 0
          AND s.ec_number IS NULL
          AND c.confidence_label = 'high'
        ORDER BY p.co2_prob DESC
        LIMIT ?
    """
    novel = df_query(conn, novel_sql, (top_n,))
    novel["ec_name"] = novel["ec_pred"].map(EC_NAMES).fillna("")

    log.info("  Top novel carboxylases (unreviewed, no experimental EC, high confidence):")
    log.info("  Total: %d sequences", len(novel))
    log.info("  Top-10:")
    for _, row in novel.head(10).iterrows():
        log.info("    %s | %-40s | %s (%.4f) | Km=%s mM | conf=%.4f",
                 row["uniprot_id"],
                 row["organism"][:40] if row["organism"] else "unknown",
                 row["ec_pred"], row["ec_prob"],
                 f"{row['km_pred_mM']:.4f}" if row["km_pred_mM"] else "N/A",
                 row["co2_prob"])

    # Breakdown by EC class
    log.info("  Novel candidates by EC class:")
    for ec, grp in novel.groupby("ec_pred"):
        log.info("    %s (%s): n=%d  mean_prob=%.4f",
                 ec, EC_NAMES.get(ec, ""), len(grp), grp["co2_prob"].mean())

    # Novel candidates with no Pfam hits (truly novel domain architecture)
    novel_no_pfam_sql = """
        SELECT
            s.cdb_id, s.uniprot_id, s.organism,
            p.ec_pred, p.co2_prob, p.km_pred_mM,
            c.confidence_label
        FROM sequences s
        JOIN predictions p ON p.sequence_id = s.id
        JOIN confidence_scores c ON c.sequence_id = s.id
        LEFT JOIN features_domains fd ON fd.sequence_id = s.id
        WHERE p.is_co2_pred = 1
          AND c.confidence_label = 'high'
          AND fd.sequence_id IS NULL
        ORDER BY p.co2_prob DESC
        LIMIT 100
    """
    try:
        novel_no_pfam = df_query(conn, novel_no_pfam_sql)
        log.info("  High-confidence carboxylases with NO Pfam hits (novel families): n=%d", len(novel_no_pfam))
        if len(novel_no_pfam) > 0:
            novel_no_pfam.to_csv(FIG_DIR / "novel_no_pfam.tsv", sep="\t", index=False)
            for _, row in novel_no_pfam.head(5).iterrows():
                log.info("    %s | %s | %s | prob=%.4f",
                         row["uniprot_id"],
                         row["organism"][:40] if row["organism"] else "unknown",
                         row["ec_pred"], row["co2_prob"])
    except Exception as e:
        log.warning("  Could not query novel_no_pfam: %s", e)

    novel.to_csv(FIG_DIR / "top_novel_carboxylases.tsv", sep="\t", index=False)

    result = {
        "task": "F_top_novel_carboxylases",
        "n_novel_high_confidence": len(novel),
        "top10": novel.head(10)[["uniprot_id","organism","ec_pred","ec_prob",
                                  "co2_prob","km_pred_mM","confidence_label"]].to_dict("records"),
        "ec_breakdown": novel.groupby("ec_pred").size().to_dict(),
    }
    json.dump(result, open(BIO_DIR / "top_novel_carboxylases.json", "w"), indent=2, cls=NumpyEncoder)
    log.info("  Saved: top_novel_carboxylases.json + top_novel_carboxylases.tsv")
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(
        description="CarboDB v5 biological top-hits analysis.")
    ap.add_argument("--tasks", nargs="+",
                    default=["A", "B", "C", "D", "E", "F"],
                    choices=["A", "B", "C", "D", "E", "F"],
                    help="Tasks: A=controls, B=high-Km, C=low-Km, D=taxonomy, E=EC dist, F=novel")
    ap.add_argument("--top-n", type=int, default=500,
                    help="Number of top candidates to output (default: 500)")
    args = ap.parse_args()

    conn    = get_conn()
    tasks   = set(args.tasks)
    summary = {}

    if "A" in tasks:
        summary["A"] = task_a_validation_controls(conn)

    if "B" in tasks or "C" in tasks:
        summary["BC"] = task_bc_novel_km_candidates(conn, args.top_n)

    if "D" in tasks:
        summary["D"] = task_d_taxonomic_analysis(conn)

    if "E" in tasks:
        summary["E"] = task_e_ec_distribution(conn)

    if "F" in tasks:
        summary["F"] = task_f_top_novel_carboxylases(conn, args.top_n)

    conn.close()

    # ── Final summary ──────────────────────────────────────────────────────
    log.info("\n" + "=" * 70)
    log.info("BIOLOGICAL ANALYSIS SUMMARY — CarboDB v5")
    log.info("=" * 70)

    if "A" in summary:
        a = summary["A"]
        log.info("Validation controls: %d/%d in DB | %d/%d EC correct | %d/%d Km within 2-fold",
                 a["n_in_database"], a["n_controls"],
                 a["n_ec_correct"], a["n_in_database"],
                 a["n_km_within_2fold"], a["n_in_database"])

    if "BC" in summary:
        bc = summary["BC"]
        log.info("Novel high-Km: top candidate = %.3f mM | Novel low-Km: top = %.6f mM",
                 bc["high_km_max_mM"], bc["low_km_min_mM"])

    if "F" in summary:
        log.info("Novel high-confidence candidates: %d", summary["F"]["n_novel_high_confidence"])

    log.info("Outputs saved to: %s", BIO_DIR)
    log.info("Done. Next: python scripts/15_metagenome_analysis.py  OR  start webapp development")


if __name__ == "__main__":
    main()
