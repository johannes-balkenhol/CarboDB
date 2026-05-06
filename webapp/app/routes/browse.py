"""
Database browser endpoints.

Three endpoints:
- GET /browse              search/filter/paginate stored sequences
- GET /stats               aggregate counters for the stats banner
- GET /db/seq/{uniprot_id} ResultDetail-shaped JSON for a single stored sequence
"""
from fastapi import APIRouter, HTTPException, Query
from typing import Optional
import sqlite3, os, json

from ..startup import EC_NAMES

router = APIRouter(tags=["browse"])

# Stats cache — banner numbers don't change between sessions, only when ingestion runs.
import time as _time
_STATS_CACHE = {"data": None, "ts": 0.0}
_STATS_TTL_SEC = 600  # 10 minutes

DB_PATH_ENV = "DB_PATH"
DEFAULT_DB = "data/primary/carbodb.sqlite"


def _db():
    """Open a read-only SQLite connection. Caller must close."""
    path = os.environ.get(DB_PATH_ENV, DEFAULT_DB)
    if not os.path.exists(path):
        raise HTTPException(503, f"Database not found at {path}")
    conn = sqlite3.connect(path, timeout=30)
    conn.row_factory = sqlite3.Row
    # Performance pragmas: 512MB page cache, in-memory temp btrees
    conn.execute("PRAGMA cache_size = -524288")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("PRAGMA mmap_size = 30000000000")
    return conn


# ────────────────────────────────────────────────────────────────────────────
# /browse — list/search/filter/paginate
# ────────────────────────────────────────────────────────────────────────────
@router.get("/browse")
def browse(
    q: Optional[str] = Query(None, description="Match uniprot_id or organism (LIKE)"),
    ec: Optional[str] = Query(None, description="EC number prefix, e.g. 4.2.1.1"),
    is_carboxylase: Optional[bool] = Query(None, description="Filter by predicted carboxylase status"),
    has_experimental_km: Optional[bool] = Query(None, description="Only sequences with measured Km"),
    reviewed: Optional[bool] = Query(None, description="SwissProt-curated entries only"),
    min_km_uM: Optional[float] = Query(None),
    max_km_uM: Optional[float] = Query(None),
    sort: str = Query("default",
                      description="default | km_asc | km_desc | length_asc | length_desc | uniprot"),
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
):
    """Search and filter the carboxylase database.

    Returns: { total, limit, offset, results: [...] }
    Each result row has the columns the frontend table renders.
    """
    where = ["s.label = 1"]   # only show known carboxylases (positives)
    params = []

    if q:
        where.append("(s.uniprot_id LIKE ? OR s.organism LIKE ?)")
        params.extend([f"%{q}%", f"%{q}%"])
    if ec:
        where.append("s.ec_number LIKE ?")
        params.append(f"{ec}%")
    if is_carboxylase is not None:
        if is_carboxylase:
            where.append("p.is_co2_pred = 1")
        else:
            where.append("(p.is_co2_pred = 0 OR p.is_co2_pred IS NULL)")
    if reviewed is not None:
        where.append("s.reviewed = ?")
        params.append(1 if reviewed else 0)
    if has_experimental_km:
        where.append("s.km_best_mM IS NOT NULL")
    if min_km_uM is not None:
        where.append("p.km_pred_mM * 1000 >= ?")
        params.append(min_km_uM)
    if max_km_uM is not None:
        where.append("p.km_pred_mM * 1000 <= ?")
        params.append(max_km_uM)

    order_by = {
        "default":     "s.km_best_mM IS NULL, p.km_pred_mM IS NULL, p.km_pred_mM",
        "km_asc":      "p.km_pred_mM IS NULL, p.km_pred_mM ASC",
        "km_desc":     "p.km_pred_mM IS NULL, p.km_pred_mM DESC",
        "km_exp_asc":  "s.km_best_mM IS NULL, s.km_best_mM ASC",
        "km_exp_desc": "s.km_best_mM IS NULL, s.km_best_mM DESC",
        "length_asc":  "s.length ASC",
        "length_desc": "s.length DESC",
        "uniprot":     "s.uniprot_id ASC",
    }.get(sort, "s.km_best_mM IS NULL, p.km_pred_mM IS NULL, p.km_pred_mM")

    where_str = " AND ".join(where)
    base_from = ("FROM sequences s "
                 "LEFT JOIN predictions p ON p.sequence_id = s.id "
                 "                       AND p.model_version = 'v5' "
                 f"WHERE {where_str}")

    conn = _db()
    try:
        total = conn.execute(f"SELECT COUNT(*) {base_from}", params).fetchone()[0]
        rows = conn.execute(
            f"SELECT s.uniprot_id, s.organism, s.ec_number AS ec_known, "
            f"       p.ec_pred AS ec_predicted, p.ec_prob AS ec_confidence, "
            f"       s.length, s.reviewed, s.source, "
            f"       p.km_pred_mM * 1000 AS km_predicted_uM, "
            f"       s.km_best_mM   * 1000 AS km_experimental_uM, "
            f"       p.is_co2_pred  AS is_carboxylase, "
            f"       p.co2_prob     AS carboxylase_probability "
            f"{base_from} "
            f"ORDER BY {order_by} "
            f"LIMIT ? OFFSET ?",
            params + [limit, offset]
        ).fetchall()
    finally:
        conn.close()

    out = []
    for r in rows:
        d = dict(r)
        # add the human-readable EC name; favour predicted, else known
        ec_use = d.get("ec_predicted") or d.get("ec_known")
        d["ec_name"] = EC_NAMES.get(ec_use, "")
        out.append(d)

    return {"total": total, "limit": limit, "offset": offset, "results": out}


# ────────────────────────────────────────────────────────────────────────────
# /stats — banner numbers
# ────────────────────────────────────────────────────────────────────────────
@router.get("/stats")
def stats():
    now = _time.time()
    cached = _STATS_CACHE["data"]
    if cached is not None and now - _STATS_CACHE["ts"] < _STATS_TTL_SEC:
        return cached
    conn = _db()
    try:
        # Single pass over sequences (uses idx_seq_label) instead of 4 passes
        agg = conn.execute(
            "SELECT COUNT(*) AS total, "
            "       SUM(CASE WHEN km_best_mM IS NOT NULL THEN 1 ELSE 0 END) AS with_km, "
            "       SUM(CASE WHEN reviewed = 1 THEN 1 ELSE 0 END) AS reviewed_n "
            "FROM sequences WHERE label = 1"
        ).fetchone()
        total = agg["total"]
        with_experimental_km = agg["with_km"] or 0
        reviewed_count = agg["reviewed_n"] or 0

        # Use distinct sequence_id since predictions table can have multi-version rows
        predicted_carboxylases = conn.execute(
            "SELECT COUNT(DISTINCT sequence_id) FROM predictions "
            "WHERE is_co2_pred = 1 AND model_version = 'v5'"
        ).fetchone()[0]

        ec_dist_rows = conn.execute(
            "SELECT ec_number, COUNT(*) AS n "
            "FROM sequences WHERE label=1 "
            "GROUP BY ec_number ORDER BY n DESC"
        ).fetchall()
    finally:
        conn.close()

    ec_distribution = [
        {"ec_number": r[0],
         "ec_name":   EC_NAMES.get(r[0], r[0]),
         "count":     r[1]}
        for r in ec_dist_rows
    ]

    result = {
        "total_sequences":       total,
        "predicted_carboxylases": predicted_carboxylases,
        "with_experimental_km":   with_experimental_km,
        "reviewed_count":         reviewed_count,
        "ec_classes_total":       len(ec_distribution),
        "ec_distribution":        ec_distribution,
    }
    _STATS_CACHE["data"] = result
    _STATS_CACHE["ts"] = now
    return result


# ────────────────────────────────────────────────────────────────────────────
# /db/seq/{uniprot_id} — ResultDetail-shaped JSON, instant from DB
# ────────────────────────────────────────────────────────────────────────────
@router.get("/db/seq/{uniprot_id}")
def get_db_seq(uniprot_id: str):
    """Return a precomputed prediction in /predict-shaped JSON.

    The shape matches what predict_sequence() returns so the frontend
    can reuse the existing <ResultDetail> component without changes.
    Note: SHAP and live BLAST nearest-neighbors are NOT included; this
    is the cheap browser path. Recompute via /predict if those are needed.
    """
    conn = _db()
    try:
        seq = conn.execute(
            "SELECT s.id, s.uniprot_id, s.cdb_id, s.ec_number, "
            "       s.length, s.organism, s.reviewed, s.source, "
            "       s.km_best_mM, s.sequence "
            "FROM sequences s WHERE s.uniprot_id = ?",
            (uniprot_id,)
        ).fetchone()
        if not seq:
            raise HTTPException(404, f"Sequence {uniprot_id} not found")

        # latest prediction
        pred = conn.execute(
            "SELECT * FROM predictions "
            "WHERE sequence_id = ? "
            "ORDER BY created_at DESC LIMIT 1",
            (seq["id"],)
        ).fetchone()

        # Pfam domains
        dom = conn.execute(
            "SELECT pfam_hits_json FROM features_domains WHERE sequence_id = ?",
            (seq["id"],)
        ).fetchone()
        pfam_hits = []
        if dom and dom["pfam_hits_json"]:
            try:
                raw = json.loads(dom["pfam_hits_json"])
                # raw may be list[str] or list[dict]; normalise
                for h in raw:
                    if isinstance(h, str):
                        pfam_hits.append({"accession": h})
                    elif isinstance(h, dict):
                        pfam_hits.append(h)
            except json.JSONDecodeError:
                pass

        # composition feature row (motifs + physicochem + AAC + dipeptide)
        comp = conn.execute(
            "SELECT * FROM features_composition WHERE sequence_id = ?",
            (seq["id"],)
        ).fetchone()
        features_computed = {}
        if comp:
            for k in comp.keys():
                if k.startswith(("aac_", "phys_", "inv_", "motif_", "dp_", "pse_")):
                    features_computed[k] = comp[k]

        # experimental Km (top-tier only)
        km_evs = conn.execute(
            "SELECT km_value_mM, substrate, source, evidence_tier, commentary "
            "FROM km_evidence "
            "WHERE sequence_id = ? "
            "ORDER BY evidence_tier ASC, km_value_mM ASC",
            (seq["id"],)
        ).fetchall()
        # Build a single "self-hit" entry so ResultDetail's neighbor card renders
        # cleanly with the experimental Km when one exists.
        top_similar = []
        if km_evs:
            ev = km_evs[0]
            top_similar.append({
                "uniprot_id":         seq["uniprot_id"],
                "organism":           seq["organism"] or "",
                "identity_pct":       100.0,
                "evalue":             0,
                "align_length":       seq["length"],
                "km_experimental_uM": ev["km_value_mM"] * 1000,
                "km_exp_substrate":   ev["substrate"] or "",
                "tier":               ev["source"],
            })

        # Compose response
        ec_pred = pred["ec_pred"] if pred else None
        return {
            "id":                       seq["uniprot_id"],
            "cdb_query_id":             seq["cdb_id"],
            "sequence_length":          seq["length"],
            "is_carboxylase":           bool(pred and pred["is_co2_pred"]),
            "carboxylase_probability":  (pred["co2_prob"] if pred else None),
            "ec_predicted":             ec_pred,
            "ec_name":                  EC_NAMES.get(ec_pred, ""),
            "ec_confidence":            (pred["ec_prob"] if pred else None),
            "km_predicted_mM":          (pred["km_pred_mM"] if pred else None),
            "km_predicted_uM":          (pred["km_pred_mM"] * 1000 if pred and pred["km_pred_mM"] is not None else None),
            "km_predicted_log10":       (pred["km_pred_log10"] if pred else None),
            "pfam_hits":                pfam_hits,
            "features_computed":        features_computed,
            "top_similar":              top_similar,
            "shap":                     {},   # not precomputed in DB; see /predict for live
            "ec_probabilities":         {},   # likewise
            "novelty_flag":             "known",
            "mode":                     "db_lookup",
            "kingdom":                  "n/a",
            "runtime_seconds":          0.0,
            "model_version":            (pred["model_version"] if pred else None),
            "organism":                 seq["organism"] or "",
            "reviewed":                 bool(seq["reviewed"]),
            "source":                   seq["source"],
        }
    finally:
        conn.close()
