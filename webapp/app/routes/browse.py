from fastapi import APIRouter, Query
from typing import Optional
import sqlite3, os

router = APIRouter(tags=["browse"])

@router.get("/browse")
def browse(
    organism: Optional[str] = Query(None),
    ec: Optional[str] = Query(None),
    kingdom: Optional[str] = Query(None),
    km_min: Optional[float] = Query(None),
    km_max: Optional[float] = Query(None),
    reviewed: Optional[bool] = Query(None),
    limit: int = Query(50, le=500),
    offset: int = Query(0)
):
    db_path = os.environ.get("DB_PATH", "data/carbodb.sqlite")
    if not os.path.exists(db_path):
        return {"total": 0, "results": [], "error": "Database not found"}

    conn = sqlite3.connect(db_path, timeout=30)
    conn.row_factory = sqlite3.Row

    where = ["s.label=1", "p.km_pred_mM IS NOT NULL"]
    params = []

    if organism:
        where.append("s.organism LIKE ?")
        params.append(f"%{organism}%")
    if ec:
        where.append("s.ec_number LIKE ?")
        params.append(f"{ec}%")
    if kingdom:
        where.append("s.kingdom=?")
        params.append(kingdom)
    if km_min is not None:
        where.append("p.km_pred_mM*1000 >= ?")
        params.append(km_min)
    if km_max is not None:
        where.append("p.km_pred_mM*1000 <= ?")
        params.append(km_max)
    if reviewed is not None:
        where.append("s.reviewed=?")
        params.append(1 if reviewed else 0)

    where_str = " AND ".join(where)

    try:
        count_row = conn.execute(
            "SELECT COUNT(*) FROM sequences s "
            "JOIN predictions p ON p.sequence_id=s.id "
            "WHERE " + where_str,
            params).fetchone()
        total = count_row[0]

        rows = conn.execute(
            "SELECT s.uniprot_id, s.organism, s.ec_number, s.length, "
            "s.reviewed, s.source, p.km_pred_mM*1000 as km_uM "
            "FROM sequences s "
            "JOIN predictions p ON p.sequence_id=s.id "
            "WHERE " + where_str + " "
            "ORDER BY p.km_pred_mM "
            "LIMIT ? OFFSET ?",
            params + [limit, offset]).fetchall()

        conn.close()
        return {"total": total, "limit": limit, "offset": offset,
                "results": [dict(r) for r in rows]}
    except Exception as e:
        conn.close()
        return {"total": 0, "results": [], "error": str(e)}


@router.get("/stats")
def stats():
    db_path = os.environ.get("DB_PATH", "data/carbodb.sqlite")
    if not os.path.exists(db_path):
        return {"error": "Database not found"}
    try:
        conn = sqlite3.connect(db_path, timeout=10)
        total = conn.execute(
            "SELECT COUNT(*) FROM sequences WHERE label=1").fetchone()[0]
        reviewed = conn.execute(
            "SELECT COUNT(*) FROM sequences WHERE label=1 AND reviewed=1").fetchone()[0]
        ec_dist = conn.execute(
            "SELECT ec_number, COUNT(*) as n FROM sequences "
            "WHERE label=1 GROUP BY ec_number ORDER BY n DESC LIMIT 10").fetchall()
        conn.close()
        return {"total_sequences": total, "reviewed": reviewed,
                "ec_distribution": {r[0]: r[1] for r in ec_dist}}
    except Exception as e:
        return {"error": str(e)}
