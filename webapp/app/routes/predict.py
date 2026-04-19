"""
routes/predict.py — Single sequence prediction endpoint
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional
import sqlite3, os, json

from ..pipeline.predict import predict_sequence
from ..startup import ModelStore

router = APIRouter(tags=["predict"])

class PredictRequest(BaseModel):
    sequence: str
    mode: str = "fast"        # fast | standard | pfam | composite
    kingdom: str = "plant"    # bacteria | plant | archaea | fungi
    seq_id: Optional[str] = "query"

@router.post("/predict")
def predict(req: PredictRequest):
    if not ModelStore.ready:
        raise HTTPException(503, "Models not loaded yet")
    if not req.sequence or len(req.sequence.strip()) < 10:
        raise HTTPException(400, "Sequence too short")
    if req.mode not in ('fast', 'standard', 'pfam', 'composite'):
        raise HTTPException(400, f"Invalid mode: {req.mode}")

    try:
        result = predict_sequence(
            sequence=req.sequence,
            mode=req.mode,
            kingdom=req.kingdom,
            seq_id=req.seq_id or "query"
        )
        # Add similar sequences from DB
        if result.get('ec_predicted') and result.get('is_carboxylase'):
            result['top_similar'] = get_similar_from_db(
                result['ec_predicted'], result.get('km_predicted_uM'))
        else:
            result['top_similar'] = []
        return result
    except ValueError as e:
        raise HTTPException(400, str(e))
    except Exception as e:
        raise HTTPException(500, f"Prediction failed: {e}")


def get_similar_from_db(ec: str, km_uM: Optional[float], limit: int = 8) -> list:
    """Get similar sequences from CarboDB for context."""
    db_path = os.environ.get("DB_PATH", "data/carbodb.sqlite")
    if not os.path.exists(db_path):
        return []
    try:
        conn = sqlite3.connect(db_path, timeout=10)
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute("""
            SELECT s.uniprot_id, s.organism, p.km_pred_mM*1000 as km_pred_uM,
                   s.reviewed
            FROM sequences s
            JOIN predictions p ON p.sequence_id = s.id
            WHERE s.label=1 AND s.ec_number=?
            AND p.km_pred_mM IS NOT NULL AND s.reviewed=1
            ORDER BY RANDOM()
            LIMIT ?
        """, (ec, limit))
        rows = cur.fetchall()
        conn.close()
        return [{'uniprot_id': r['uniprot_id'],
                 'organism': r['organism'],
                 'km_predicted_uM': round(r['km_pred_uM'], 1) if r['km_pred_uM'] else None,
                 'km_experimental_uM': None,
                 'reviewed': bool(r['reviewed'])} for r in rows]
    except Exception as e:
        return []