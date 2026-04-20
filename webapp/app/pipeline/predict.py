
import os, sys, time, logging, subprocess, json, tempfile
from pathlib import Path
import numpy as np

from ..startup import ModelStore, EC_NAMES, KM_EC_CLASSES

log = logging.getLogger(__name__)

# Use local copy of annotation script (standalone, no CarboDB repo needed)
_local = Path(__file__).parent / "annotate.py"
_parent_scripts = Path(__file__).resolve()
for _ in range(6):
    _parent_scripts = _parent_scripts.parent
    candidate = _parent_scripts / "scripts" / "11_annotate_sequence.py"
    if candidate.exists():
        break
SCRIPT_11 = _local if _local.exists() else candidate


def predict_sequence(sequence, mode="fast", kingdom="plant", seq_id="query"):
    t = time.time()

    # Write sequence to temp FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".faa",
                                     delete=False) as f:
        f.write(f">{seq_id}\n{sequence.strip()}\n")
        fasta_path = f.name

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        out_path = f.name

    try:
        cmd = [sys.executable, str(SCRIPT_11),
               "--fasta", fasta_path,
               "--out", out_path,
               "--kingdom", kingdom]

        if mode in ("fast", "pfam"):
            cmd.append("--no-esm2")

        result = subprocess.run(cmd, capture_output=True, timeout=120,
                                env={**os.environ,
                                     "PFAM_HMM": os.environ.get("PFAM_HMM", "data/Pfam-A.hmm")})

        if result.returncode != 0:
            raise ValueError(result.stderr.decode()[:500])

        raw = json.load(open(out_path))
        d = raw[0] if isinstance(raw, list) else raw

        # Normalise output to webapp format
        return {
            "is_carboxylase": d.get("is_carboxylase", False),
            "carboxylase_probability": d.get("carboxylase_probability", 0.0),
            "ec_predicted": d.get("ec_predicted", "unknown"),
            "ec_name": d.get("ec_name", ""),
            "ec_confidence": d.get("ec_probabilities", {}).get(d.get("ec_predicted",""), 0.0),
            "ec_probabilities": d.get("ec_probabilities", {}),
            "km_predicted_mM": d.get("km_predicted_mM"),
            "km_predicted_uM": d.get("km_predicted_mM") * 1000 if d.get("km_predicted_mM") else None,
            "sequence_length": d.get("sequence_length", 0),
            "pfam_hits": d.get("pfam_hits", []),
            "novelty_flag": "known" if d.get("pfam_hits") and d.get("carboxylase_probability",0) > 0.8 else
                            "borderline" if d.get("carboxylase_probability",0) > 0.5 else "novel",
            "features_used": d.get("features_used", []),
            "mode": mode,
            "kingdom": kingdom,
            "runtime_seconds": round(time.time() - t, 2),
            "top_similar": [],
        }
    finally:
        for p in [fasta_path, out_path]:
            try: os.unlink(p)
            except: pass
