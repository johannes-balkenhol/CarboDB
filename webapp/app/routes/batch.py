from fastapi import APIRouter, UploadFile, File, Form, BackgroundTasks, HTTPException
from fastapi.responses import FileResponse
from typing import Optional
import uuid, os, json
from datetime import datetime

router = APIRouter(tags=["batch"])

JOBS_DIR = os.environ.get("JOBS_DIR", "jobs")
MAX_FAST = int(os.environ.get("MAX_BATCH_FAST", 5000))
MAX_STANDARD = int(os.environ.get("MAX_BATCH_STANDARD", 500))


@router.post("/batch")
async def submit_batch(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    mode: str = Form("fast"),
    kingdom: str = Form("plant"),
    email: Optional[str] = Form(None)
):
    if mode not in ('fast', 'standard'):
        raise HTTPException(400, "Batch mode must be fast or standard")
    content = await file.read()
    text = content.decode('utf-8', errors='ignore')
    n_seqs = text.count('>')
    if n_seqs == 0:
        raise HTTPException(400, "No sequences found in FASTA file")
    max_seqs = MAX_STANDARD if mode == 'standard' else MAX_FAST
    if n_seqs > max_seqs:
        raise HTTPException(400, f"Too many sequences: {n_seqs} > max {max_seqs}")
    job_id = str(uuid.uuid4())[:8]
    job_dir = os.path.join(JOBS_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)
    input_path = os.path.join(job_dir, "input.faa")
    with open(input_path, 'w') as f:
        f.write(text)
    meta = {
        "job_id": job_id, "status": "queued", "mode": mode,
        "kingdom": kingdom, "email": email, "n_sequences": n_seqs,
        "processed": 0, "created_at": datetime.utcnow().isoformat(),
        "input_file": input_path,
        "result_file": os.path.join(job_dir, "results.tsv")
    }
    with open(os.path.join(job_dir, "job.json"), 'w') as f:
        json.dump(meta, f)
    background_tasks.add_task(run_batch_job, job_id, input_path, mode, kingdom)
    est_min = round(n_seqs * (45 if mode == 'standard' else 3) / 60)
    return {"job_id": job_id, "status": "queued",
            "n_sequences": n_seqs, "estimated_minutes": est_min}


@router.get("/jobs/{job_id}")
def get_job(job_id: str):
    job_path = os.path.join(JOBS_DIR, job_id, "job.json")
    if not os.path.exists(job_path):
        raise HTTPException(404, "Job not found")
    with open(job_path) as f:
        meta = json.load(f)
    n = meta.get('n_sequences', 1)
    meta['progress_pct'] = round(meta.get('processed', 0) / n * 100) if n > 0 else 0
    return meta


@router.get("/jobs/{job_id}/results.tsv")
def download_results(job_id: str):
    job_path = os.path.join(JOBS_DIR, job_id, "job.json")
    if not os.path.exists(job_path):
        raise HTTPException(404, "Job not found")
    with open(job_path) as f:
        meta = json.load(f)
    if meta.get('status') != 'completed':
        raise HTTPException(400, f"Job not completed. Status: {meta.get('status')}")
    result_file = meta.get('result_file')
    if not result_file or not os.path.exists(result_file):
        raise HTTPException(404, "Results file not found")
    return FileResponse(result_file, media_type='text/tab-separated-values',
                        filename=f"carbodb_batch_{job_id}.tsv")


@router.get("/jobs/{job_id}/seq/{seq_id}")
def get_seq_detail(job_id: str, seq_id: str):
    """Return full per-sequence predict response saved by run_batch_job."""
    job_path = os.path.join(JOBS_DIR, job_id, "job.json")
    if not os.path.exists(job_path):
        raise HTTPException(404, "Job not found")
    seq_json = os.path.join(JOBS_DIR, job_id, f"seq_{seq_id}.json")
    if not os.path.exists(seq_json):
        raise HTTPException(404, "Per-sequence detail not found")
    with open(seq_json) as f:
        return json.load(f)


def run_batch_job(job_id: str, input_path: str, mode: str, kingdom: str):
    from ..pipeline.predict import predict_sequence
    from ..pipeline.blast_similar import run_blast_similar
    job_dir = os.path.join(JOBS_DIR, job_id)
    meta_path = os.path.join(job_dir, "job.json")

    def update_meta(updates):
        with open(meta_path) as f:
            meta = json.load(f)
        meta.update(updates)
        with open(meta_path, 'w') as f:
            json.dump(meta, f)

    update_meta({"status": "running", "started_at": datetime.utcnow().isoformat()})
    result_path = os.path.join(job_dir, "results.tsv")
    header = ("seq_id\tlength\tis_carboxylase\tprob_binary\tec_predicted\tec_confidence\t" "km_predicted_mM\tkm_predicted_uM\tpfam_hits\tnovelty_flag\truntime_seconds\t" "nearest_uniprot\tnearest_pident\tnearest_organism\tnearest_km_exp_uM\tnearest_tier\n")
    processed = 0
    try:
        with open(input_path) as fin, open(result_path, 'w') as fout:
            fout.write(header)
            seqs = {}
            sid, buf = None, []
            for line in fin:
                line = line.strip()
                if line.startswith('>'):
                    if sid:
                        seqs[sid] = ''.join(buf)
                    sid = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line)
            if sid:
                seqs[sid] = ''.join(buf)
            for seq_id, sequence in seqs.items():
                try:
                    r = predict_sequence(sequence, mode=mode,
                                        kingdom=kingdom, seq_id=seq_id)
                    # Save full per-sequence response so the UI Details panel can load
                    # rich data (SHAP, features_computed, top_similar) without re-running predict.
                    try:
                        seq_json_path = os.path.join(job_dir, f"seq_{seq_id}.json")
                        with open(seq_json_path, "w") as fseq:
                            json.dump(r, fseq)
                    except Exception as exc:
                        print(f"Failed to save per-seq JSON for {seq_id}: {exc}")
                    # pfam_hits is now list[dict] with {accession, name, e_value, bitscore}.
                    # Collapse to accessions-only for the TSV (rich data lives in single-predict).
                    _hits = r.get('pfam_hits', []) or []
                    pfam_str = ';'.join(
                        h.get('accession', '') if isinstance(h, dict) else str(h)
                        for h in _hits
                    )
                    # Nearest-neighbor via BLAST against experimental-Km DB (top 1 only)
                    nn_uid = nn_pident = nn_org = nn_km = nn_tier = ''
                    if r.get('ec_predicted') and r.get('is_carboxylase'):
                        try:
                            nn = run_blast_similar(
                                sequence=sequence,
                                ec_predicted=r['ec_predicted'],
                                limit=1,
                                manifest_path="data/blast_ec_dbs_exp/manifest.json",
                            )
                            if nn:
                                h0 = nn[0]
                                nn_uid = h0.get('uniprot_id', '') or ''
                                nn_pident = f"{h0.get('identity_pct', ''):.1f}" if h0.get('identity_pct') is not None else ''
                                nn_org = (h0.get('organism') or '').replace('\t', ' ')
                                nn_km = h0.get('km_experimental_uM')
                                nn_km = f"{nn_km:.2f}" if nn_km is not None else ''
                                nn_tier = h0.get('tier', '') or ''
                        except Exception as exc:
                            print(f"BLAST failed for {seq_id}: {exc}")
                    km_mM_str = f"{r['km_predicted_mM']:.4f}" if r.get('km_predicted_mM') is not None else ''
                    km_uM_str = f"{r['km_predicted_uM']:.2f}" if r.get('km_predicted_uM') is not None else ''
                    fout.write(
                        f"{seq_id}\t{r['sequence_length']}\t"
                        f"{r['is_carboxylase']}\t{r['carboxylase_probability']:.4f}\t"
                        f"{r['ec_predicted']}\t{r['ec_confidence']:.4f}\t"
                        f"{km_mM_str}\t{km_uM_str}\t"
                        f"{pfam_str}\t{r.get('novelty_flag','')}\t"
                        f"{r.get('runtime_seconds','')}\t"
                        f"{nn_uid}\t{nn_pident}\t{nn_org}\t{nn_km}\t{nn_tier}\n"
                    )
                    fout.flush()
                except Exception as e:
                    fout.write(f"{seq_id}\t0\tERROR\t\t\t\t\t\t\t\t{str(e)[:100]}\n")
                processed += 1
                if processed % 10 == 0:
                    update_meta({"processed": processed})
        update_meta({"status": "completed", "processed": processed,
                     "completed_at": datetime.utcnow().isoformat(),
                     "result_file": result_path})
    except Exception as e:
        update_meta({"status": "failed", "error_message": str(e)})
