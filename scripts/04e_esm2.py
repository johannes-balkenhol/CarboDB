#!/usr/bin/env python3
"""
04e_esm2.py
===========
CarboxyDB — ESM-2 embeddings on A100 GPU (wbbi206 login node).

Run directly in a screen session — NOT via SLURM (A100 is on login node).

Model: esm2_t33_650M_UR50D  (1280-dim, 650M parameters)
Input: data/primary/master.fasta
Output: data/features/esm2/esm2_{chunk:04d}.tsv  (one file per 1000 seqs)

Usage:
    screen -S esm2
    cd ~/Projects_shared/CarboDB_v3
    python scripts/04e_esm2.py
    # Ctrl+A D to detach

    # Monitor progress:
    ls data/features/esm2/ | wc -l   # target: 2380 files
    tail -f logs/04e_esm2_*.log

Resume:
    Already-processed chunks are skipped automatically.
    Just rerun the same command.

Expected runtime: ~30h on A100 80GB at batch_size=64
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, setup_logging

log = setup_logging("04e_esm2")

BATCH_SIZE    = 64     # safe for A100 80GB
MAX_SEQ_LEN   = 1022   # ESM-2 token limit
CHUNK_SIZE    = 1000   # sequences per output file


def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs = []
    cdb_id = seq_lines = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cdb_id and seq_lines:
                    seqs.append((cdb_id, "".join(seq_lines)))
                cdb_id    = line[1:].split("|")[0]
                seq_lines = []
            elif line:
                if seq_lines is not None:
                    seq_lines.append(line)
    if cdb_id and seq_lines:
        seqs.append((cdb_id, "".join(seq_lines)))
    return seqs


def load_model():
    log.info("Loading ESM-2 model: %s", CFG.ESM2_MODEL)
    import esm
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval()
    if torch.cuda.is_available():
        model = model.cuda()
        log.info("  GPU: %s (%.1f GB VRAM)",
                 torch.cuda.get_device_name(0),
                 torch.cuda.get_device_properties(0).total_memory / 1e9)
    else:
        log.warning("  No GPU found — running on CPU (will be very slow)")
    batch_converter = alphabet.get_batch_converter()
    return model, alphabet, batch_converter


def embed_batch(model, batch_converter, seqs_batch: list[tuple[str,str]],
                device: str) -> np.ndarray:
    """
    Embed a batch of sequences.
    Returns numpy array of shape (n_seqs, 1280).
    """
    # Truncate long sequences
    data = [(cdb_id, seq[:MAX_SEQ_LEN]) for cdb_id, seq in seqs_batch]

    _, _, tokens = batch_converter(data)
    tokens = tokens.to(device)

    with torch.no_grad():
        results = model(tokens, repr_layers=[33], return_contacts=False)

    # Mean-pool over sequence positions (exclude BOS/EOS tokens)
    token_representations = results["representations"][33]
    embeddings = []
    for i, (_, seq) in enumerate(data):
        # Tokens: [BOS] seq [EOS] [PAD...]
        seq_len = min(len(seq), MAX_SEQ_LEN)
        emb = token_representations[i, 1:seq_len+1].mean(dim=0)
        embeddings.append(emb.cpu().float().numpy())

    return np.array(embeddings)


def process_chunk(chunk_seqs: list[tuple[str,str]], chunk_id: int,
                  model, batch_converter, device: str) -> pd.DataFrame:
    rows = []
    n_batches = (len(chunk_seqs) + BATCH_SIZE - 1) // BATCH_SIZE

    for b in range(n_batches):
        batch = chunk_seqs[b*BATCH_SIZE : (b+1)*BATCH_SIZE]
        try:
            embs = embed_batch(model, batch_converter, batch, device)
            for (cdb_id, _), emb in zip(batch, embs):
                row = {"cdb_id": cdb_id}
                for j, v in enumerate(emb):
                    row[f"esm2_{j}"] = float(v)
                rows.append(row)
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                log.error("OOM on batch %d/%d chunk %04d — reduce BATCH_SIZE",
                          b+1, n_batches, chunk_id)
                torch.cuda.empty_cache()
                # Process one by one as fallback
                for item in batch:
                    try:
                        embs = embed_batch(model, batch_converter, [item], device)
                        cdb_id = item[0]
                        row = {"cdb_id": cdb_id}
                        for j, v in enumerate(embs[0]):
                            row[f"esm2_{j}"] = float(v)
                        rows.append(row)
                    except Exception as e2:
                        log.warning("  Skip %s: %s", item[0], e2)
            else:
                log.error("Error in batch %d: %s", b+1, e)

    return pd.DataFrame(rows)


def main():
    PATHS.FEAT_ESM2.mkdir(parents=True, exist_ok=True)

    # Load sequences
    log.info("Reading master.fasta...")
    all_seqs = read_fasta(PATHS.MASTER_FASTA)
    log.info("  %d sequences total", len(all_seqs))

    # Load model
    model, alphabet, batch_converter = load_model()
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Split into chunks
    n_chunks = (len(all_seqs) + CHUNK_SIZE - 1) // CHUNK_SIZE
    log.info("Processing %d chunks of %d sequences", n_chunks, CHUNK_SIZE)

    t_start = time.time()
    done = 0

    for chunk_id in range(1, n_chunks + 1):
        out_path = PATHS.FEAT_ESM2 / f"esm2_{chunk_id:04d}.tsv"

        if out_path.exists():
            done += 1
            continue

        chunk_seqs = all_seqs[(chunk_id-1)*CHUNK_SIZE : chunk_id*CHUNK_SIZE]

        df = process_chunk(chunk_seqs, chunk_id, model, batch_converter, device)

        if not df.empty:
            df.to_csv(out_path, sep="\t", index=False)

        done += 1
        elapsed = time.time() - t_start
        rate    = done / elapsed * 3600  # chunks per hour
        eta_h   = (n_chunks - done) / max(rate / 3600, 1e-6) / 3600
        log.info("  Chunk %04d/%04d done | %.0f chunks/h | ETA %.1fh",
                 chunk_id, n_chunks, rate, eta_h)

        # Clear GPU cache periodically
        if chunk_id % 10 == 0 and torch.cuda.is_available():
            torch.cuda.empty_cache()

    log.info("ESM-2 embedding complete. %d chunks written.", done)
    log.info("Next: python scripts/04_merge_features.py")


if __name__ == "__main__":
    main()
