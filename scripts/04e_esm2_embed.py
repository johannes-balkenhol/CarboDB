#!/usr/bin/env python3
"""
04e_esm2_embed.py -- CarboDB v5.2: ESM-2 embeddings (esm2_t33_650M_UR50D, layer 33, mean over residues)

Two inputs:
  --missing-from-db DB     sequences that have no features_esm2 row in the v2.2 DB (the 894)
  --fasta FILE             any FASTA (e.g. GMGC candidates); IDs = first word of the header
Output (stored as files, not in SQL): <out>/esm2.npy (float32, N x 1280, memory-mapped),
  <out>/row_index.tsv (row -> id), <out>/run.json ; written shard by shard, resumable.

SAFETY: before embedding anything, --check-db DB re-embeds 40 sequences that already have a v5
embedding (20 random + 20 longer than 1022 aa) and compares. If the minimum cosine < 0.9999 the
run stops: the new vectors would not be compatible with v5.

Half the machine: torch uses 2 CPU threads, GPU batch limited by --toks (default 12000 ~ <40 GB).
"""
import argparse, hashlib, json, os, random, socket, sqlite3, sys, time
import numpy as np

def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)

def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line: seq.append(line)
    if name: yield name, "".join(seq)

class Embedder:
    def __init__(self, toks, trunc):
        import torch, esm
        torch.set_num_threads(2)
        self.torch = torch
        self.model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.model.eval().cuda()
        self.bc = self.alphabet.get_batch_converter()
        self.toks, self.trunc = toks, trunc
    def __call__(self, items):                     # items: list of (id, seq)
        out = {}
        items = sorted(items, key=lambda x: len(x[1]))
        batch = []
        def run(b):
            _, _, t = self.bc([(i, s[:self.trunc]) for i, s in b])
            with self.torch.no_grad():
                r = self.model(t.cuda(), repr_layers=[33])["representations"][33]
            for k, (i, s) in enumerate(b):
                L = min(len(s), self.trunc)
                out[i] = r[k, 1:L + 1].mean(0).float().cpu().numpy()
        for it in items:                           # sorted by length: the current item is the longest
            L = min(len(it[1]), self.trunc) + 2
            if batch and L * (len(batch) + 1) > self.toks:   # padded tokens if we add it
                run(batch); batch = []
            batch.append(it)
        if batch: run(batch)
        return out

def check_against_db(emb, db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    rows = con.execute("SELECT s.cdb_id, s.sequence, e.embedding_blob FROM features_esm2 e "
                       "JOIN sequences s ON s.id=e.sequence_id WHERE s.length<=1022 "
                       "ORDER BY random() LIMIT 20").fetchall()
    rows += con.execute("SELECT s.cdb_id, s.sequence, e.embedding_blob FROM features_esm2 e "
                        "JOIN sequences s ON s.id=e.sequence_id WHERE s.length>1022 "
                        "ORDER BY random() LIMIT 20").fetchall()
    con.close()
    new = emb([(c, s) for c, s, _ in rows])
    cos = []
    for c, s, blob in rows:
        old = np.frombuffer(blob, dtype=np.float32)
        n = new[c]
        cos.append((float(np.dot(old, n) / (np.linalg.norm(old) * np.linalg.norm(n))), len(s), c))
    cos.sort()
    log(f"check vs v5: min cosine {cos[0][0]:.6f} (len {cos[0][1]}), median {cos[len(cos)//2][0]:.6f}")
    return cos[0][0]

def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--missing-from-db"); g.add_argument("--fasta")
    ap.add_argument("--check-db", required=True, help="v2.2 DB with v5 embeddings for the compatibility check")
    ap.add_argument("--out", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--toks", type=int, default=12000); ap.add_argument("--trunc", type=int, default=1022)
    ap.add_argument("--shard", type=int, default=20000); ap.add_argument("--min-cos", type=float, default=0.9999)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    if a.missing_from_db:
        con = sqlite3.connect(f"file:{a.missing_from_db}?mode=ro", uri=True)
        items = con.execute("SELECT s.cdb_id, s.sequence FROM sequences s LEFT JOIN features_esm2 e "
                            "ON e.sequence_id=s.id WHERE e.sequence_id IS NULL ORDER BY s.id").fetchall()
        con.close(); src = f"missing features_esm2 rows in {a.missing_from_db}"
    else:
        items = list(read_fasta(a.fasta)); src = a.fasta
    log(f"{len(items):,} sequences from {src}")

    emb = Embedder(a.toks, a.trunc)
    mc = check_against_db(emb, a.check_db)
    if mc < a.min_cos:
        sys.exit(f"STOP: embeddings differ from v5 (min cosine {mc:.6f} < {a.min_cos}). Check pooling/truncation in v5 04e_esm2.py.")

    shards = [items[i:i + a.shard] for i in range(0, len(items), a.shard)]
    for k, sh in enumerate(shards):
        f = os.path.join(a.out, f"shard_{k:05d}.npy")
        if os.path.exists(f): continue
        t = time.time(); e = emb(sh)
        np.save(f + ".part.npy", np.stack([e[i] for i, _ in sh]).astype(np.float32))
        os.replace(f + ".part.npy", f)
        log(f"shard {k+1}/{len(shards)} ({len(sh):,} seqs, {time.time()-t:.0f} s)")

    log("merging shards")
    arr = np.lib.format.open_memmap(os.path.join(a.out, "esm2.npy"), mode="w+", dtype=np.float32, shape=(len(items), 1280))
    r = 0
    for k in range(len(shards)):
        x = np.load(os.path.join(a.out, f"shard_{k:05d}.npy")); arr[r:r + len(x)] = x; r += len(x)
    arr.flush(); del arr
    with open(os.path.join(a.out, "row_index.tsv"), "w") as f:
        f.write("row\tid\tlength\n")
        for i, (sid, s) in enumerate(items): f.write(f"{i}\t{sid}\t{len(s)}\n")
    h = hashlib.sha256(open(os.path.join(a.out, "esm2.npy"), "rb").read()).hexdigest()
    json.dump({"stage": "04e_esm2_embed", "name": a.name, "model": "esm2_t33_650M_UR50D", "layer": 33,
               "pooling": "mean over residues 1..L (no BOS/EOS)", "truncate": a.trunc, "n": len(items),
               "shape": [len(items), 1280], "dtype": "float32", "sha256": h, "check_min_cosine_vs_v5": mc,
               "source": src, "host": socket.gethostname(), "finished": time.strftime("%Y-%m-%dT%H:%M:%S")},
              open(os.path.join(a.out, "run.json"), "w"), indent=1)
    for k in range(len(shards)): os.remove(os.path.join(a.out, f"shard_{k:05d}.npy"))
    log(f"done: {a.out}/esm2.npy  {len(items):,} x 1280  sha256 {h[:12]}")

if __name__ == "__main__":
    main()
