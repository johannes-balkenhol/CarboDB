# Deployment

How to install CarboDB on a fresh server. Two stages:

1. **Stage 1 — Single-host, single-process** (what wbbi206 currently runs).
   Public-facing but with no concurrency tolerance. Good for ≤10 concurrent
   users, scientific demos, internal-team-only deployments.

2. **Stage 2 — Multi-process with job queue** (planned, not yet built).
   Required for true public release where concurrent /predict requests
   shouldn't queue against each other.

---

## Stage 1 — what's currently deployed

### Host requirements

| Resource | Minimum | Recommended | Why |
|---|---|---|---|
| OS | Ubuntu 22.04 LTS | same | what wbbi206 runs |
| RAM | 32 GB | 64 GB+ | ESM-2 inference + SQLite mmap |
| Disk | 350 GB | 500 GB | data is 285 GB; need headroom for indices, intermediates |
| CPU | 16 cores | 64+ cores | parallel HMMER + InterProScan |
| GPU | none required | any 24 GB+ NVIDIA | makes ESM-2 ~30× faster |
| Network | 100 Mbps | 1 Gbps | UniProt + AlphaFold pulls |

### Step 1 — clone both repos

```bash
mkdir -p ~/Projects_shared && cd ~/Projects_shared
git clone git@github.com:johannes-balkenhol/CarboDB.git CarboDB_v3
git clone git@github.com:johannes-balkenhol/CarboDB-App.git CarboDB-App-v2
```

(The `_v3` and `-v2` suffixes are historical; we kept them after a
project-naming evolution. Don't change them — paths are hardcoded in
several places, including `start_app.sh`.)

### Step 2 — install conda environment

```bash
cd ~/Projects_shared/CarboDB_v3
# Install miniforge if not already
# Then:
conda env create -f environment.yml
conda activate carboxylase
```

The `environment.yml` is large (~150 packages) because it pins exact
versions of bioinformatics tools (CD-HIT 4.8.1, HMMER, BLAST,
ESM-2 dependencies). Don't bump these without testing.

Required additional pip packages (some are not in the conda env):
```bash
pip install httpx>=0.25 fastapi uvicorn[standard] python-multipart
```

### Step 3 — install bioinformatics tools

The conda env should provide most. Verify:
```bash
which hmmscan       # HMMER
which interproscan.sh
which cd-hit
which blastp makeblastdb
```

If `interproscan.sh` is missing, install separately from
https://www.ebi.ac.uk/interpro/download/InterProScan/ (it's a 12 GB Java
distribution that doesn't fit nicely in conda).

### Step 4 — get the data

The 50 GB SQLite is too large for git. Options:

**Option A: rsync from existing wbbi206 deployment.**
```bash
mkdir -p data/primary data/models data/dbs
rsync -avh --progress wbbi206:~/Projects_shared/CarboDB_v3/data/primary/carbodb.sqlite data/primary/
rsync -avh --progress wbbi206:~/Projects_shared/CarboDB_v3/data/models/ data/models/
rsync -avh --progress wbbi206:~/Projects_shared/CarboDB_v3/data/dbs/pfam/ data/dbs/pfam/
```

**Option B: rebuild from scratch.** Run scripts 01–10 in order. Takes
~5 days end-to-end. See `docs/DATA_INGESTION.md`.

The webapp also needs `webapp/models/`:
```bash
rsync -avh --progress wbbi206:~/Projects_shared/CarboDB_v3/webapp/models/ webapp/models/
```

For BLAST nearest-neighbor search:
```bash
rsync -avh --progress wbbi206:~/Projects_shared/CarboDB_v3/webapp/data/blast_dbs/ webapp/data/blast_dbs/
```

(These are pre-built per-EC BLAST databases. Rebuilding takes ~6 hours via
`webapp/scripts/build_ec_blast_dbs.py`.)

### Step 5 — install frontend dependencies

```bash
cd ~/Projects_shared/CarboDB-App-v2/frontend
npm install
```

This pulls Vue 3, Vite, Pinia, lucide icons, and a few more (~190 MB).

### Step 6 — start the app

```bash
cd ~/Projects_shared/CarboDB-App-v2
./start_app.sh restart
```

Open http://YOUR_HOST:5173 in a browser. Should see the home page within
~10 seconds (uvicorn lifespan loads ESM-2; vite first-render is fast).

### Step 7 — verify

```bash
# Backend health
curl -s http://localhost:8090/api/v1/stats | python3 -m json.tool | head -20

# DB lookup of a known sequence
curl -s http://localhost:8090/api/v1/db/seq/P00875 | python3 -m json.tool | head -20

# External proxy
curl -sI http://localhost:8090/api/v1/external/P00875/structure | head -8
# Expected: HTTP/1.1 200, content-type: chemical/x-pdb
```

### Stage 1 limitations

- **uvicorn is single-worker, single-process.** A second simultaneous
  /predict request waits in the asyncio queue while the first holds the
  ESM-2 forward pass. Effective concurrency ≈ 1 for /predict.
- **No reverse proxy / TLS.** Port 5173 is plain HTTP. Browsers complain
  about mixed-content if any page links out.
- **No daemon.** `start_app.sh` uses `nohup`; if the shell session is
  killed, processes survive but aren't supervised. A reboot loses state.
- **No backup.** SQLite is on local disk only.

These are all fine for the current "internal-and-collaborator demo"
deployment. They're not fine for a public-facing site with traffic.

---

## Stage 2 — plan for production

Not yet implemented. This is the **rough plan to discuss with Johannes**
before any of it is built.

### Concurrency: jobs queue

```
                    ┌──────────────────────────────────┐
  HTTP /predict ──► │ FastAPI worker (4×, behind nginx)│
                    │   - validates input              │
                    │   - enqueues job to Redis        │
                    │   - returns job_id immediately   │
                    └─────────────┬────────────────────┘
                                  ▼
                    ┌──────────────────────────────────┐
                    │ Redis queue (RQ or Celery)       │
                    └─────────────┬────────────────────┘
                                  ▼
                    ┌──────────────────────────────────┐
                    │ Worker process (1× per GPU/CPU)  │
                    │   - dequeues, runs prediction    │
                    │   - writes result to Redis       │
                    │   - publishes "done" event       │
                    └──────────────────────────────────┘
                                  ▲
                                  │
  HTTP /predict/{job_id} ─────────┘  (poll or SSE for updates)
```

This converts /predict from synchronous to async. Frontend either polls
the job_id endpoint or listens via Server-Sent Events for completion.

**Why Redis + RQ specifically:** simpler than Celery, fewer moving parts,
worker-restart-safe, supports priority queues out of the box. Celery is
acceptable but heavier than we need.

### Reverse proxy + TLS

nginx or Caddy in front:

```nginx
server {
    listen 443 ssl http2;
    server_name carbodb.example.org;
    ssl_certificate     /etc/letsencrypt/live/.../fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/.../privkey.pem;

    # Static frontend (built bundle)
    root /var/www/carbodb-frontend;
    try_files $uri $uri/ /index.html;

    # API proxy
    location /api/ {
        proxy_pass http://127.0.0.1:8090;
        proxy_set_header X-Forwarded-For $remote_addr;
        proxy_read_timeout 180s;        # /predict can take 100+ s
    }
}
```

Build the frontend with `npm run build` and serve the `dist/` directory.

### Process supervision: systemd units

```ini
# /etc/systemd/system/carbodb-api.service
[Unit]
Description=CarboDB FastAPI backend
After=network.target

[Service]
Type=simple
User=carbodb
WorkingDirectory=/srv/CarboDB_v3
Environment=DB_PATH=data/primary/carbodb.sqlite
Environment=PFAM_HMM=data/dbs/pfam/Pfam-A.hmm
Environment=MODELS_DIR=webapp/models
Environment=ESM2_DEVICE=cpu
ExecStart=/srv/conda/envs/carboxylase/bin/uvicorn webapp.app.main:app \
          --host 127.0.0.1 --port 8090 --workers 4
Restart=on-failure
RestartSec=5

[Install]
WantedBy=multi-user.target
```

Per-worker memory: ESM-2 is ~3 GB; running 4 workers means 12 GB
just for that. Either accept it (with 64 GB RAM machine) or shift
ESM-2 to a separate worker pool that the FastAPI workers query over IPC.

### Database: SQLite is fine, until it's not

50 GB SQLite handles single-user well. With 4 concurrent uvicorn workers
sharing the same file via `mmap_size`, we should be fine for read traffic.
**Writes** would be a problem (SQLite serializes writers) — but this app
is read-only at runtime (the only writes are to `external_annotations_cache`,
which is low-volume and tolerant to occasional contention).

If write contention shows up: switch to WAL mode (`PRAGMA journal_mode = WAL`)
or migrate the cache to Redis.

### Data migration: when models retrain

Retraining pipeline (manual today):
1. Re-run scripts 06–10 on a staging copy of the DB
2. New model artifacts land in `data/models/<v6>.json`
3. New `webapp/models/feature_names_*.json` mirrors written
4. `predictions` table re-populated by script 10 (∼12 hours on full DB)
5. SQLite swap: copy new `carbodb.sqlite` to production, atomic rename

The webapp loads models at startup. To deploy a new model:
- Stop systemd unit
- Replace files
- Start systemd unit

There's no model-hot-swap path. Acceptable for v1 but should be considered
for v6+ if retraining cadence increases.

### Monitoring

Plan, in order of priority:
1. **uptime check** — ping `/api/v1/stats` every minute, alert on >2 failures
2. **error rate** — count 5xx in nginx logs, alert at >1% over 5 min
3. **job queue depth** — Redis LLEN, alert if >50 jobs queued
4. **disk usage** — alert at <10 GB free (the `external_annotations_cache`
   could grow indefinitely if many users; 30-day TTL helps but isn't enforced)

Tooling: Prometheus + Grafana is overkill for now. A 50-line bash script
running from cron, mailing on failure, is fine for v1.

### Hardware split

Right now everything runs on `wbbi206`. The MOTD warns this is for
interactive use only. For real production:

- **API host** (4–8 cores, 32 GB RAM, 100 GB disk) — uvicorn workers + nginx
- **Compute host** (32+ cores, 64+ GB RAM, GPU optional) — RQ workers
  with HMMER, InterProScan, ESM-2
- **Storage host** (or shared NFS/SAN) — the 50 GB SQLite + Pfam HMMs +
  BLAST DBs + InterPro DBs

Or collapse all three onto a single beefy machine if the budget allows
(prefer 1 big machine over 3 small ones for SQLite-based architectures —
mmap across hosts is painful).

### What to ask Johannes before building Stage 2

1. Expected concurrent users — 5? 50? 500? This drives worker count and
   queue capacity.
2. Acceptable response time for /predict — synchronous would mean ≥100 s
   wait; queued means user gets a job_id immediately and polls. Both are
   defensible UX choices.
3. Whether GPU is available. ESM-2 on GPU collapses a 30 s step to 1 s;
   transforms the response budget entirely.
4. Whether we need SSL/TLS termination or if it lives behind a campus VPN.
5. Backup strategy. SQLite snapshots? S3 mirror?
6. Logging compliance. Are we required to log requests? Is there a privacy
   policy for the sequences users submit (some labs care a lot about this)?

---

## Backup & disaster recovery

### What to back up

- `data/primary/carbodb.sqlite` — irreplaceable without ~5 days of pipeline
- `data/models/*.json` — irreplaceable without retraining
- `webapp/models/*.json` — derived from above, but small and fast
- Both git repos — but those are on GitHub, lower priority
- `data/intermediate/*.tsv` — useful but rebuildable

### What to skip

- `data/features/` (intermediate per-feature outputs; rebuildable)
- `webapp/jobs/` (transient batch artifacts)
- `tmp/`, `__pycache__/`, `node_modules/`

### Backup procedure (suggested)

```bash
# Daily SQLite snapshot
sqlite3 data/primary/carbodb.sqlite ".backup '/backup/carbodb-$(date +%F).sqlite'"

# Weekly model + webapp config
rsync -av data/models/ webapp/models/ webapp/jobs/ /backup/weekly/
```

### Restore procedure

1. Stop systemd / start_app.sh
2. Copy SQLite back to `data/primary/carbodb.sqlite`
3. Copy models back to `data/models/` and `webapp/models/`
4. Restart
5. Verify with `curl /api/v1/stats`

---

## Known operational gotchas

- **First /stats call after backend restart is slow (40 s)** because
  the cache is cold and the query scans `predictions`. Subsequent calls
  are 15 ms. Document this for users; do not let monitoring alert on
  the cold call.

- **Vite HMR sometimes drops new files.** If a Vue component is *added*
  (not just edited), vite occasionally fails to register it without a
  hard restart. `start_app.sh restart` always works.

- **InterProScan can hang on certain sequences** (rare but documented,
  see DATA_INGESTION.md). The 60-min batch timeout protects the pipeline
  but a single stuck request will hold a uvicorn worker. With multiple
  workers (Stage 2) this is less critical.

- **AlphaFold proxy depends on AlphaFold's URL pattern.** They bumped from
  v4 → v6 between training and now. The proxy walks `[6, 5, 4]` versions;
  add new versions to `external.py:ALPHAFOLD_VERSIONS_TO_TRY` if EBI bumps
  again.

- **`/tmp` cleanup** — wbbi206's `/tmp` is volatile. Don't store anything
  there that needs to survive a reboot. Batch jobs should write to
  `webapp/jobs/` not `/tmp/`.
