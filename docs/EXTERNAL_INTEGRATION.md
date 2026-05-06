# External Annotation Integration — Setup Guide

This adds a lazy-loaded "Extended annotation" panel to the sequence Details
view, fetching UniProt + AlphaFold data on demand and rendering an interactive
3D structure with NGL Viewer.

**Files in this package:**
- `backend/external.py`         — FastAPI router, drop into `webapp/app/routes/`
- `frontend/StructureViewer.vue` — NGL Viewer Vue component
- `frontend/ExtendedDetails.vue` — Lazy-loaded extended details panel
- `INTEGRATION.md`              — this file

**Estimated integration time:** 30-45 minutes for the first run.


## Backend integration (3 steps)

### Step 1: install httpx if missing

```bash
cd ~/Projects_shared/CarboDB_v3
pip install httpx
# or, to track it:
echo "httpx>=0.25" >> webapp/requirements.txt
```

### Step 2: drop the route file into place

```bash
cp /path/to/external.py webapp/app/routes/external.py
```

### Step 3: register the router in webapp/app/main.py

Open `webapp/app/main.py` and find where the other routers are registered.
Look for a block that probably looks like:

```python
from .routes import predict, batch, browse
app.include_router(browse.router,  prefix="/api/v1")
app.include_router(predict.router, prefix="/api/v1")
app.include_router(batch.router,   prefix="/api/v1")
```

Add the external router:

```python
from .routes import predict, batch, browse, external
app.include_router(browse.router,   prefix="/api/v1")
app.include_router(predict.router,  prefix="/api/v1")
app.include_router(batch.router,    prefix="/api/v1")
app.include_router(external.router, prefix="/api/v1")
```

The router itself already has `prefix="/external"` defined, so the final
endpoint path is `/api/v1/external/{uniprot_id}` — matching what the Vue
component expects.

### Step 4: restart uvicorn and verify

```bash
# Kill old uvicorn (if running)
pkill -f "uvicorn webapp.app.main"
sleep 2

# Restart via your usual mechanism (start_app.sh, nohup, etc.)
~/Projects_shared/CarboDB-App-v2/start_app.sh restart

# Test that the cache table got created
sqlite3 ~/Projects_shared/CarboDB_v3/data/primary/carbodb.sqlite \
  "SELECT name FROM sqlite_master WHERE name='external_annotations_cache';"
# Should print: external_annotations_cache

# Test fetch (should take ~1-2s on first call, hit UniProt + AlphaFold)
time curl -s 'http://localhost:8090/api/v1/external/P00875' | head -50

# Second call should be instant from cache
time curl -s 'http://localhost:8090/api/v1/external/P00875' | head -50

# Cache stats
curl -s 'http://localhost:8090/api/v1/external/_admin/cache_stats'
# {"total":1,"ok_count":1,"partial_count":0,"oldest":1715000000,"newest":1715000000}

# Test invalid uniprot_id (should 404 without hitting upstream)
curl -i 'http://localhost:8090/api/v1/external/AAAAAA'
# HTTP/1.1 404 Not Found  -- "uniprot_id AAAAAA not in CarboDB"

# Test rate limit (run this 65 times in a loop, last few should 429)
for i in $(seq 1 65); do
  curl -s -o /dev/null -w "%{http_code}\n" \
    'http://localhost:8090/api/v1/external/P00875';
done | sort | uniq -c
# Should see ~60 200s and ~5 429s
```

If all four tests pass, the backend is good. Move on to frontend.


## Frontend integration (4 steps)

### Step 1: load NGL Viewer via CDN

NGL is loaded as a global, not bundled, to avoid Vite/webpack complications.
Edit `~/Projects_shared/CarboDB-App-v2/frontend/index.html` and add inside `<head>`:

```html
<script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.39/dist/ngl.js"></script>
```

This is a ~250 KB file, loaded once per page. NGL exposes itself as the
global `NGL` namespace, which the StructureViewer component checks for.

### Step 2: drop the Vue components into place

```bash
cp /path/to/StructureViewer.vue \
   ~/Projects_shared/CarboDB-App-v2/frontend/src/components/

cp /path/to/ExtendedDetails.vue \
   ~/Projects_shared/CarboDB-App-v2/frontend/src/components/
```

### Step 3: wire ExtendedDetails into ResultDetail.vue

Open `~/Projects_shared/CarboDB-App-v2/frontend/src/components/ResultDetail.vue`
(this is the component that renders inside the modal in `DatabaseView.vue`).

At the top of the `<script>` block, add the import:

```javascript
import ExtendedDetails from "./ExtendedDetails.vue";
```

Add it to the components registration:

```javascript
export default {
  name: "ResultDetail",
  components: {
    // ... existing components ...
    ExtendedDetails,
  },
  // ...
}
```

Then near the bottom of the `<template>` (after the SHAP section, after the
nearest-neighbors / BLAST section, before whatever closing `</div>` wraps the
panel), insert:

```html
<!-- Extended UniProt + AlphaFold annotation, lazy-loaded -->
<ExtendedDetails
  v-if="result && result.uniprot_id"
  :uniprot-id="result.uniprot_id"
  :api-base="apiBase || '/api/v1'"
  :pfam-hits="result.pfam_hits || []"
/>
```

The `v-if` guard means the section is hidden for novel sequences submitted
via /predict (they have no UniProt ID). It only shows for stored DB sequences.

### Step 4: restart vite and test

```bash
~/Projects_shared/CarboDB-App-v2/start_app.sh restart
```

Then in the browser:
1. Open the database browser, click any RuBisCO row (P00875 spinach is good)
2. Scroll to the bottom of the Details panel
3. You should see a button: "🔬 Show extended annotation"
4. Click it. After 1-2 seconds, the panel populates with function text,
   GO terms, active sites, taxonomic lineage, cross-references, and the
   AlphaFold structure colored by pLDDT.
5. Try the "Color by" dropdown — switch between pLDDT, Pfam, motifs, rainbow.

If the structure shows up but coloring by Pfam doesn't work, see "Pfam
positions" below.


## Troubleshooting

### Backend: `pip install httpx` fails

If you can't install on the HPC due to package permissions:

```bash
pip install --user httpx
# or in conda env:
conda install -c conda-forge httpx
```

### Backend: cache table not created

If `external_annotations_cache` doesn't exist after restart, the import of
`external.py` failed silently. Check uvicorn logs:

```bash
tail -100 ~/Projects_shared/CarboDB_v3/webapp/logs/webapp.log
```

Most common cause: the import path in `main.py` is wrong. The router file
must be importable as `webapp.app.routes.external`.

### Frontend: "NGL Viewer library not loaded"

The CDN script tag in `index.html` either isn't there or hasn't been
re-loaded. Hard-refresh the browser (Ctrl+Shift+R). Check the browser console
for failed network requests to `cdn.jsdelivr.net`.

### Frontend: structure viewer is blank/black

NGL needs the canvas's parent element to have a non-zero height. The
component sets `height: 420px` in scoped CSS but if your parent CSS
overrides this, the canvas collapses. Inspect the `.sv-canvas` element in
DevTools and confirm it has a height.

### Pfam positions not available

CarboDB's stored `pfam_hits` (in `features_domains` table) may not include
`start`/`end` residue positions, depending on whether they were captured
during the HMMER run (`scripts/04b_hmmer.sh`). The ExtendedDetails component
falls back to UniProt's domain features when that's the case, but if neither
source has positions, the "Pfam" coloring option in the structure viewer
will be disabled. You can:

- Fix forward: re-run `04b_hmmer.sh` capturing `--domtblout` and store
  start/end columns in the DB
- Live with it: pLDDT and rainbow coloring still work fine

### "Rate limit exceeded" on legitimate browsing

The 60 requests/minute/IP limit is generous for normal user browsing
(roughly: opening 60 different sequences in a minute, or refreshing
the same sequence 60 times). If you find yourself hitting it during
normal use, raise `RATE_LIMIT_REQUESTS` in `external.py`:

```python
RATE_LIMIT_REQUESTS = 200  # was 60
```

For public deployment, leave it at 60 — it's primarily an anti-abuse measure.

### UniProt or AlphaFold returns garbage / changes their API

Both services have stable v1+ APIs and have not broken backwards compatibility
in years, but if it happens:

- The UniProt JSON parsing is concentrated in `_summarize_uniprot()` —
  edit there if fields move
- The AlphaFold metadata fetch is in `_fetch_alphafold_meta()` — same
- Fallback behaviour: if the parser breaks but the fetch succeeds, you'll
  get an empty UniProt summary. The cache will fill up with empties.
  In that case, manually clear:

```sql
DELETE FROM external_annotations_cache;
```


## Production considerations (when you go public)

The current code is ready for an internal/HPC deployment. For public-facing
production, consider:

1. **CDN for NGL.** Loading from jsdelivr is fine for low-traffic. For
   high traffic, host the NGL JS yourself or use a paid CDN.

2. **Cache size.** The cache grows by ~30 KB per unique sequence. At
   500k sequences that's 15 GB inside `carbodb.sqlite`. Consider moving
   the cache to a separate `external_cache.sqlite` if you're concerned
   about main DB bloat. Change `DB_PATH` in `external.py` to a separate
   file path. The cache table creation is idempotent so it'll work
   automatically.

3. **Rate limit by API key, not IP.** The current sliding-window limiter
   uses client IP, which can be defeated by anyone with multiple IPs.
   For a real public deployment add an API key system. Out of scope
   for v1.

4. **Reverse proxy headers.** When deployed behind nginx/Caddy/Cloudflare,
   make sure the proxy is configured to pass `X-Forwarded-For` so that
   `_get_client_ip()` returns the real user IP. Otherwise everything will
   look like one IP and rate limiting becomes a global cap.

5. **Background cache prewarming.** If you want the Details panel to feel
   instant for popular sequences, write a one-shot script that pre-fetches
   the top N most-clicked uniprot_ids. Out of scope for v1 — the lazy
   approach is fine and minimizes upstream API load.


## What's NOT in this package

These were considered and intentionally deferred:

- **InterPro REST integration.** Adds richer domain info but is a separate
  upstream call. The current Pfam data from `features_domains` is good
  enough to start.
- **STRING DB protein-protein interactions.** Cool but rarely actionable.
- **PDBe / RCSB experimental structures.** Most carboxylases of interest
  have AlphaFold; the few with experimental structures can be reached
  via the cross-references row.
- **Multi-fragment AlphaFold support.** Sequences >2700 aa get split by
  AlphaFold into F1, F2, etc. The current code only fetches F1. Affects
  pyruvate carboxylase and biotin carboxylase fully assembled forms,
  but their catalytic domains (which are what motif analysis cares about)
  are usually in F1.
- **SHAP-colored structures.** Coloring residues by per-position SHAP
  contribution would be the most scientifically interesting view. Requires
  per-residue SHAP, which CarboDB doesn't compute today (it computes
  feature-level SHAP, not residue-level).


## Verifying everything works end-to-end

Quick smoke test once both backend and frontend are deployed:

```bash
# 1. Backend route registered
curl -s http://localhost:8090/api/v1/external/_admin/cache_stats
# {"total": 0, ...}

# 2. Cache populates on first fetch
curl -s http://localhost:8090/api/v1/external/P00875 | jq '.uniprot.protein_name'
# "Ribulose bisphosphate carboxylase large chain"

# 3. PDB proxy returns binary
curl -sI http://localhost:8090/api/v1/external/P00875/structure
# HTTP/1.1 200 OK
# content-type: chemical/x-pdb

# 4. Database browser at http://132.187.22.206:5173/
# - Click any RuBisCO row
# - Scroll to "Show extended annotation" button
# - Click — structure renders, GO terms appear, etc.
```

If all four pass, the integration is complete. Commit and ship.
