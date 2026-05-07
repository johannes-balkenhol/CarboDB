# Webapp — backend, frontend, and how they connect

A close-up of the running web application. Read this when you want to
add a new endpoint, extend a Vue component, or debug request flow.

---

## Backend (FastAPI / uvicorn)

Repository: `~/Projects_shared/CarboDB_v3/webapp/`

### Entry point

`webapp/app/main.py` — creates the FastAPI app, sets up CORS, registers
all routers. The `lifespan` hook calls `load_all_models()` from
`startup.py` once at boot, attaching the loaded XGBoost + ESM-2 to a
`ModelStore` singleton.

### Route modules (under `webapp/app/routes/`)

| File | Endpoints |
|---|---|
| `predict.py` | `POST /api/v1/predict` — synchronous live prediction for a single sequence |
| `batch.py` | `POST /api/v1/batch` (multipart FASTA upload) and `GET /api/v1/batch/{job_id}` (polling) |
| `browse.py` | `GET /api/v1/browse`, `GET /api/v1/stats`, `GET /api/v1/db/seq/{uniprot_id}` |
| `external.py` | `GET /api/v1/external/{uniprot_id}` (combined UniProt+AlphaFold metadata), `GET /api/v1/external/{uniprot_id}/structure` (PDB proxy), `GET /api/v1/external/_admin/cache_stats` (diagnostic) |

### Pipeline modules (under `webapp/app/pipeline/`)

These do the actual prediction work, called from the routes.

| File | Purpose |
|---|---|
| `predict.py` | Orchestrator — takes a sequence, runs the 3-model cascade, returns the JSON payload |
| `feature_extraction.py` | Wrapper that calls composition/HMMER/InterProScan/ESM-2 |
| `feature_extraction_v5.py` | Newer v5 feature stack (current production) |
| `feature_extraction_old.py` | Older v3/v4 stack (kept for benchmark replay; safe to remove eventually) |
| `composition.py` | AAC/dipeptide/PAAC/physicochem features — pure numpy |
| `annotate.py` | HMMER and InterProScan subprocess calls |
| `blast_similar.py` | BLAST against per-EC db with experimental Km, returns top-3 nearest hits |
| `carbodb_config.py` | Path constants and feature-name lists |

### Startup module (`webapp/app/startup.py`)

The `ModelStore` class is module-level — loaded once at process start, never
reloaded. Contains:
- `xgb_binary`, `xgb_ec`, `xgb_km`, `xgb_km_final` (XGBoost JSON-loaded)
- `esm_model`, `esm_alphabet`, `esm_device` (ESM-2 from facebookresearch/esm)
- `feature_names_binary`, `feature_names_km` (column-order dictionaries)
- `ec_label_map`, `ec_inv_map` (39-class index ↔ EC string mapping)
- `EC_NAMES` (display names like "Ribulose bisphosphate carboxylase (RuBisCO)")

If you ever see "ModelStore.ready is False" in logs, the lifespan hook
didn't complete — check the trace; usually a missing model file under
`webapp/models/` or `data/models/`.

### How a /predict request flows

1. POST hits `predict.py:predict_endpoint`
2. Sequence validated (length, alphabet)
3. `pipeline/predict.py:predict_sequence(sequence, kingdom, mode)` called
4. `feature_extraction_v5.extract_features()`:
   - composition (~ms)
   - HMMER subprocess (5–15 s)
   - InterProScan subprocess (30–60 s) — the bottleneck
   - ESM-2 forward pass (20–40 s on CPU, 1 s on GPU)
5. Vector assembled into a 1793-feature numpy array
6. `xgb_binary.predict_proba` — gives `is_carboxylase` boolean + probability
7. If positive: `xgb_ec.predict_proba` — top class + confidence
8. If CO₂-active EC: `xgb_km.predict` — log10(Km in mM)
9. SHAP from precomputed per-EC contributions
10. `blast_similar.find_nearest_with_km()` — top-3 BLAST hits with experimental Km
11. Result assembled into ResultDetail-shaped JSON, returned

Total: ~100 s on wbbi206 CPU for a 470-aa protein. Returned as a single
synchronous response — connection held open the whole time.

### How a /batch request flows

1. POST with multipart file → `batch.py` writes FASTA to
   `webapp/jobs/{job_id}/in.fasta`
2. Spawns a background asyncio task that iterates sequences and calls
   `predict_sequence()` per sequence, writing each result as JSON to
   `webapp/jobs/{job_id}/results/{seq_idx}.json`
3. Frontend polls `/api/v1/batch/{job_id}` for progress; backend returns
   `{progress: "3/12", done: false, results_so_far: [...]}` until complete

Per-sequence cache files in `webapp/jobs/<job_id>/results/` mean that
clicking "Details" on a finished batch row doesn't re-run the model — it
just reads the cached JSON.

### How a /browse + /db/seq + /external trio renders the Database tab

1. User opens Database tab → frontend fetches `/api/v1/stats` and
   `/api/v1/browse` in parallel.
2. `/stats` returns banner counts (cached 10 min in process memory after
   first call).
3. `/browse?limit=50` returns first page of rows. User filters/sorts →
   another `/browse` call with new params.
4. User clicks a Details button → frontend fetches `/api/v1/db/seq/{uid}` →
   gets ResultDetail-shaped JSON → opens modal with `<ResultDetail>`.
5. User clicks "🔬 Show extended annotation" inside the modal → frontend
   fetches `/api/v1/external/{uid}` → renders ExtendedDetails component
   with function text, GO terms, taxonomy, cross-refs.
6. ExtendedDetails mounts StructureViewer → fetches
   `/api/v1/external/{uid}/structure` (returns AlphaFold PDB) → NGL
   Viewer renders 3D model with pLDDT coloring.

Steps 1–3 take milliseconds (warm). Step 4 takes ~0.5 s. Step 5 takes
~1 s cold, ~80 ms cached. Step 6 takes ~1 s plus NGL render time.

### SQLite access pattern

`browse.py:_db()` opens a fresh connection per request. Setup:

```python
conn = sqlite3.connect(path, timeout=30)
conn.row_factory = sqlite3.Row
conn.execute("PRAGMA cache_size = -524288")     # 512 MB page cache
conn.execute("PRAGMA temp_store = MEMORY")
conn.execute("PRAGMA mmap_size = 30000000000")  # 30 GB mmap
```

These PRAGMAs are critical. Without them /browse takes 17 s; with them ~12 s
(the remaining cost is real query work). The full-scan join over
`predictions` is what dominates — see ROADMAP for the rewrite plan.

The `external_annotations_cache` table is created lazily by `external.py`
on first import. Schema:
```sql
CREATE TABLE external_annotations_cache (
    uniprot_id TEXT PRIMARY KEY,
    payload    TEXT NOT NULL,        -- JSON
    fetched_at INTEGER NOT NULL,     -- unix ts
    status     TEXT NOT NULL         -- 'ok' | 'partial'
);
```

---

## Frontend (Vue 3 / Vite)

Repository: `~/Projects_shared/CarboDB-App-v2/frontend/`

### Entry and routing

- `src/main.js` — Vue app boot, creates router + Pinia, mounts `App.vue`
- `src/App.vue` — global layout (header + `<RouterView>`)
- `src/router/index.js` — 5 routes (see ARCHITECTURE.md)

### Views (`src/views/`)

| View | Route | Purpose |
|---|---|---|
| `HomeView.vue` | `/` | Landing page, project intro |
| `AnalysisView.vue` | `/analysis` | Runs /predict on user input or /batch on uploaded FASTA |
| `DatabaseView.vue` | `/database` | Browse precomputed predictions with filters and pagination |
| `AllSearchesView.vue` | `/searches` | History of past Analysis runs (Pinia store) |
| `AboutView.vue` | `/about` | Static about |
| `HmmerSearchView.vue` | (not in nav) | Legacy, mostly unused |

### Components (`src/components/`)

The component tree for the Database flow:

```
DatabaseView
└── (table of /browse results, each row has Details button)
    └── CommonModal (opens on Details click)
        └── ResultDetail (props: result)
            ├── (verdict, Km comparison, Pfam, motifs, physicochem all inline)
            ├── ResultDetailItem  (one of several sub-blocks)
            └── ExtendedDetails (props: uniprot-id, api-base, pfam-hits)
                └── StructureViewer (props: uniprot-id, api-base, pfam-hits)
                     └── (NGL Viewer canvas, color dropdown, reset, download)
```

The same `ResultDetail` component is also used by AnalysisView, with the
same `result` prop shape. The only behavioral difference is `result.mode`:
- `mode === 'live'` — from /predict, has SHAP + 3 BLAST hits
- `mode === 'db_lookup'` — from /db/seq, no SHAP, single self-hit
- `mode === 'batch'` — from /batch's per-row, has SHAP

### How `<ExtendedDetails>` is wired in

Currently in `ResultDetail.vue`, before the closing `</template>`:

```vue
<ExtendedDetails
  v-if="result && result.mode === 'db_lookup' && (result.uniprot_id || result.id)"
  :uniprot-id="result.uniprot_id || result.id"
  :api-base="'/api/v1'"
  :pfam-hits="result.pfam_hits || []"
/>
```

The v-if guard is **intentional**:
- For `mode === 'live'`, the user's submitted sequence may not be in
  UniProt at all, so /external/{id} would 404. The frontend would still
  show the button, the user would click, and they'd get an error.
- For `mode === 'db_lookup'`, every sequence has a real UniProt ID
  (the database is *built from* UniProt), so /external/{id} always
  has something to return.

**Future extension** (see ROADMAP): for `mode === 'live'`, we could show
ExtendedDetails for the **nearest BLAST hit** instead — that hit always
has a UniProt ID, and is highly relevant context for the user's sequence.
This is a one-prop change.

### Pinia stores

`src/stores/searches.js` — keeps the user's recent Analysis history in
localStorage. Each entry has `{id, timestamp, sequence_preview, result}`.
Used by AllSearchesView.

`src/stores/validation.js` — FASTA validation helpers (regex checks for
allowed amino acid letters, FASTA header parsing, multi-vs-single
auto-detect). Used by FileUpload.

### NGL Viewer + Vue 3 reactivity gotcha

`StructureViewer.vue` uses NGL Viewer (loaded via CDN script tag in
`index.html`). NGL is built on Three.js, which uses Object3D properties
that are read-only by JavaScript spec.

**If you ever assign an NGL `Stage`, `Component`, or `Representation`
object to a `ref()` or `reactive()` field, Vue wraps it in a Proxy and
NGL/Three.js will throw 1000+ errors per second during render** — we hit
this on May 6 and recovered. Pattern that works:

```js
import { markRaw } from 'vue'

let stage = null   // plain let, NOT a ref
const canvasEl = ref(null)

onMounted(() => {
  stage = new NGL.Stage(canvasEl.value)
  // OR if you need it to be reactive for some reason:
  // myReactive.stage = markRaw(new NGL.Stage(...))
})
```

This applies to anything from NGL or Three. Keep them out of Vue's
reactivity graph.

### Vite dev proxy

Vite serves the frontend on :5173. Requests to `/api/*` are proxied to
:8090 (the FastAPI backend). Proxy config is in
`frontend/vite.config.js`. In production this would be replaced by an
nginx or Caddy rule (see DEPLOYMENT.md).

### Frontend → backend contract

Every `result` object the frontend renders has this shape (abridged):

```typescript
{
  id: string,                    // UniProt ID for db_lookup, FASTA header for live
  cdb_query_id?: string,         // CDB internal ID (db_lookup only)
  uniprot_id?: string,           // explicit UniProt ID (db_lookup only, May 7+)
  mode: "live" | "db_lookup" | "batch",
  
  is_carboxylase: boolean,
  carboxylase_probability: number,    // 0-1
  ec_predicted: string,               // "4.1.1.39"
  ec_name: string,                    // "Ribulose bisphosphate carboxylase (RuBisCO)"
  ec_confidence: number,
  ec_probabilities?: {[ec: string]: number},  // top-N
  
  km_predicted_mM: number,
  km_predicted_uM: number,
  km_predicted_log10: number,
  
  sequence_length: number,
  pfam_hits: Array<{accession: string, name?: string, evalue?: string, bitscore?: number}>,
  features_computed: {[name: string]: number},   // motif counts, physicochem, etc.
  shap?: {[ec: string]: Array<[name: string, value: number]>},
  top_similar: Array<{
    uniprot_id: string, organism: string, identity_pct: number,
    evalue: number, align_length: number,
    km_experimental_uM?: number, km_exp_substrate?: string, tier: string
  }>,
  
  organism?: string,
  reviewed?: 0 | 1,
  source?: "brenda" | "swissprot" | "trembl",
  novelty_flag?: "known" | "novel" | "very-novel",
  kingdom?: string,
  features_used?: string[],
  runtime_seconds?: number,
  warnings?: string[]
}
```

Adding a new field: add it to the backend response, then to the
frontend rendering. The shape lives implicitly across both repos —
when this gets bigger, make a shared TypeScript type definition in
the frontend.

---

## How to debug a stuck request

**Backend stuck:**
```bash
ps aux | grep uvicorn
tail -f ~/Projects_shared/CarboDB_v3/webapp/logs/webapp.log
```

If uvicorn is alive but a request hangs, it's almost always inside HMMER
or InterProScan. Check `ps aux | grep -E 'hmmscan|interproscan'`. Stuck
java processes are common; kill them and uvicorn will return a 504/timeout.

**Frontend not reloading:**
- vite HMR sometimes misses additions of new component files. Hard refresh
  the browser (Ctrl+Shift+R).
- If still wrong, restart vite: `~/Projects_shared/CarboDB-App-v2/start_app.sh restart`

**Wrong field name in /db/seq response:**
- We had a `cdb_query_id` vs `uniprot_id` confusion that bit us May 7.
  `cdb_query_id` looks like a UniProt ID variable but holds the CDB
  internal ID. Check `docs/ROADMAP.md` for the rename plan.

**3D structure won't render:**
- Open browser DevTools (F12) → Network tab → look for the
  `/api/v1/external/{uid}/structure` request. Should be 200 with
  ~300 KB of `chemical/x-pdb`. If 404, see external.py for AlphaFold
  version fallback list (currently `[6, 5, 4]`).
- If 200 but blank canvas, see the Vue 3 reactivity gotcha above.

**Browse is slow:**
- Cold first call after backend restart is slow because cache is empty.
  Second call should be fast. If it stays slow, check:
  `sqlite3 carbodb.sqlite 'PRAGMA cache_size'` should return `-524288`
  (negative = KB).

**External annotation 404:**
- The validator in `external.py:_validate_uniprot_id()` checks
  `sequences.uniprot_id` before any upstream call. If the ID isn't in
  the database, you get a 404 with a clear message. This is by design —
  prevents abuse of the proxy as a generic UniProt scraper.
