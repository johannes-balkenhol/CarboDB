# Architecture

How CarboDB is shaped, top to bottom. Read this after `HANDOFF.md`.

---

## System map

```
                                                  ┌────────────────────────────┐
  ┌─────────────────┐   HTTP/JSON                 │ External APIs              │
  │ Browser (user)  │ ←───────────────────────→   │  • UniProt REST            │
  │ Vue 3 + Vite    │                             │  • AlphaFold DB            │
  │ port 5173       │                             │  • BRENDA SOAP (ingest)    │
  └────────┬────────┘                             └──────────────┬─────────────┘
           │ Vite dev proxy / nginx (production)                 │ httpx
           ▼                                                     ▼
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ FastAPI (uvicorn) — port 8090 — webapp/app/main.py                       │
  │                                                                          │
  │   /api/v1/predict       single seq → JSON (synchronous, ~100 s)          │
  │   /api/v1/batch         multipart FASTA → job ID → poll                  │
  │   /api/v1/browse        DB search/filter/paginate                        │
  │   /api/v1/stats         banner aggregates (cached 10 min)                │
  │   /api/v1/db/seq/{uid}  precomputed prediction in /predict shape         │
  │   /api/v1/external/{uid}    UniProt+AlphaFold metadata (cached 30 days)  │
  │   /api/v1/external/{uid}/structure  AlphaFold PDB proxy                  │
  └────────┬────────────────────────────────────┬────────────────────────────┘
           │                                    │
           ▼                                    ▼
  ┌──────────────────────────┐       ┌──────────────────────────────────────┐
  │ Pipeline (in-process)    │       │ data/primary/carbodb.sqlite (50 GB)  │
  │ webapp/app/pipeline/     │       │                                      │
  │  • feature_extraction.py │       │  sequences           2.38 M rows     │
  │  • predict.py            │       │  predictions         2.38 M rows     │
  │  • blast_similar.py      │       │  km_evidence         2.97 K rows     │
  │  • annotate.py           │       │  features_*          ~2.38 M each    │
  │                          │       │  external_annotations_cache (live)   │
  │ Models loaded once at    │       │                                      │
  │ startup (xgb + ESM-2)    │       │  See docs/DATABASE.md for full schema│
  └──────────────────────────┘       └──────────────────────────────────────┘
                  │
                  ▼
  ┌──────────────────────────────────────────────────────────────────────────┐
  │ data/models/                                                             │
  │  binary_v5.json     — XGBoost binary classifier                          │
  │  ec_v5.json         — XGBoost EC multiclass (39 classes)                 │
  │  km_v5_final.json   — XGBoost Km regressor (log10 mM)                    │
  │  webapp/models/     — feature_names + ec_label_map mirrors               │
  └──────────────────────────────────────────────────────────────────────────┘
```

## Three-step ML cascade

CarboDB's classifier is **three sequential XGBoost models**, not one:

1. **Binary** — is this protein sequence a carboxylase at all?  
   Output: probability ∈ [0,1].
   AUROC: 0.9999 on held-out test.

2. **EC class** — *if* binary > threshold, which of 39 EC classes?  
   Output: softmax probabilities over 39 classes.
   Top-1 accuracy: 99.82%.

3. **Km regression** — *if* assigned to a CO₂-active EC class, predict
   the Michaelis constant for the natural substrate.  
   Output: log₁₀(Km in mM).  
   R² = 0.9253, Pearson r = 0.9628, RMSE = 0.4072 (log10 mM units).

This cascade matters because:
- Each model is trained on a feature stack that fits its task (binary
  takes everything; EC takes Pfam-heavy; Km takes ESM-2-heavy + dipeptides).
- A "this is a carboxylase but I don't know which one" answer is impossible —
  the EC step always picks one. We expose `ec_confidence` for honesty.
- Km is **predicted within an EC class context** — the regressor implicitly
  uses the EC assignment by being trained on labeled CO₂-active sequences.

## Feature stack — 1793 features

The full feature vector for one sequence:

| Group | Size | Source | Shape |
|---|---|---|---|
| Amino acid composition (AAC) | 20 | composition.py | per-letter freq |
| Dipeptide composition | 400 | composition.py | (20×20) flatten |
| Pseudo-AAC (PAAC) | varies | composition.py | hydrophobicity-weighted |
| Physicochemical | ~10 | composition.py | MW / pI / GRAVY / aromaticity / instability / log10 length |
| **Pfam HMM hits** | binary on ~17K | HMMER vs Pfam-A | accession + e-value + bitscore |
| **InterPro** | binary | InterProScan | PANTHER, Gene3D, TIGRFAM, SCOP, CDD, HAMAP, PROSITE |
| **ESM-2 embeddings** | 1280 | ESM-2 t33_650M | mean-pooled per-residue → 1280-vec |
| Expert motifs | ~7 | hand-coded regex | RuBisCO K-K, G-K; CA H-H, His cluster; PEPC R-R; Biotin M-K, A-M-K |

The Pfam e-value and bitscore are scalar features, not just presence
booleans — so the model can learn "how textbook is this enzyme" as a
continuous Km signal (this is what Panel A of the motif analysis matches).

## Request lifecycle — `/api/v1/predict` (live prediction)

```
POST /api/v1/predict        sequence={...}
  │
  │ webapp/app/routes/predict.py
  ▼
  predict_sequence(sequence, kingdom, mode)
  │
  │ webapp/app/pipeline/predict.py
  │   1. composition features    (~ms)
  │   2. HMMER scan vs Pfam-A    (~5–15 s for one sequence, single-threaded)
  │   3. InterProScan            (~30–60 s, calls java subprocess, often the bottleneck)
  │   4. ESM-2 forward pass      (CPU: ~20–40 s for ~500 aa, GPU: ~1 s)
  │   5. assemble 1793-vec
  │   6. xgb_binary.predict_proba
  │   7. if positive → xgb_ec.predict_proba
  │   8. if CO2-active EC → xgb_km.predict
  │   9. SHAP explanation (training-time precomputed where available)
  │   10. BLAST against per-EC db with experimental Km
  │
  ▼ JSON response with the full ResultDetail-shaped payload
```

**Total: ~100 seconds** for a 470-aa RuBisCO on wbbi206 (CPU). InterProScan
dominates. ESM-2 is the second cost. This is why /batch exists for >1
sequence and why Stage 2 of the deployment plan is critical.

## Request lifecycle — `/api/v1/db/seq/{uid}` (precomputed lookup)

For sequences already in the database, the prediction is precomputed and
stored in `predictions`. No model runs. Response shape matches /predict
(same JSON keys, same ResultDetail.vue rendering), but lacks SHAP and live
BLAST. Total: **~0.5 seconds**, mostly network + JSON serialization.

This is the "Database → click row → Details" path. It's the fast path. If
the user wants SHAP explanations they have to click "Re-analyze" which
re-runs through /predict.

## Request lifecycle — `/api/v1/external/{uid}` (annotation)

Pure proxy + cache. Hits UniProt REST + AlphaFold metadata API, merges into
one JSON, caches in `external_annotations_cache` for 30 days, returns
combined payload. PDB structure is a separate proxy endpoint
`/external/{uid}/structure` that fetches AlphaFold PDB and streams it.

Why proxy on the server side rather than from browser: avoids CORS issues,
hides UniProt query rate from each user (we have one shared rate-limit
budget), and the SQLite cache means repeated views of the same protein
don't hit upstream at all.

Versions tried for AlphaFold: `[6, 5, 4]` in order. AlphaFold bumped from
v4 to v6 between training and now; old hardcoded URLs would 404. The proxy
walks down the list until it gets a 200.

## Frontend architecture

Vue 3 with Composition API (`<script setup>`). Vite for dev/build. Pinia for
stores (kept light — only `searches.js` and `validation.js` exist).
Single-page-app with Vue Router.

```
src/
├── main.js                 entry; mounts App.vue with router + pinia
├── App.vue                 layout shell — header + <RouterView>
├── router/
│   └── index.js            5 routes — /, /analysis, /database, /about, /searches
├── views/
│   ├── HomeView.vue        landing
│   ├── AnalysisView.vue    /predict + /batch UI
│   ├── DatabaseView.vue    /browse UI with filters and paginated table
│   ├── AllSearchesView.vue history of past Analysis runs (Pinia)
│   ├── HmmerSearchView.vue legacy, mostly unused
│   └── AboutView.vue       static about
├── components/
│   ├── ResultDetail.vue        the giant shared "details panel"
│   ├── ResultDetailItem.vue    sub-block within ResultDetail
│   ├── ResultList.vue          batch results table
│   ├── ResultListItem.vue      single batch row
│   ├── ExtendedDetails.vue     UniProt + AlphaFold annotation (May 6)
│   ├── StructureViewer.vue     NGL Viewer wrapper with 4 color schemes
│   ├── SearchMenu.vue          input + mode dropdown
│   ├── FileUpload.vue          FASTA upload
│   ├── CommonButton.vue        styled button
│   ├── CommonModal.vue         modal shell
│   ├── Navigation.vue          top tab bar
│   └── TheHeader.vue           page-level header
└── stores/
    ├── searches.js         Pinia store — keeps Analysis history in localStorage
    └── validation.js       FASTA validation helpers
```

`ResultDetail.vue` is the single most important frontend component. It's
~750 lines, used by both Analysis and Database, accepts a `result` prop
that always conforms to one ResultDetail-shape. The Database mode adds
`mode === 'db_lookup'` and the Analysis mode adds SHAP + BLAST sections.
ExtendedDetails is wired in conditionally — currently only renders for
`db_lookup` mode (this is intentional and a good place to extend, see
ROADMAP.md).

## Key design decisions

- **One SQLite file, not Postgres.** 50 GB of mostly read-only data.
  Write-once, read-many. SQLite is fine if you set proper PRAGMAs
  (`cache_size = -524288`, `mmap_size = 30000000000`, `temp_store = MEMORY`).
  Postgres adds operational complexity we don't need.

- **Models loaded ONCE at startup.** ModelStore is a class-level singleton
  in `webapp/app/startup.py`. XGBoost JSON loads in <1 s, ESM-2 takes ~10 s.
  Reloading per-request would multiply our 100 s by 10× — keep the lifespan
  hook honest.

- **Synchronous /predict, not async/queued.** Easier to reason about, fine
  for low traffic. When we go multi-user public, this becomes a job queue
  (Redis + RQ or Celery + Redis). See DEPLOYMENT.md Stage 2.

- **/db/seq is the fast path; /predict is the slow path.** The Database tab
  uses /db/seq exclusively, which is why it's instant. The Analysis tab
  always goes through /predict because the user is asking for a new analysis,
  often with SHAP. Don't conflate the two paths.

- **External annotations are a separate concern.** UniProt+AlphaFold is
  fetched only when the user explicitly clicks "Show extended annotation".
  Lazy load → cheap response → no upstream rate-limit pressure for users
  who never click it.

- **No ORM.** Plain `sqlite3` with `row_factory = sqlite3.Row`. The schema
  is small enough to fit in two screens; an ORM would be more weight than
  benefit and would obscure the SQLite-specific PRAGMAs that matter.

## What's not in the architecture

- **No live retraining.** Models are static JSON files. Retraining is a
  manual `scripts/08_train_models.py` run. Future work: scheduled retrain
  pipeline (see ROADMAP).
- **No GPU in the running webapp.** ESM-2 runs on CPU. GPU dispatch is
  possible (env var `ESM2_DEVICE=cuda` if `wbbi203` had a free slot) but
  the model loads on CPU by default to avoid colocation issues.
- **No authentication.** Public endpoint. Anyone can submit any sequence.
  Rate-limited per-IP only on the external-annotation endpoint
  (60 req/min sliding window). /predict has no rate limit yet.
- **No structured logging.** Plain `logging.basicConfig` with stdout
  formatting. Fine for dev. Production should switch to JSON logs.
