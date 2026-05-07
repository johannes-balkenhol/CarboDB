# CarboDB — Co-Developer Handoff

**Date:** 2026-05-07  
**Audience:** colleague joining as shared owner of the codebase  
**Read this first.** Everything else in `docs/` is a deeper dive into one slice.

---

## What CarboDB is, in one paragraph

A protein database + ML pipeline + web app for CO₂-fixing carboxylases.
We ingest BRENDA + UniProt, annotate features (composition, Pfam, InterPro,
ESM-2 embeddings, CD-HIT clusters, BLAST), train three XGBoost models
(binary carboxylase / EC class / Km regression), apply them to all 2.38M
sequences, and serve everything via a FastAPI + Vue 3 web app at
http://132.187.22.206:5173. Public-facing, predictions live, model
performance AUROC 0.9999 / EC top-1 99.82% / Km R² 0.9253.

## Two repositories, one running system

| Repo | What | Where |
|---|---|---|
| **CarboDB** | backend, ML, data pipeline, SQLite | `~/Projects_shared/CarboDB_v3` |
| **CarboDB-App** | Vue 3 frontend | `~/Projects_shared/CarboDB-App-v2` |

The frontend repo's path has a `-v2` suffix; the directory inside it is just
`frontend/`. There's no `-App-v3` planned — we kept v2 and rewrote internally.

Server: `wbbi206`, user `job37yv`, conda env `carboxylase` (Python 3.11).
Total project footprint: 285 GB (data 284 GB, webapp 28 MB, frontend 191 MB).

## How to start the app

```bash
~/Projects_shared/CarboDB-App-v2/start_app.sh restart
```

This kills any running uvicorn + vite, restarts both, and prints the URLs.
- Backend (FastAPI): http://132.187.22.206:8090
- Frontend (Vue/Vite): http://132.187.22.206:5173

The script is short — read it (`cat start_app.sh`) to understand what it
does. It sets `DB_PATH`, `PFAM_HMM`, `MODELS_DIR`, `JOBS_DIR`, `ESM2_DEVICE`
as env vars, then forks `uvicorn` + `npm run dev` into the background.

Logs:
- Backend: `~/Projects_shared/CarboDB_v3/webapp/logs/webapp.log`
- Frontend: `/tmp/vite.log`

## Where to look first

In rough order of "useful for new contributor":

| Question | File |
|---|---|
| Big-picture architecture, data flow, model layout | `docs/ARCHITECTURE.md` |
| How the data was ingested, scripts 01–23 | `docs/DATA_INGESTION.md` |
| Backend route catalog + frontend route catalog | `docs/WEBAPP.md` |
| Install on a fresh server | `docs/DEPLOYMENT.md` |
| What's done, in progress, planned | `docs/ROADMAP.md` |
| Motif analysis status + Km/AA evaluation strategy | `docs/MOTIF_ANALYSIS_v2.md` |
| Current SQLite schema | `docs/DATABASE.md` (existing, still accurate) |
| API endpoints (older, partly stale) | `docs/API.md` |
| Frontend spec (older, says Flask — actually FastAPI) | `docs/FRONTEND_SPEC.md` |
| External annotation integration (added May 6) | `docs/EXTERNAL_INTEGRATION.md` |

The older docs (API.md / FRONTEND_SPEC.md / PIPELINE_PLAN.md) are correct
on big-picture concepts but stale on a few details (FRONTEND_SPEC says Flask;
DATABASE.md says ~30 GB; PIPELINE_PLAN.md was last touched Mar 30). The new
HANDOFF/ARCHITECTURE/WEBAPP/DEPLOYMENT/ROADMAP docs are authoritative as of
2026-05-07.

## What was just shipped this week (May 5–7)

1. **Motif bundles v2 cleaned** — removed Helianthus 18 mM HCO₃⁻ contamination
   and 19 stuck XGBoost values from the high-Km set. New tarball at
   `data/motifs_v2_clean.tar.gz`.
2. **Database browser rewrite** — search/filter/paginate/sort + `/api/v1/stats`
   banner aggregates. `webapp/app/routes/browse.py` was fully rewritten.
3. **Performance fixes** — SQLite cache PRAGMAs (512 MB) + `/stats` cache
   (10-min TTL). `/stats` warm path is now 15 ms (was 60 s).
4. **External annotation integration** — `/api/v1/external/{uniprot_id}` proxy
   pulls UniProt + AlphaFold, lazy-loaded behind a button on the Details panel.
   Renders function text, GO terms, taxonomy, cross-refs, **3D structure
   via NGL Viewer**, color modes (pLDDT / Pfam / motifs / rainbow).
5. **Bug fixes** — uniprot_id field naming in `/db/seq` response;
   AlphaFold v6 PDB URL handling; FastAPI HEAD method support for the
   structure proxy.

Full git history with descriptions is in `docs/ROADMAP.md` (Done section).

## What needs the most attention right now

In rough priority order:

1. **The HCO₃⁻ contamination story is bigger than v2 cleanup.** ~20–30 of the
   600 BRENDA Km entries for RuBisCO are problematic (HCO₃⁻ measurements
   relabeled as CO₂, mutant entries broadcast at organism level, 26 distinct
   UniProt IDs all sharing exactly 18.0 mM). v2 cleaned this for the motif
   bundles but the **trained model still sees these in km_evidence**, and a
   v6 retrain is needed to fix it properly. See `docs/ROADMAP.md` →
   "BRENDA ingestion redesign + v6 retrain" — this is multi-day work
   that should be scoped together first.

2. **19 stuck XGBoost predicted-Km values** are visible in the database.
   E.g., 15.5278882980347 mM appears identically across 400 sequences;
   2.68330836296082 across 139 sequences. These are tree-ensemble leaf-grid
   degenerate outputs at sparse training tail, not contamination. Webapp
   should show a warning indicator on Detail panel for these (not done yet).

3. **Form I vs Form II/III phylogeny confounds the motif analysis.** Current
   low-Km/high-Km comparison mixes Form I plant RuBisCO with Form II/III
   archaeal/bacterial. Within-Form-I split is the next analysis to run; bundle
   ready in `data/motifs_v3_form_split/`.

4. **Live prediction is a 100 s blocking call.** /predict for one sequence
   runs HMMER + InterProScan + ESM-2 sequentially in-process. There's no job
   queue, no batching, no GPU dispatch. Fine on wbbi206 for one user at a
   time; will not survive multi-user public access. See deployment plan in
   `docs/DEPLOYMENT.md` § Stage 2.

5. **Browse query is still 12 s.** SQLite PRAGMAs got it down from 17 s. The
   real fix is a query rewrite — full join scan over predictions is the cost.

## What to do *first* on day 1

1. Read this file end-to-end (you're doing it).
2. Read `docs/ARCHITECTURE.md` for the system shape.
3. SSH into wbbi206 → `~/Projects_shared/CarboDB-App-v2/start_app.sh restart`
   → open http://132.187.22.206:5173 → click around the three tabs (Home,
   Analysis, Database) → poke the example sequences.
4. Open `docs/DATA_INGESTION.md` and `scripts/`. Read scripts 01 → 11 in
   order to understand how the database was built.
5. Pick **one item from `docs/ROADMAP.md` → In-Progress / Planned** and
   discuss with Johannes before starting.

## Communication

- Both repos on GitHub: `johannes-balkenhol/CarboDB` and
  `johannes-balkenhol/CarboDB-App`
- Commit message convention: `<type>(<area>): <imperative summary>` —
  examples in `git log --oneline -20`. Types in use: `feat`, `fix`, `perf`,
  `docs`, `chore`, `polish`, `refactor`. Areas: `external`, `browse`,
  `database`, `analysis`, `batch`, `motifs`, etc. Bodies are wrapped at 72.
- Backup files convention: `*.before_*` (gitignored) — a quick restore-able
  snapshot before risky edits. Delete after the change verifies green.

## Two warnings about the host

- `wbbi206` is described in the MOTD as "interactive jobs only, do not occupy
  >50% CPU/RAM". Long-running predictions for every-arriving-user *will*
  violate this. We're currently below the threshold but a multi-user public
  deployment must move computation to `julia2` or a job queue.
- `/tmp` is volatile. Batch-prediction artifacts and uploaded FASTAs land
  there briefly; nothing important should rely on `/tmp` persistence.

## When in doubt

Check `docs/ROADMAP.md` for the live to-do list. Anything that looks weird
in the running system probably has a known cause documented there. If it
isn't, write it down and add it to the Roadmap — that's how we keep the
project legible to future-us.
