# Roadmap

What's done, what's in progress, what's planned. Updated 2026-05-07.

This is a living list. When you finish or start something, update this
file in the same commit.

---

## Done (recent — last 2 weeks)

### Backend

- ✅ **`/api/v1/external/{uniprot_id}` proxy** — UniProt + AlphaFold
  metadata, 30-day SQLite cache, stale-cache fallback, sliding-window
  rate limit (60/min/IP). Validates UniProt ID against `sequences`
  table before any upstream call. *(commit b7d2612)*
- ✅ **`/api/v1/external/{uniprot_id}/structure` PDB proxy** — fetches
  AlphaFold PDB, walks model versions `[6, 5, 4]` for compatibility,
  serves with `chemical/x-pdb` content-type. Both GET and HEAD. *(b7d2612)*
- ✅ **`/api/v1/browse` rewrite** — search/filter/paginate/sort.
  *(536114d)*
- ✅ **`/api/v1/stats` banner aggregates with EC distribution** *(536114d)*
- ✅ **`/api/v1/db/seq/{uniprot_id}` instant ResultDetail-shaped JSON
  for stored sequences** *(536114d)*
- ✅ **SQLite cache PRAGMAs + `/stats` 10-min process cache** —
  `/stats` warm 15 ms (was 60 s); `/browse` 12 s (was 17 s). *(bc38428)*
- ✅ **`uniprot_id` field added to `/db/seq` response** — fixes the
  cdb_query_id confusion. *(2026-05-07)*

### Frontend

- ✅ **Database browser page rewrite** — 5-card stats banner,
  6 quick-pick examples, filters, paginated sortable table. *(0602ba1)*
- ✅ **ExtendedDetails component** — lazy-loaded UniProt + AlphaFold
  panel with function text, GO terms, taxonomy lineage, cross-refs.
  *(c737519)*
- ✅ **StructureViewer component with NGL Viewer** — 3D model with four
  color schemes (pLDDT / Pfam / motifs / rainbow), reset view, PDB
  download. *(c737519)*
- ✅ **NGL Viewer Vue 3 reactivity workaround** — using plain `let` not
  `ref()` for Stage/Component/Representation objects. *(c737519)*

### Data / pipeline

- ✅ **Motif bundles v2 cleaned** — Helianthus 18 mM HCO₃⁻ contamination
  removed (26 sequences); 19 stuck XGBoost predicted-Km values excluded
  from high-Km set. Tarball at `data/motifs_v2_clean.tar.gz` (19 MB).
- ✅ **Motif bundle generator script** — reproducible build with seed=42.
  `scripts/build_motif_bundles_v2.py`.
- ✅ **Per-EC BLAST databases for nearest-neighbor search** —
  `webapp/scripts/build_ec_blast_dbs.py`. *(453f51e, 882a549)*

---

## In progress (started, not finished)

### Documentation

- 🟡 **This batch of handoff docs** — HANDOFF / ARCHITECTURE / WEBAPP /
  DEPLOYMENT / DATA_INGESTION / MOTIF_ANALYSIS_v2 / ROADMAP. *(this PR)*
  Older docs (API.md, FRONTEND_SPEC.md, DATABASE.md, PIPELINE_PLAN.md)
  remain but should be cross-referenced from the new master docs.

### Motif analysis

- 🟡 **Within-Form-I motif bundle generated** — to test whether the
  Km-discriminating residues survive controlling for phylogeny. See
  `data/motifs_v3_form_split/` and `MOTIF_ANALYSIS_v2.md`. *(this PR)*
  **Awaiting:** colleague's analysis run on the new bundle.

### Webapp UX

- 🟡 **Detail panel "in-progress" UX** — during /batch processing,
  unfinished sequences currently show "0% / Not carboxylase / NaN" which
  looks like a confident wrong answer. Should show "Pending" or spinner.
  Fix is small (one v-if in ResultListItem.vue) but not yet applied.

---

## Planned (in priority order)

### High priority — core science

#### 1. BRENDA ingestion redesign + v6 retrain

The May 5 motif cleanup is post-hoc — it removes contamination from the
*motif bundles* but the *trained model* still saw them. Real fix requires:

- Redesign `scripts/01_brenda_download.py` to:
  - Separate WT vs mutant entries
  - Detect HCO₃⁻ measurements vs CO₂ measurements (substrate field
    contains "HCO3" or "bicarbonate" → flag)
  - Convert HCO₃⁻ → CO₂ at measurement pH using Henderson-Hasselbalch
    (pKa ≈ 6.35) when conditions allow
  - Exclude ambiguous entries
- Re-train models 8-stage (model v6)
- Re-predict all 2.38 M sequences (script 10)
- Validate against held-out experimental Km

Estimated: ~7 days of work, including the model training time.

**Action:** scope this on a call with Johannes before starting.

#### 2. Within-Form-I motif analysis

Once colleague reports back on the within-Form-I bundle, do:

- Map any Form-I-confirmed Km-discriminating residues onto 1RCX in PyMOL
- Cross-reference with SHAP feature importance from the v5 model
- Check whether dipeptides (YK, QP, TD, WT) align with residue
  changes within Form I (the report's Concerning Question 1)
- Co-occurrence analysis: are the changes independent or do they
  always travel together?

#### 3. Tara Oceans metagenome scan

Pending: link from Johannes for the gene catalog (likely OM-RGC.v2 or
the Tara Microbiome Atlas). Once we have it:

- Scan with binary classifier (carboxylase yes/no)
- For positives, classify EC
- For CO₂-active EC, predict Km
- Cross-reference with sample metadata (depth, temperature, latitude,
  oligotrophic vs eutrophic)
- Generate publication figure (extends `scripts/22_metagenome_scan.py`
  and `23_tara_figure.py`)

### High priority — webapp engineering

#### 4. Live-prediction backend pipeline

Convert /predict from synchronous to async with a job queue. See
`docs/DEPLOYMENT.md` § Stage 2. Required before public release.

Substeps:
- Redis + RQ setup
- Refactor `predict_sequence()` to be a worker task
- Add `/predict/{job_id}` polling endpoint
- Frontend: change AnalysisView to submit + poll, with cancellable UX
- ESM-2 on GPU as a separate worker pool

#### 5. Browse query rewrite

The 12 s /browse default time is dominated by a full-scan join over
`predictions`. Two avenues:
- Materialized view: pre-compute the join result into a flat table,
  refresh nightly (or on retrain).
- Index on the actual filter columns: many filter combos work, others
  don't, depending on which indexes get picked.

Aim: <2 s for typical filter+sort.

#### 6. AnalysisView ExtendedDetails

Currently ExtendedDetails only renders for `mode === 'db_lookup'`.
For Analysis page results with a strong BLAST hit (e.g. >95% identity),
it would be valuable to show ExtendedDetails for the **nearest hit**
not the user's sequence — that hit always has a UniProt ID and is
highly relevant context.

Implementation: two-line change in ResultDetail.vue's v-if, plus a prop
`:uniprot-id="result.top_similar?.[0]?.uniprot_id"` for live mode.

### Medium priority — webapp polish

#### 7. Stuck-value warning indicator

The 19 known degenerate XGBoost predicted-Km values should show a
warning badge on the Detail panel: "⚠ Predicted Km is in a known
degenerate region of the model (n=400 sequences share this value).
Treat as low-precision."

#### 8. Sort by experimental Km

`webapp/app/routes/browse.py` already supports `km_exp_asc` and
`km_exp_desc` sort modes. Frontend table header doesn't yet expose
the click. ~5-line change in DatabaseView.vue.

#### 9. Stats definitions doc

`docs/STATS_DEFINITIONS.md` — canonical SQL queries behind the 5 banner
numbers (Total, Predicted carboxylases, With experimental Km, EC classes,
SwissProt-curated). For when we get questions like "why does this number
disagree with X".

#### 10. AlphaFold-coupled per-residue SHAP coloring

Currently StructureViewer can color by Pfam domain spans or by motif
positions, but not by per-residue feature contribution. Per-residue SHAP
is computationally heavy (one ESM-2 occlusion pass per residue) and
requires precomputation. Useful but expensive — defer until #4 is done
and we have a queue to run it on.

#### 11. Tooltip system upgrade

Native `title=` tooltips are simple but limited (no rich content, no
positioning control). Replace with floating-vue or similar. Affects
tooltips throughout the app.

#### 12. Motif library expansion

Currently 7 hand-coded motifs (RuBisCO K-K, G-K; CA H-H, His cluster;
PEPC R-R; Biotin M-K, A-M-K). Plan: extend to all 39 EC classes via
either:
- Hand curation with a domain expert (slower, higher quality)
- MEME/STREME on the per-EC FASTAs from the motif bundles (faster, may
  need filtering)

### Lower priority — nice-to-haves

#### 13. Old-aesthetic restoration

Pre-rewrite UI had a different look that some users liked. Bring back
selected elements (e.g., gradient header, EC badges) without losing
the new features. Needs design refs.

#### 14. Field-naming cleanup

`cdb_query_id` is a misleading name in the /db/seq response — it holds
the CDB internal ID, not a UniProt ID. Renaming requires coordination
across both repos. Can be done as a single PR touching backend response
shaping + all frontend references.

#### 15. SHAP coverage for all 39 EC classes

`features_computed.shap` currently shows "No SHAP data available for
this EC class" for some EC classes. Either precompute SHAP for every
class, or hide the section when missing.

#### 16. Automated retraining pipeline

Currently retraining is a manual `scripts/08_train_models.py` invocation.
For v7 and beyond, automate:
- Trigger: new BRENDA release detected
- Steps 01–10 run as a Snakemake or Nextflow pipeline
- Artifacts get versioned
- Backed-up old SQLite kept for ~30 days
- Optional smoke-test stage before promoting to prod

#### 17. Structured logging

Replace `logging.basicConfig` with JSON-formatted logs to stdout, ready
for ingestion by a log shipper (Loki, ELK, whatever). Required for
production observability.

#### 18. CI / tests

There's a `__tests__/HelloWorld.spec.js` skeleton in the frontend but
no real test coverage. At minimum:
- Backend: smoke test for each /api/v1/* endpoint with known inputs
- Frontend: Vitest for the validation store and a few critical components
- GitHub Actions workflow that runs both on every push to `main`

---

## Cancelled / decided not to do

- **Postgres migration** — discussed, decided against. SQLite + PRAGMAs
  is enough for our read-heavy workload and avoids operational complexity.
- **Authentication** — public site, no auth needed. Rate-limiting is the
  right tool; per-IP sliding window already on /external.
- **Microservices split** — the prediction pipeline + DB + frontend would
  not benefit from being microservices at our scale. Keep it monolithic.

---

## How to add to this list

When you start something:
- Move it from Planned to In Progress
- Add your name and start date

When you finish:
- Move it to Done with the commit hash
- Note any caveats (e.g. "fixed for db_lookup, live mode still TBD")

When you decide not to do something:
- Move it to Cancelled with a one-line reason
