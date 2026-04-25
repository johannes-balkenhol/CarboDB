# CarboDB Webapp — Remaining Work Plan

**Author:** session of April 20, 2026
**State at plan time:** Backend `ef625f5`, frontend `a00901e`. Analysis details panel works end-to-end with SHAP/motifs/phys. Four tasks remain.

---

## Current state snapshot (what the diagnostics revealed)

**Nearest-neighbor SQL** (`webapp/app/routes/predict.py`): currently does `ORDER BY RANDOM() LIMIT 8` over `reviewed=1` sequences with a predicted Km. Does not rank by sequence similarity, does not use experimental Km, does not distinguish BRENDA vs SwissProt tiers. Returns random reviewed entries.

**Batch endpoint** (`webapp/app/routes/batch.py`): **exists and is well-designed.** Takes `multipart/form-data` file upload, creates a job directory, runs annotation in a BackgroundTask, writes `results.tsv`. Has `/jobs/{job_id}` status polling and TSV download. But: the frontend AnalysisView.vue's `predictBatch()` sends JSON `{fasta: ...}` to `/batch`, which doesn't match the multipart/file contract. **So batch is broken at the frontend integration layer, not the backend.**

**DB data for nearest-neighbor** (`4.1.1.39` / RuBisCO as sample):
- BRENDA source: **153,936 sequences, 916 with experimental Km**, reviewed=0
- TrEMBL source: 497 sequences, 0 with experimental Km, reviewed=0
- **No SwissProt-tagged rows for RuBisCO in the `source` column.** The `reviewed` flag is 0 everywhere. This matters for the tier logic.

**DB indices**: good for EC-based lookups (`idx_seq_ec`, `idx_seq_label_ec`), good for km filtering (`idx_km_tier`, `idx_km_val`). Missing: no index on `km_best_mM` directly — may need one for the "has experimental Km" filter. Quick to add if query is slow.

**DatabaseView.vue**: 578 lines, already has substantial structure — stats header, search filters including query/ec_class, calls `browseDatabase()` and `getDatabaseStats()`. Does NOT call either on mount — only on user action. The masterplan's "show default results on load" is a ~5-line fix.

**Browse.js / Batch.js**: Browse.js is complete (`browseDatabase`, `getDatabaseStats`). Batch.js has three functions (`submitBatchJob` using FormData, `getJobStatus`). AnalysisView's inline batch code bypasses Batch.js entirely — another broken integration.

---

## Task 1 — Nearest-neighbor with curation tiers

### Scope
Replace the `ORDER BY RANDOM()` logic in `get_similar_from_db()` with tier-ranked nearest-neighbor. Return **two** neighbors when appropriate:
- **Primary** — best experimental-Km match (gold standard)
- **Secondary** — closer sequence match with predicted Km only, if it's significantly closer than the primary

Each returned neighbor carries a visible `tier` label.

### In scope
- New SQL query in `webapp/app/routes/predict.py` that joins `sequences` + `km_evidence` + `predictions` + `ec_evidence`
- Tier derivation: `brenda_experimental` > `swissprot_reviewed` > `trembl_reviewed` > `carbodb_predicted`
- Similarity ranking: use BLAST identity from `features_blast` table (column name TBD — needs one more diagnostic) if available; else fall back to length similarity as a cheap proxy
- Frontend: update `ResultDetail.vue` to show the tier badge and render experimental vs predicted Km distinctly
- Add the "Similar sequences" section population (currently shows "-" because `top_similar` is empty)

### Out of scope
- Full all-vs-all BLAST computation (too slow per-request — use pre-computed `features_blast`)
- New clustering / MMseqs work

### Files touched
- Backend: `webapp/app/routes/predict.py` (replace `get_similar_from_db`)
- Frontend: `frontend/src/components/ResultDetail.vue` (tier badge, two-neighbor layout)

### Input/output contracts

**New response shape for `top_similar`:**
```json
[
  {
    "rank": 1,
    "uniprot_id": "P04718",
    "organism": "Nicotiana tabacum",
    "ec_number": "4.1.1.39",
    "km_experimental_uM": 10.5,
    "km_experimental_mM": 0.0105,
    "km_predicted_uM": 9.8,
    "sequence_identity_pct": 94.2,
    "tier": "brenda_experimental",
    "tier_label": "BRENDA — experimental Km",
    "evidence_tier_numeric": 1
  },
  {
    "rank": 2,
    "uniprot_id": "...",
    "tier": "carbodb_predicted",
    "tier_label": "CarboDB — predicted Km",
    "sequence_identity_pct": 98.1,
    "km_predicted_uM": 10.1,
    "km_experimental_uM": null
  }
]
```

### Acceptance test
```bash
curl -s -X POST http://localhost:8090/api/v1/predict \
  -H 'Content-Type: application/json' \
  -d '{"sequence":"MSPQTETKAG...","mode":"fast","kingdom":"plant"}' \
  | python3 -c "
import json, sys
d = json.load(sys.stdin)
ts = d['top_similar']
assert len(ts) >= 1, 'no neighbors returned'
t0 = ts[0]
assert t0['tier'] == 'brenda_experimental', f'primary should be BRENDA, got {t0[\"tier\"]}'
assert t0['km_experimental_uM'] is not None, 'primary should have experimental Km'
print('PASS:', t0['uniprot_id'], t0['tier'], t0['km_experimental_uM'], 'µM')
"
```

### Commit boundaries
- Commit 1 (backend): `feat(api): nearest-neighbor with curation tiers` — predict.py patch only
- Commit 2 (frontend): `feat(ResultDetail): show curation tier on similar sequences` — component update

### Risks
1. `features_blast` table might not have entries for all sequences. If a user's query sequence is novel, there's no pre-computed similarity to DB. Fallback: length-similarity ranking or "no similar sequence found" tier.
2. The `source` column only has `brenda` / `trembl` values, not `swissprot`. The 3-tier naming in the plan may need revision — could be 2-tier: `brenda_experimental` / `carbodb_predicted`. Needs one more diagnostic: `SELECT DISTINCT source, reviewed FROM sequences`.
3. `km_evidence.evidence_tier` numeric values unknown. Assume 1=experimental (BRENDA), 2=curated-inferred, 3=predicted, but verify.

### Time estimate
2 hours backend + 1 hour frontend

---

## Task 2 — Usability patch (SHAP → my-sequence matching + plain-language summary)

### Scope
Make the SHAP feature importance panel interpretable to non-experts. Three additions to `ResultDetail.vue`:

1. **"Your seq" column** in SHAP feature rows showing whether user's sequence matches each top feature
2. **Plain-language summary block** at top of Feature Importance section explaining the prediction in 2-4 bullets
3. **Tooltips on feature names** via `title` attribute with feature family descriptions

### In scope
- `featureMatch(feat, result)` helper that interprets SHAP feature names (`pfam_PF02788`, `dp_KK`, `esm2_1083`, `aac_Y`, `motif_4111`) against the user's `pfam_hits` + `features_computed` and returns `{status, value, label}`
- Rendering: a compact `[✓ 3.3e-43]` / `[✗ —]` / `[n/a]` badge between feature name and bar
- Plain-language summary generator: computed property that picks top 3 SHAP features with clear interpretations, produces bullets with Pfam names mapped to human-readable labels
- Pfam-name dictionary for the top ~30 recurring accessions across carboxylase EC classes

### Out of scope
- Per-ESM-2-dimension interpretations (no meaningful content available)
- Interactive SHAP explorer / waterfall charts (future work)
- Training a natural-language model to generate summaries (rule-based is sufficient)

### Files touched
- Frontend only: `frontend/src/components/ResultDetail.vue` (one file, ~200 lines added)

### Input/output contracts
Purely a rendering change — no API changes. Component computes everything from existing `result.shap`, `result.features_computed`, `result.pfam_hits`.

### Acceptance test
Visual regression on the RuBisCO case:
- Summary block at top of Feature Importance section says "RuBisCO large-subunit domain" and "RuBisCO N-terminal domain" in bullets
- "EC class" SHAP tab row 1 shows `Pfam PF02788 [✓ 3.3e-43]` — green checkmark, actual e-value
- "EC class" SHAP tab row 3 shows `ESM-2 dim 1083 [n/a]` — gray, no interpretability at per-dim level
- "Km" tab row 7 shows `dipeptide NC [✓ 0.018]` — checkmark + frequency
- Hover over feature name → tooltip appears

### Commit boundaries
One commit: `feat(ResultDetail): interpret SHAP features against user sequence`

### Risks
1. Pfam-name dictionary needs curated entries. I have the top ~15 from earlier training-doc work (PF00016 RuBisCO_large, PF02788 RuBisCO_large_N, PF00194 CA_superfamily, PF00289 biotin_carboxylase_N, etc.). The rest fall through to "Pfam PF#####" as-is.
2. Summary generator needs to handle edge cases: user sequence hits ZERO of the top SHAP features. In that case summary should say "unusual prediction" and advise checking confidence.
3. Token cost: this is ~200 lines of frontend, easy to handle in one commit with the base64 approach we've been using.

### Time estimate
1-2 hours

---

## Task 3 — Batch search (fix broken integration, not build new thing)

### Scope
The backend `/batch` endpoint is complete. The frontend is broken: AnalysisView's `predictBatch()` sends the wrong format to the wrong contract. Fix the frontend to use the actual endpoint properly.

### In scope
- Rewrite `AnalysisView.vue` `predictBatch()` to:
  - Use the existing `submitBatchJob()` from `frontend/utils/commands/Batch.js` (FormData, multipart)
  - Accept FASTA text from the textarea: package it as a Blob with filename, not raw JSON
  - Show a "Queued (job {job_id})" state after submission
  - Poll `/jobs/{job_id}` every ~3 seconds for progress
  - When `status == "completed"`, fetch `/jobs/{job_id}/results.tsv`, parse rows, populate `batchResults`
- Progress indicator: show `progress_pct`, estimated time remaining
- Cancel button: stub (just stops polling on the frontend — full cancel is out of scope)
- Handle `status == "failed"` — display the error from the job JSON

### Out of scope
- Backend changes (it's fine)
- Email notification system (already stubbed in backend, not wired)
- Persistent job history in UI (see job IDs if user keeps them)

### Files touched
- Frontend: `frontend/src/views/AnalysisView.vue` (major rewrite of the batch path, ~80 lines changed)

### Input/output contracts
Use the documented backend contract from `batch.py`:
- `POST /batch` with `multipart/form-data`: `file`, `mode`, `kingdom`, optional `email` → `{job_id, status: "queued", n_sequences, estimated_minutes}`
- `GET /jobs/{job_id}` → `{status, processed, n_sequences, progress_pct, ...}`
- `GET /jobs/{job_id}/results.tsv` → TSV file

### Acceptance test
1. User pastes 5 FASTA sequences into the batch textarea, clicks Analyze
2. UI shows "Queued (job abc123) · 5 sequences · ~15s"
3. Progress bar animates 0% → 100%
4. Results table populates with 5 rows, each clickable → Details panel
5. Download TSV button downloads `carbodb_batch_abc123.tsv`

### Commit boundaries
One commit: `feat(batch): wire AnalysisView to real /batch+jobs endpoints`

### Risks
1. `run_batch_job` background task uses the same `annotate.py` subprocess we fixed earlier. Depends on uvicorn being launched from the project root (not webapp/). Same fix — already applied.
2. TSV output format from `batch.py` — we don't know what columns it emits. Need one diagnostic: run a small batch, inspect the TSV. If it's just prediction columns (no features_computed or SHAP), details panel for batch results will be thinner than for single-sequence — acceptable, since batch is about scale.
3. Polling overhead: if user leaves the page open, polling continues. Add cleanup on component unmount.

### Time estimate
2-3 hours (depends on what the TSV format is)

---

## Task 4 — Database browsing (improve existing view)

### Scope
`DatabaseView.vue` already exists with search filters and stats header. Two gaps from the masterplan:

1. "Show default results on load" — 5-line fix
2. "Find nearest neighbor in database" — interpretation question (see below)

### Ambiguity to resolve first
You said "find nearest neighbor in database" — I can read this two ways:

**A.** Add a "paste a sequence, find its closest DB matches" feature on the Database page itself (mini standalone tool).
**B.** Just make the existing browse-by-filter work correctly and add a "find similar" button next to each row that takes you to the Analysis page with that sequence pre-loaded.

**Recommended:** B. It reuses the existing nearest-neighbor logic from Task 1, avoids feature creep on the Database page, and nudges users toward the richer Analysis view. A is a larger ticket and probably doesn't add enough over just running Analysis on a pasted sequence.

### In scope (assuming B)
- `onMounted` in DatabaseView: call `browseDatabase({limit: 50})` to show default rows
- Add "Analyze similar" button on each row: navigates to `/analysis` with the sequence pre-filled (via Vue Router query param or Pinia store)
- Empty-state message: "No matches — try clearing filters"
- Loading skeleton while fetching

### Out of scope
- Server-side pagination (existing `limit`/`offset` in `/browse` is enough for v1)
- Column sorting (nice-to-have, skip for now)
- Export filtered results (future task)

### Files touched
- Frontend only: `frontend/src/views/DatabaseView.vue`

### Acceptance test
1. Navigate to `/database` → table shows 50 rows immediately (no user action needed)
2. Filter by EC class = 4.1.1.39 → table updates to only RuBisCO entries
3. Click "Analyze similar" on a row → Analysis page loads with that sequence in the textarea

### Commit boundaries
One commit: `feat(database): default results on mount + analyze-similar navigation`

### Risks
1. Route navigation with sequence payload: if the sequence is long (>2000 aa), URL query param is too big. Use Pinia store or router `state`. Already pattern-matched in `stores/searches.js` (existing file).
2. `browseDatabase({limit: 50})` — verify `/browse` accepts `limit` param (based on Browse.js it does).

### Time estimate
1-2 hours

---

## Recommended execution order

**Reasoning:** maximize user-visible value per unit of work, minimize task dependencies.

1. **Task 3 (batch) — FIRST.** Reason: batch is completely broken right now (sends wrong format, hits wrong contract). Highest functional gap. Enables you to demo "analyze 20 sequences at once" to reviewers.
2. **Task 1 (nearest-neighbor tiers) — SECOND.** Reason: single highest-impact feature for scientific credibility. Answers "is this prediction trustworthy?" in a way the UI currently can't. Needs one more DB diagnostic before execution.
3. **Task 2 (SHAP usability) — THIRD.** Reason: polish on what's already working. Non-blocking. Can land after publication-deadline pressure eases.
4. **Task 4 (database view) — LAST.** Reason: existing DatabaseView is partially usable. The "default results on load" piece is trivial. The "analyze similar" depends on Task 1 being done.

Dependencies:
- Task 4 "Analyze similar" button depends on Task 1 completing for end-to-end value.
- Tasks 2 and 3 are independent of everything else.

---

## Session-hygiene recommendation

This conversation has run long and paste-corruption bugs have cost us multiple rounds. For each task above, I recommend **starting a fresh Claude conversation** with:

1. A short paste of the current git state (`git log --oneline -5` from both repos)
2. This plan file as context
3. The specific task number to execute

Each task's spec above is self-contained enough to kick off cold. The backend + frontend patches we've made today (commits `ef625f5`, `c0376c4`, `a00901e`) are already in git so a new session just needs to see the HEAD state.

---

## Pre-task-1 diagnostic

Before Task 1 can execute, run this and include in the opening message:

```bash
cd ~/Projects_shared/CarboDB_v3

# Which tiers/sources actually exist in the DB
sqlite3 data/primary/carbodb.sqlite "
SELECT source, reviewed, COUNT(*) FROM sequences GROUP BY source, reviewed;"

# What evidence_tier values exist in km_evidence
sqlite3 data/primary/carbodb.sqlite "
SELECT evidence_tier, source, COUNT(*) FROM km_evidence GROUP BY evidence_tier, source;"

# features_blast schema
sqlite3 data/primary/carbodb.sqlite ".schema features_blast"

# Does carbodb store sequence identity per pair, or just nearest-neighbor?
sqlite3 data/primary/carbodb.sqlite "SELECT * FROM features_blast LIMIT 1;"
```

These three outputs answer the open questions in Task 1's risk section.


---

# Future work (noted April 20, 2026 — out of scope for current sessions)

## Task 5 — Structure prediction / AlphaFold integration for input sequences

### Vision (from user)
Give users a 3D structure view of their submitted sequence — not just the nearest neighbors. Three cases to handle:

1. **Known UniProt ID** (rare — user pasted a reference sequence): fetch from EBI AlphaFold API at `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}`, returns PDB URL. Render with NGL Viewer or 3Dmol.js embedded in ResultDetail.vue.

2. **Novel sequence matching a known UniProt 100%**: detect via BLAST (Task 1), show structure of the matched UniProt entry with a note "structure of closest match".

3. **Truly novel sequence**: run ESMFold or ColabFold to predict structure on demand. ESMFold is simpler (single-sequence, fast, ~30s-5min CPU), ColabFold is more accurate but needs MSA + template search (minutes to hours). Most reasonable for a webapp is ESMFold for inputs <400 aa; cache predictions in a `structures/` directory keyed by sequence hash.

### Scope for a future session
- Backend: new endpoint `/api/v1/structure` that takes sequence, returns PDB string
  - If sequence is in DB and matches a known UniProt → fetch AlphaFold
  - Else → run ESMFold subprocess (cached)
- Frontend: expandable "3D Structure" section in ResultDetail.vue and per-neighbor
  - 3Dmol.js viewer (more accessible than NGL for bio users)
  - Color by confidence (pLDDT)
- Batch implications: on-demand structure prediction is too slow for 500-sequence batches. Option: for batch, only show AlphaFold for top hits (neighbors with known UniProt IDs).

### Risks / dependencies
- ESMFold requires ~5 GB model download and PyTorch with CUDA for reasonable speed
- Disk usage for structure cache grows quickly — need a cache eviction policy
- Not viable as synchronous API call for long sequences; may need job queue similar to batch

### Time estimate
1-2 full days, needs GPU or patient users

## Task 6 — General protein features page

### Vision (from user)
An additional page/inlet showing richer protein features, expanding beyond what ResultDetail already has. Not yet specified concretely — potential additions:
- Secondary structure prediction (DSSP if structure available, or JPred / PSIPRED)
- Transmembrane topology (Phobius, DeepTMHMM, TMbed)
- Signal peptide prediction (SignalP)
- Subcellular localization (DeepLoc, TargetP)
- Predicted functional sites / active-site residues
- InterPro domains beyond Pfam
- Conservation analysis (if MSA available)

### Questions to resolve before scoping
- Which of these features are actually priorities for your users? (Define 2-3, not all 7.)
- Do any require external API credentials / rate limits?
- Integrate inline in ResultDetail.vue (expandable sections) OR as a separate "Protein Features" page tab?

### Time estimate
Unknown until features are selected — anywhere from 1 day (single added feature from open-source CLI) to 2 weeks (multiple features with different tooling).
