# Session Summary — 2026-05-07

This is the answer to your three asks. Each section is short; the
detail lives in the docs/ folder.

---

## Ask 1 — finish the app + handoff docs

### What got finished today

- **Field-naming bug fixed:** `/db/seq` now returns `uniprot_id` explicitly
  alongside `cdb_query_id`. Frontend uses `uniprot_id || result.id` for
  ExtendedDetails (no more 404 on CDB internal IDs).
- **NGL Viewer Vue 3 reactivity bug:** documented in WEBAPP.md so the
  colleague doesn't hit it again. Use plain `let stage = null`, never
  `ref()` for NGL/Three.js objects.
- **Extended annotation panel verified working** for db_lookup mode:
  function text, GO terms in 3 columns, taxonomy lineage, cross-refs,
  3D structure with pLDDT coloring (mean pLDDT 98.6).

### Handoff docs delivered (in `docs/`)

| File | Lines | Purpose |
|---|---|---|
| `HANDOFF.md` | 162 | First read — overview, where to look, what's pending |
| `ARCHITECTURE.md` | 239 | Data flow, 3-step ML cascade, request lifecycles |
| `DATA_INGESTION.md` | 332 | Pipeline scripts 01–23, what they do, how to re-run |
| `WEBAPP.md` | 343 | Backend + frontend route catalog, debugging hints |
| `DEPLOYMENT.md` | 372 | Stage 1 (current) + Stage 2 (production plan) |
| `ROADMAP.md` | 279 | 18 items: done / in-progress / planned, prioritized |
| `MOTIF_ANALYSIS_v2.md` | 309 | What v2 found, v3 plan, AA-shift evaluation strategy |

Total ~2,036 lines / ~85 KB. Sized for the colleague (shared owner)
audience — direct and technical, no marketing fluff. Older docs
(API.md, FRONTEND_SPEC.md, DATABASE.md, PIPELINE_PLAN.md) are kept and
referenced from these.

### What you asked about that ended up in ROADMAP.md

All your concerns made it into ROADMAP.md as planned items, in priority order:
- **#1 BRENDA ingestion redesign + v6 retrain** (HCO₃⁻ contamination + automated full pipeline for updates)
- **#2 Within-Form-I motif analysis** (using the v3 bundle below)
- **#3 Tara Oceans metagenome scan** (waiting for your link)
- **#4 Live-prediction backend pipeline** (job queue, Redis+RQ, async /predict)
- **#5 Browse query rewrite** (12s → <2s)
- **#10 AlphaFold-coupled per-residue SHAP coloring** (the "color the motifs" item — depends on #4 for compute)

I added a few items you didn't mention but should have:
- **#6 AnalysisView ExtendedDetails extension** (show external annotation
  for the nearest BLAST hit when user submits a novel sequence)
- **#7 Stuck-value warning indicator** on Detail panel for the 19 known
  degenerate XGBoost outputs
- **#9 STATS_DEFINITIONS.md** (canonical SQL behind banner numbers — for
  when someone asks why a number disagrees with their own count)
- **#14 cdb_query_id → uniprot_id rename** across both repos
- **#16 Automated retraining pipeline** (Snakemake or Nextflow for v7+)

You asked: "did I miss anything?" — three things that aren't in the
roadmap but should be on your radar:

1. **No authentication, no /predict rate limit.** Currently anyone on the
   internet can submit any sequence and pin a CPU for 100 s. Fine for
   internal demo, not fine when the link gets shared. /external has a
   60/min/IP limit; /predict has none. Adding a simple per-IP
   sliding-window limiter to /predict is ~30 lines.
2. **The 50 GB SQLite has no backup.** wbbi206 is interactive-use-only
   per the MOTD. A bad shutdown could corrupt it. Daily
   `sqlite3 carbodb.sqlite ".backup ..."` cron is ~5 minutes to set up.
3. **No `vite build` for production.** Currently the frontend serves via
   `npm run dev`, which is fine for a demo but slow to first-paint and
   not optimized. `vite build` produces a static `dist/` that nginx can
   serve in milliseconds. Step in DEPLOYMENT.md Stage 2.

---

## Ask 2 — new motif bundle from what we learned

### Delivered

`scripts/build_motif_bundles_v3_form_split.py` — generates three FASTAs:

- `form_I_low_km_pred.fasta` — Form I only, predicted Km < 0.01 mM
- `form_I_high_km_pred.fasta` — Form I only, predicted Km 0.1–5 mM,
  with the 19 stuck XGBoost values excluded
- `form_II_III_pooled.fasta` — Forms II + III together for contrast
- `form_classification_summary.tsv` — sanity check the assignments
- `README.md` — protocol the colleague should follow

Form classification uses InterPro PANTHER family `PTHR42704` (Form I)
plus `IPR017443` (Form II) + regex on `raw_ipr_json` for Form III/IV-RLP.
Where InterPro is silent, the sequence is labeled "Form ?" and excluded
from bundles (~5% of sequences). Reproducibility seed = 42, max 10
sequences per (form, genus) cell.

### Why this bundle, not "more sequences"

The colleague's v2 report itself flags the main concern (their
"Concerning Question 1"): the low-Km vs high-Km contrast might be
phylogenetic (Form I plant RuBisCO) rather than Km-specific. Sending
"more sequences" without addressing that confound would just produce
the same problem at higher resolution.

The right next experiment is **within-Form-I split**. If the v2
positions (251, 255, 258, 86, 391, 445, 461, 467) stay significant
within Form I, the signal is genuinely Km-related. If they vanish,
v2 was phylogenetic. Either outcome is publishable and clarifies the
biology.

### AA-shift evaluation strategy (in MOTIF_ANALYSIS_v2.md, summary)

For each surviving Form-I residue:

1. **Confirm Km association is Form-I-specific** — repeat Fisher's exact
   within Form I, Bonferroni over 470 positions
2. **Map onto 1RCX in PyMOL** — visualize the "second shell" around the
   catalytic core
3. **Categorize by mechanism** — charge change / hydrophobicity change /
   backbone flexibility / steric / conservative
4. **Check motif overlap** — compare against the 7 expert motifs, Pfam
   spans, FIMO hits, loop 6 (332–338)
5. **Co-occurrence pattern table** — independent (engineerable
   per-residue) vs coupled (likely phylogenetic at sub-Form level)
6. **Tie back to SHAP** — match dipeptides (YK, QP, TD, WT) to specific
   substitutions; check Pfam e-value as conservation-proxy; probe ESM-2
   dim 1083, 1059, 448 for residue correlations

### How to run

On the server:
```bash
cd ~/Projects_shared/CarboDB_v3
python scripts/build_motif_bundles_v3_form_split.py \
    --db   data/primary/carbodb.sqlite \
    --out  data/motifs_v3_form_split

tar czf data/motifs_v3_form_split.tar.gz -C data motifs_v3_form_split
```

Tarball goes to the colleague with the v3 README.

---

## Ask 3 — Tara Oceans metagenome scan

**Pending: the link.**

Once you paste it, the plan is:

1. **Identify the catalog format.** Most likely OM-RGC.v2 (~47 M genes,
   Salazar 2019) or the newer Tara Microbiome Atlas. Format dictates
   whether we batch through the existing webapp /batch endpoint or
   write a dedicated streaming scanner.
2. **Re-use `scripts/22_metagenome_scan.py`** as a starting point — it
   already handles the model cascade application; needs minor adaptation
   to whatever the new catalog's headers look like.
3. **For each gene:** binary classify → if positive, classify EC →
   if CO₂-active EC, predict Km.
4. **Cross-reference with sample metadata** — depth, temperature, latitude,
   oligotrophic vs eutrophic.
5. **Generate publication figure** — extends `scripts/23_tara_figure.py`.

Estimated runtime depends on catalog size. For ~47 M genes, full scan
on the wbbi206 CPU is ~3-5 days; on a GPU node ~12 hours.

---

## Files presented with this summary

```
docs/
├── HANDOFF.md
├── ARCHITECTURE.md
├── DATA_INGESTION.md
├── WEBAPP.md
├── DEPLOYMENT.md
├── ROADMAP.md
└── MOTIF_ANALYSIS_v2.md

scripts/
└── build_motif_bundles_v3_form_split.py
```

Plus this `SESSION_SUMMARY.md`.

To put them in place on the server:

```bash
cd ~/Projects_shared/CarboDB_v3

# Copy docs
cp /path/to/handoff_outputs/docs/*.md docs/

# Copy script
cp /path/to/handoff_outputs/scripts/build_motif_bundles_v3_form_split.py scripts/

# Generate the v3 bundle
python scripts/build_motif_bundles_v3_form_split.py
tar czf data/motifs_v3_form_split.tar.gz -C data motifs_v3_form_split

# Commit
git add docs/HANDOFF.md docs/ARCHITECTURE.md docs/DATA_INGESTION.md \
        docs/WEBAPP.md docs/DEPLOYMENT.md docs/ROADMAP.md \
        docs/MOTIF_ANALYSIS_v2.md \
        scripts/build_motif_bundles_v3_form_split.py
git commit -m "docs: handoff bundle for shared-owner colleague + v3 motif bundle generator"
git push
```

Send the colleague: motifs_v2_clean.tar.gz (already exists) plus
motifs_v3_form_split.tar.gz (after running the script) plus a link to
docs/MOTIF_ANALYSIS_v2.md.
