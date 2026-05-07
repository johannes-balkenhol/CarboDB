# Data Ingestion Pipeline

How the 50 GB SQLite database was built. Reference for the colleague when
re-running, debugging, or extending any stage.

The canonical reference for stage-by-stage *intent* is the existing
`docs/PIPELINE_PLAN.md` (kept for historical context). This document
adds **what each script actually does, what it produces, and how to
re-run it** without breaking downstream stages.

---

## Quick decision tree

> "I want to..."

- **Add new BRENDA Km data** → re-run script 01 → 03 → 04 (composition only) → 05
- **Add new UniProt sequences** → re-run 02 → 03 → 04 (full) → 05
- **Train a new model on existing features** → 07 → 08 → 09 → 10
- **Re-run feature extraction for one feature type** → re-run only that 04* sub-script
- **Just rebuild the SQLite from existing intermediates** → 05 alone

The pipeline is designed to be re-runnable at any step. Intermediates are
under `data/intermediate/` (TSV/JSONL); the final SQLite at
`data/primary/carbodb.sqlite`.

---

## Stage 0 — bootstrap (00_setup_project.py)

Creates `data/`, `models/`, `figures/`, `tmp/`, `logs/`. Verifies that
`hmmscan`, `interproscan.sh`, `cd-hit`, `blastp`, `makeblastdb` are on
PATH. Run once on any fresh checkout.

---

## Stage 1 — BRENDA download (01_brenda_download.py)

Downloads three things from BRENDA via SOAP:

1. **EC-class list of CO₂-active enzymes.** Discovered dynamically by
   scanning every Km entry in BRENDA for "CO₂" / "carbon dioxide" /
   "HCO3" substrate strings. **Do not hardcode this list** — new BRENDA
   releases add EC classes. Output: `data/intermediate/co2_ec_list.json`,
   currently 39 EC classes.

2. **Gold-standard Km values.** All Km measurements with CO₂/HCO₃⁻ as
   substrate, with organism, sequence (when available), and commentary.
   Output: `data/intermediate/brenda_km_gold.tsv`, ~3,000 unique
   sequences after deduplication.

3. **Negative pool.** Random sample of ~3.3 M BRENDA entries from
   non-CO₂-active EC classes for the binary classifier negatives.

**Credentials required.** Export `BRENDA_EMAIL` and `BRENDA_PASSWORD`
(SHA-256 of password, not plaintext) before running. Without these the
SOAP calls fail with auth error.

**Re-run cost:** ~6 hours on the BRENDA SOAP API. They throttle.

**Known issue:** BRENDA stores some HCO₃⁻ Km measurements with the
substrate field labeled "CO2 in form of HCO3-". Script 01 does not
currently distinguish these from real CO₂ measurements. The May 5
motif bundle cleanup excludes them post-hoc; the v6 retrain (see ROADMAP)
should fix this in script 01 itself.

---

## Stage 2 — UniProt download (02_uniprot_download.py)

REST queries against `https://rest.uniprot.org/uniprotkb/search` for:

- **Reviewed (SwissProt) positives** matching `ec:{co2_ec}` for each
  CO₂-active EC from stage 1. ~50K sequences.
- **Reviewed (SwissProt) negatives** — random sample of
  `reviewed:true AND NOT ec:{co2_ec}` to match positive count.
- **Unreviewed (TrEMBL) positives** — `reviewed:false AND ec:{co2_ec}`.
  ~600K sequences. Where the bulk of the 2.38M comes from.
- **Unreviewed (TrEMBL) negatives** — random sample matched to positives.

Outputs separate FASTAs and metadata TSVs per source under
`data/intermediate/uniprot/`.

**Re-run cost:** ~2-4 hours. UniProt REST is faster than BRENDA SOAP.

---

## Stage 3 — merge & dedupe (03_merge_all_sources.py)

Joins BRENDA + SwissProt + TrEMBL into one master table, deduplicating
by UniProt ID. Assigns:

- `cdb_id` (`CDB000001` … `CDB2380446`) — internal stable ID, separate
  from UniProt accession (UniProt can rename; CDB ID is forever).
- `label` (1 = positive, 0 = negative)
- `source` (`brenda` / `swissprot` / `trembl`)
- `km_best_mM` — picks the most-trustworthy Km value (BRENDA wild-type
  > BRENDA mutant > SwissProt) from `km_evidence`

Output: `data/intermediate/master.tsv` and `data/intermediate/master.fasta`.

**Re-run cost:** minutes.

---

## Stage 4 — feature annotation (04*.py + 04*.sh)

Split into 6 sub-scripts because the feature types have wildly different
runtimes. Each writes its own output TSV/Parquet under
`data/features/<type>/`.

| Script | What | Backend | Time on full DB |
|---|---|---|---|
| `04a_composition.py` | AAC, dipeptide, PAAC, physicochem | numpy | minutes |
| `04b_hmmer.sh` | Pfam HMM scan | HMMER hmmscan | days (parallelized) |
| `04c_interproscan.sh` | InterPro families/superfamilies | InterProScan (Java) | days (parallelized) |
| `04d_ankh.sh` | Ankh embeddings (alternative to ESM-2) | Ankh-base, GPU | days |
| `04e_esm2.py` | ESM-2 t33_650M mean-pooled embedding | facebook ESM-2, GPU | day-ish |
| `04f_blast_benchmark.py` | per-EC BLAST DBs for benchmarking | blastp | hours |

**The `*_wbbi*.sh` scripts** (`run_hmmer_wbbi206.sh`, etc.) are just
wrappers that submit each task to a specific compute node with the
right resource limits. Use these on Artemis cluster, not the bare 04b/c/d
scripts.

**Single-sequence inference** (the webapp /predict path) re-runs HMMER
and InterProScan **per request, in-process**. The DB-precomputed features
are NOT consulted for novel sequences. This is why /predict is slow.

**Re-run cost for one annotation type:**
- composition: ~30 min for full DB
- HMMER: ~24 h with 64 cores
- InterProScan: ~48 h with 64 cores
- ESM-2: ~12 h on a GPU; weeks on CPU
- Ankh: ~12 h on a GPU
- BLAST DB rebuild: ~6 h

---

## Stage 5 — SQLite database build (05_build_database.py)

Joins all intermediates from stages 1–4 into the final 50 GB SQLite at
`data/primary/carbodb.sqlite`. Creates indices, populates evidence tables,
generates `db_metadata`, computes `confidence_scores` and `best_evidence`.

The schema is in `scripts/schema.sql` (29 KB, the canonical source) and
is summarized in `docs/DATABASE.md`.

**Re-run cost:** ~3-4 hours (mostly index creation on the 2.38 M rows).

**Idempotent:** safe to re-run; drops and recreates tables.

**Tables created (16 + 1 cache):**
- `sequences` — master table (id, cdb_id, uniprot_id, ec_number, label,
  source, sequence, length, organism, reviewed, km_best_mM, km_log10_mM,
  seq_valid, created_at)
- `predictions` — model outputs (one row per sequence_id)
- `km_evidence` — every Km measurement (multiple rows per sequence_id)
- `ec_evidence` — every EC assignment with provenance
- `features_composition` — 1 row per seq with all numeric features
- `features_domains` — Pfam + InterPro hits as JSON
- `features_esm2` — 1280-vector blob per sequence
- `features_ankh` — alternative embedding (parallel to ESM-2, currently unused by webapp)
- `features_interpro` — InterPro families/PANTHER/Gene3D/SCOP/CDD/etc.
- `features_blast` — pre-computed BLAST hits per sequence
- `features_expert_motifs` — hand-coded regex motifs
- `features_fimo` — MEME/FIMO motif hits (used by motif analysis)
- `confidence_scores` — derived "trust this prediction" scores
- `best_evidence` — joined view: best Km per sequence
- `id_map` — UniProt ID alternates and synonyms
- `db_metadata` — global key-value (build date, version, etc.)
- `external_annotations_cache` — runtime cache for /external/{uid} (auto-created)

---

## Stage 6 — CD-HIT cluster + train/val/test split (06_cluster_and_split.py)

Clusters at 90% sequence identity using CD-HIT. Assigns each cluster a
single split (train / val / test) so models never see homologs of test
sequences during training. Output: `data/intermediate/splits.tsv`.

**Re-run cost:** ~4 hours for full DB.

---

## Stage 7 — feature matrix prep (07_build_feature_matrix.py)

Joins composition + Pfam + ESM-2 + dipeptide into a 1793-column matrix
in chunks (memory-efficient), with labels. Outputs `train.npz`,
`val.npz`, `test.npz` plus a `feature_names.json`.

**Re-run cost:** ~30 min.

---

## Stage 8 — train models (08_train_models.py)

Three XGBoost trainings in sequence:
1. Binary classifier (positive vs negative carboxylase)
2. EC multiclass (39 classes, only on positives)
3. Km regressor (only on sequences with experimental Km, log10 target)

Outputs: `data/models/binary_v5.json`, `ec_v5.json`, `km_v5.json`,
plus `km_v5_final.json` (with class-weight rebalancing) and
`km_v5_weighted.json` (alternative weighting scheme; experimental).

**Re-run cost:** ~6 hours total on wbbi206.

The training script also writes:
- `webapp/models/feature_names_binary.json`
- `webapp/models/feature_names_km.json`
- `webapp/models/ec_label_map.json`

These are mirrors that the webapp loads at startup. Keep them in sync.

---

## Stage 9 — benchmark (09_benchmark.py + 13_*.py variants)

Multiple benchmark scripts because we kept asking different questions:

- `09_benchmark.py` — basic held-out test metrics (AUROC, F1, R²)
- `13_benchmark_identity.py` — accuracy as a function of nearest-neighbor
  identity in training set (the "OOD curve")
- `13_hard_tests.py` — adversarial cases (mutants, distant homologs)
- `13_publication_benchmark.py` — final figure-ready numbers

Outputs go to `figures/` and `data/benchmarks/`.

---

## Stage 10 — predict everything (10_predict_all.py)

Runs the full feature stack + 3-model cascade across **every sequence
in the database** and writes results into the `predictions` table. This
is what makes /db/seq fast — predictions are precomputed, not live.

**Re-run cost:** ~12 hours on full DB.

**Re-run trigger:** any time models are retrained.

---

## Stage 11 — annotate one sequence (11_annotate_sequence.py)

Helper used by the webapp `pipeline/annotate.py`. Takes a raw sequence,
returns the same payload the webapp uses for live predictions. Useful
for ad-hoc CLI testing.

```bash
python scripts/11_annotate_sequence.py --sequence "MAQ..." --kingdom plant
```

---

## Stages 12-23 — analyses and figures

These produce paper-ready figures and exploratory analyses. Not part of
the daily web app loop.

| Script | Output |
|---|---|
| `12_shap_analysis.py` | per-EC SHAP feature importances (used by Analysis page) |
| `13_*_benchmark.py` | publication metrics |
| `14_biological_analysis.py` | kingdom/phylum-level patterns |
| `15_publication_figure.py` | overview figure (model performance) |
| `16_figure_feature_importance.py` | SHAP figure |
| `17_figure_biological.py` | taxonomic patterns |
| `18_carnivorous_figure.py` | carnivorous-plant case study |
| `19_taxonomic_deepdive.py` | extended taxonomy |
| `20_evolutionary_analysis.py` | phylogenetic patterns |
| `21_ecological_findings_figure.py` | ecological findings |
| `22_metagenome_scan.py` | scan a metagenome catalog with model |
| `23_tara_figure.py` | Tara Oceans figure |

---

## How to re-run safely

```bash
cd ~/Projects_shared/CarboDB_v3

# Always:
conda activate carboxylase

# To re-run a single stage, e.g. composition features only:
python scripts/04a_composition.py --threads 32

# To re-run the entire pipeline from scratch (don't unless really needed):
for s in scripts/0[12345]_*.py scripts/0[6789]_*.py; do
  python "$s" || break
done
```

**Backup before destructive re-runs.** The 50 GB SQLite is hard to rebuild
in one sitting:
```bash
cp data/primary/carbodb.sqlite data/primary/carbodb.sqlite.YYYY-MM-DD.bak
```

**Where things land**:
- Intermediate TSV/JSON: `data/intermediate/`
- Per-feature outputs: `data/features/<type>/`
- Final SQLite: `data/primary/carbodb.sqlite`
- Models: `data/models/` and `webapp/models/`
- Figures: `figures/`
- Logs: project-root `logs/` and `webapp/logs/`

---

## Known issues

1. **BRENDA HCO₃⁻ contamination is not separated at ingest time.** Script 01
   merges CO₂ and HCO₃⁻ Km measurements. The May 5 motif cleanup handled
   this post-hoc; v6 retrain will need this fixed in script 01 with proper
   pH-dependent CO₂↔HCO₃⁻ conversion via Henderson-Hasselbalch. ~20-30 of
   600 RuBisCO Km entries affected. Discuss with Johannes before retrain.

2. **InterProScan sometimes hangs on specific sequences.** The
   `04c_interproscan.sh` wrapper has a 60-min per-batch timeout, but
   individual sequences within a hung batch get retried. Watch the logs.

3. **ESM-2 GPU memory.** With t33_650M and batch_size=8, peaks at ~22 GB
   GPU RAM. wbbi203 has the GPU. Smaller batch_size if OOM.

4. **CD-HIT 4.8.1 is the version in the conda env.** Version 4.6 has
   subtly different cluster representative selection — don't downgrade.

5. **Form I/II/III not labeled in `sequences`.** Forms can be derived from
   `features_interpro.panther_family` (PTHR42704 = Form I large chain) or
   from `raw_ipr_json` text matching. The within-Form-I motif analysis
   (see MOTIF_ANALYSIS_v2.md) uses this.
