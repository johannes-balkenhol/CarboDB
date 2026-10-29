# CarboDB v5.2 — run plan (the upgraded rerun of what v5 did end to end)

*Lives at `CarboDB/docs/RUNPLAN_v5.2.md`. Companion to ROADMAP.md (why/when) and
RESULTS_CATALOG.md (which figure needs what). This file says **what runs, in which order,
with which change versus v5, and how we know it worked**. Nothing outside this plan is run
for Paper 1; changes need a line in the roadmap decisions log.*

Updated 2026-10-29.

---

## 0. Principle: a controlled upgrade, not a new pipeline

v5 ran end to end (March–April 2026, scripts 01–24). v5.2 is the **same pipeline, same
inputs, same stages**, with a closed list of improvements. To know which change moved which
number, the improvements are applied as a **ladder**, each step benchmarked:

| Step | What changes | Expected effect | Why it matters for the paper |
|---|---|---|---|
| L0 | v5 features + v5 labels through the **new** harness | numbers ≈ v5 (±0.01) | proves the new code reproduces the old results |
| L1 | + Km cleaning (use/convert/exclude) | Km R² may drop (honest) | the real data-quality effect |
| L2 | + grouped evaluation (cluster + measurement group) | Km R² drops toward honest CV (~0.7) | the number reviewers will trust |
| L3 | + fixes (894 rows, EC on all positives) | small | completeness |
| L4 | + new features (Pfam bit scores, InterPro full, nearest-measured Km, kingdom, BLAST) | Km within-EC up? | the improvement, block by block (ablation) |
| L5 | + RuBisCO features (Form, Q149, conservation score) in the RuBisCO Km model | small R², big interpretability | Fig 5 |

**v5.2 final = the best ladder step per target**, chosen by the ablation rule:
*every block must improve its target on the grouped test, or it is out.*

---

## 0b. House rules (apply to every stage)

**Hardware: never more than half the machine**, even if a job then takes longer.
| Resource | Limit | How |
|---|---|---|
| CPU | ≤ 64 threads in total across all running jobs | `NJOBS × threads ≤ 64`; jobs queue, they never stack on top of each other |
| GPU (A100 80 GB) | ≤ 40 GB memory | batch size chosen to fit; one GPU job at a time |
| RAM | ≤ 1 TB | big matrices as memory-mapped `.npy` |
| Priority | `nice -n 10` | the other half stays responsive for colleagues |
| Running | `nohup … &` + PID file + log in `CarboDB/docs/` | resumable: every job skips finished chunks |
| Scratch | local SSD `/var/tmp/job37yv/` | results copied to `/storage` as soon as each chunk is done; scratch is never the only copy |

**Storage: where every kind of data lives** (follows the group standard raw / primary / secondary)
```
CarboDB_data/
├── raw/          inputs, never changed (read-only, checksums, SOURCE.md per folder)
│   ├── brenda/  uniprot/  gmgc/              downloads, frozen
│   ├── dbs/pfam/  dbs/interpro/              tool databases with version
│   └── labels/motifs_input_v1/  labels/km_rules_v1/   curated inputs
├── primary/      per-sequence results of tools (regenerable, but expensive)
│   ├── annotations/pfam/domtblout/*.gz       raw tool output, kept
│   ├── annotations/{interpro,diamond,rubisco}/
│   └── embeddings/{esm2,ankh}/*.npy + row_index.tsv   big arrays as files
├── secondary/    derived for analysis (cheap to regenerate)
│   ├── splits/  matrices/  models/  benchmarks/  shap/  predictions/  figures/
└── releases/     what others use
    ├── carbodb_v5.2.sqlite            the database (app, paper tables)
    └── carbodb_v5.2_manifest.json     code commit, input checksums, model hashes
```
**SQL vs files:** everything that is queried row by row (sequences, evidence, annotations with
scores, predictions, benchmark results, TP/TN lists) lives **in SQLite**. Big arrays
(embeddings, feature matrices, models, raw tool output) live **as files**; the database
holds a **pointer** to each one in the `artifacts` table (path, sha256, shape, row order,
run). Nothing big is copied into SQL just to have it there; nothing tabular lives only in a file.

**Never lose data silently:** every stage reports counts in → counts out; any filter
(source, label, quality) is written down with its count. (Example caught 2026-10-29:
222 TrEMBL Km rows would have been dropped by a `source='brenda'` filter.)

## 1. Inputs and environment (frozen)

| Input | Release | Status |
|---|---|---|
| BRENDA | 2026-03-27 download | frozen; checksum into `CarboDB_data/raw/` |
| UniProt SwissProt + TrEMBL + negatives | March 2026 | frozen |
| Pfam-A | as used in v5 (copied 2026-10-29) | ✅ `raw/dbs/pfam` + sha256 |
| InterProScan | 5.72-103.0 | frozen |
| ESM-2 | esm2_t33_650M_UR50D | frozen |
| GMGC10 marine | downloaded May 2026, 500 K Pfam candidates | frozen |
| Environment | `carboxylase` conda env | export to `env/environment.lock.yml` |

No re-download for Paper 1 (it would change every number and cost 1–2 weeks of annotation).

---

## 2. Stages (mirrors v5 scripts; ✅ done · 🔄 running · ⬜ to do)

### S1 — Sequences (v5: 01–03)
- **Keep:** 2,380,446 sequences, CDB IDs, label 0/1/2, EC from sources.
- **Change:** QC rule for the 78 label-1 sequences with non-standard letters (X/U/B…):
  kept in DB, flagged `seq_valid=0`, excluded from training/benchmark.
- **Check:** counts per label identical to v5 except the flagged ones.
- Compute: minutes. ⬜

### S2 — Km labels (v5: 01 BRENDA Km → km_evidence) — **decision gate G1**
- **Sources: only BRENDA Km is used** (decision 2026-10-29). The 222 rows labelled `trembl` in v5 are checked: values that trace back to a BRENDA record (same EC + organism + value) count as BRENDA values transferred to another sequence (`match_level = organism_ec`); all others are kept visible but excluded (`origin_unverified`).
- **Change (F1):** every Km row gets `condition_class` use / convert / exclude, `co2_form`,
  `is_mutant`, `exclude_reason`, `measurement_group_id`, `match_level` (UniProt-exact vs
  organism-broadcast). Rules drafted from the 2,971 rows, agreed with Aparna.
- **Check:** `v_km_training` row count + number of distinct measurement groups per EC;
  no 18 mM Helianthus block; every exclusion has a reason.
- Owner: Johannes + Aparna. Compute: minutes. ⬜ (needs the Km export)

### S3 — Annotation (v5: 04a–04f)
| Block | v5 | v5.2 | Status |
|---|---|---|---|
| Composition / physchem (489) | ✅ | reuse; recompute 894 | ⬜ minutes |
| ESM-2 (1,280) | ✅ | reuse; recompute 894 | ⬜ minutes (A100) |
| InterProScan | ✅ (5 counts used) | reuse; **expand to ~120 binary** (PANTHER, Gene3D, HAMAP…) | ⬜ extract only |
| Pfam | 19 binary flags | **rescan with scores** → `domain_hits`; new features: **best bit score per carboxylase family, normalised by model length** (= conservation to the family) | 🔄 pass 1 Fri, pass 2 Sat |
| Kingdom / taxon | 4 one-hot (Km only) | NCBI lineage for all (reuse `fetch_lineage.py`) | ⬜ hours |
| BLAST/DIAMOND | dropped in v5 | **DIAMOND vs training positives** (pident, evalue, best EC) | ⬜ hours |
| Nearest measured Km | — | Km + identity of the most similar sequence **with experimental Km** | ⬜ computed per CV fold (see S5) |
| RuBisCO block | — | **Form I/II/III** (from HAMAP hits already in InterPro), **Q149** (hmmalign to Form I profile), optional conservation score | ⬜ < 1 h |
| Motifs (Aparna) | — | NRGEIA presence (biotin carboxylases), Q149 above | ⬜ after her handover |
| Ankh | ✅ (unused) | stays unused; optional one benchmark row | — |
- **Check:** 0 sequences without composition/ESM-2/Pfam rows; `run.json` per block.

### S4 — Splits (v5: 06)
- **Keep:** CD-HIT 40 % clusters (45,915), 80/10/10, stratified by EC + Km bin + kingdom.
- **Change:** Km split **grouped by measurement group and cluster** (no copy of the same
  literature value in train and test); split stored in `sequences.split`.
- **Check:** zero measurement groups and zero clusters shared across train/test.

### S5 — Feature matrices (v5: 07)
- **Change:** a manifest `config/features.yaml` lists every block → one matrix per target
  and ladder step; feature order saved with each model.
- **Leakage rules:** ESM-2 PCA (100 comp.) fit on train only; nearest-measured Km and the
  conservation score computed **inside each CV fold** from train sequences only.
- Compute: hours (1 TB RAM, on SSD).

### S6 — Training (v5: 08)
| Target | Model(s) | Notes |
|---|---|---|
| T1 carboxylase yes/no | XGBoost full + ablations (no-ESM, Pfam+composition, ESM-only, motifs-only) + RF + MLP | ablation models = app modes (Aman) |
| T2 EC class | XGBoost on **all 503,841 positives** (v5: tiers 1+2) | 27 classes |
| T3a Km global | XGBoost, log10 Km | grouped CV |
| T3b **Km within EC** | one model per trainable EC (≥ N measurement groups, N decided in G1) | the "good result" from v5; RuBisCO model gets Form + Q149 + conservation |
- Compute: ~1–2 days A100 + CPU. Output: `models` table + model files + feature-order files.

### S7 — Benchmark (v5: 09, 13) — **gate G2: L0 reproduces v5 first**
- Components alone: Pfam rule, PROSITE/InterPro, PANTHER, BLAST-NN (EC + Km transfer), ESM-2 only, motifs only.
- Full pipeline and ablations; learners XGB vs RF vs MLP.
- External Km: UniKP, CatPred (GreenAse if access).
- Km global vs within-EC, per EC.
- Strata: identity bins, fragments 75/50/25 % (regenerate missing BLAST queries), kingdom, no-Pfam, OOD.
- **TP/TN lists** per method → `benchmark_calls`.
- Compute: 2–4 days incl. CatPred setup.

### S8 — Interpretation (v5: 12)
- SHAP per tier (binary, EC, Km global, Km per EC), aggregated per feature block.
- RuBisCO: SHAP of Form / Q149 / conservation in the RuBisCO Km model → Fig 5.

### S9 — Predict + release DB (v5: 10, 05, gmgc/05b)
- **GMGC embeddings were not saved in v5** → ESM-2 for the ~500 K candidates (~15 h A100),
  queued after the Pfam rescan so the 64-thread limit holds; stored in `primary/embeddings/esm2_gmgc/`.
- Predict 2.38 M (EC/Km only for predicted carboxylases; Km within-EC where available,
  else global; `km_in_domain` flag) + GMGC 500 K candidates.
- Build `carbodb_v5.2.sqlite` in **one pass** on SSD (schema v2.3) → `releases/` + manifest
  (code commit, input checksums, model hashes). Integrity checks = 0 missing, 0 FK errors.
- Compute: ~1 day.

### S10 — Figures + manuscript (RESULTS_CATALOG)
- One script per figure, reading only the release DB → Fig 1–6 + supplement tables.

---

## 3. Gates (stop and decide)

| Gate | Question | Who | Before |
|---|---|---|---|
| G1 | Km cleaning rules + minimum measurement groups per EC for a within-EC model | Johannes + Aparna | S2, S6 |
| G2 | Does L0 reproduce v5 within tolerance? If not, find out why before going on | Johannes | S7 ladder |
| G3 | Which blocks stay (ablation rule) | Johannes, data decides | final models |
| G4 | Release integrity (0 missing rows, 0 FK errors, counts match) | Johannes | S9 → Aman |

---

## 4. Calendar (compute in background, half the machine: 64 CPU, 1 TB RAM, A100, SSD)

| Week | Stages | Gate |
|---|---|---|
| 28 Oct – 1 Nov | S0 freeze ✅, Pfam rescan 🔄, Stage A DB for Aman | — |
| 2 – 8 Nov | S1, **S2 (Km rules)**, S3 remaining blocks | G1 |
| 9 – 15 Nov | S4, S5, L0 reproduction | G2 |
| 16 – 29 Nov | S6 training, ladder L1–L5, Aparna handover (Q149 re-test on cleaned Km) | G3 |
| 30 Nov – 9 Dec | S7 full benchmark incl. CatPred, S8 SHAP | — |
| 7 – 14 Dec | S9 predict + release DB → Aman | G4 |
| 14 – 18 Dec | S10 figures, manuscript to co-authors | — |
