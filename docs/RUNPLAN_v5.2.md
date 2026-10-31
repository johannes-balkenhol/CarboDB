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
| Motifs (Aparna) | — | see **S3m** below | ⬜ after her handover |
| Ankh | ✅ (unused) | stays unused; optional one benchmark row | — |
- **Check:** 0 sequences without composition/ESM-2/Pfam rows; `run.json` per block.

### S3m — Motifs (Aparna's work, integrated as data, not as a side analysis)
Question (decision 2026-10-29): **which motifs associate with Km globally, and with Km within an EC class?**
Not: which motifs discriminate carboxylases (that is answered by Pfam/ESM already).
1. **Motif set frozen** → `motif_catalog` (id, family, consensus/PWM file, origin = MEME run / literature, version).
   Inputs: `raw/labels/motifs_input_v1/` (RuBisCO, ACC, PCC, PEPC 150 each, VKGC 80; gold 116/100/50).
   Set = Aparna's motifs + known markers (NRGEIA biotin carboxylase, RuBisCO Gln149 position).
2. **Scan all carboxylase sequences** (FIMO, fixed p-threshold, same for all) → `motif_hits`
   (sequence_id, motif_id, start, p, score); Gln149 via hmmalign to the Form I profile → residue per sequence.
   Compute: < 1 h at 32 threads.
3. **Association test** on cleaned Km (S2, BRENDA only), unit = measurement group:
   - global: motif present/score vs log10 Km;
   - within EC: same, per EC with ≥ N groups (N from G1); within RuBisCO also within Form I/II/III;
   - controls: phylogeny (cluster/genus-grouped permutation or leave-genus-out), FDR over all motifs.
   Output: `secondary/motifs/km_association.tsv` → Table S + Fig 5.
4. **As features:** motifs that pass step 3 enter the Km models as one block (S6), kept only if the
   ablation improves within-EC Km (G3). "Motifs only" is also one benchmark row (S7).
Owner: Aparna (motif set, step 3 interpretation) + Johannes (scan, DB, models). Gln149 re-test on
cleaned Km settles ρ = −0.34 vs −0.27.

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
- External Km: UniKP, CatPred (GreenAse if access) — see **S7b**.
- Km global vs within-EC, per EC.
- Strata: identity bins, fragments 75/50/25 % (regenerate missing BLAST queries), kingdom, no-Pfam, OOD.
- **TP/TN lists** per method → `benchmark_calls`.
- Compute: 2–4 days incl. CatPred setup.

### S7b — Km comparison with CatPred (and UniKP) on a fair test set
**Why it needs care:** CatPred is trained on CatPred-DB (41,174 Km values from BRENDA +
SABIO-RK). Our Km labels are BRENDA too, so a naive test set contains measurements CatPred has
already seen, which would make CatPred look better than it is. UniKP has the same issue.

**Test sets (all report the same metrics, side by side):**
| Set | Definition | What it shows |
|---|---|---|
| T-all | our grouped Km test split (cleaned, BRENDA only, use/convert) | headline numbers; CatPred may have seen part of it → marked "possibly seen" |
| **T-fair** (main) | T-all minus every sequence ≥ 40 % identical (DIAMOND) to any **CatPred-DB Km training** sequence, and minus measurements present in CatPred-DB (same UniProt + same value) | unseen by **both** models → the comparison for the paper |
| T-identity | T-all binned by identity to the nearest training sequence of each tool (<30 / 30–50 / 50–70 / >70 %) | how each tool degrades for novel enzymes |
| T-per-EC | T-fair per EC with ≥ 10 measurement groups (RuBisCO, CA, PEPC, …) | where the within-EC model helps |

**Inputs for CatPred:** sequence + substrate SMILES. CO₂ → `O=C=O`; values converted from HCO₃⁻
are compared as CO₂. Both tools predict log10 Km; units checked against CatPred-DB (mM).

**Metrics:** R², Pearson, Spearman, fraction within 2× and 10×, MAE in log10; bootstrap
95 % CI **over measurement groups** (not rows). CatPred's uncertainty estimate is reported as an
extra column (are its confident predictions better?).

**Expected story (to test, not assume):** a general predictor (CatPred) vs a carboxylase-
specific one (CarboDB global and within-EC). If CatPred matches us on raw R², the paper's case
rests on EC detection, the curated resource and interpretability (Gln149); if the within-EC
model is better on T-fair per EC, that is the headline.

**Practical:** CatPred code + pretrained checkpoints (github: maranasgroup/catpred) and
CatPred-DB (Hugging Face: RosettaCommons/CatPred-DB) into `raw/external/catpred/` with
version + checksum; run on GPU (few hundred to few thousand sequences: minutes). Setup 1–2 days.
Results → `benchmark_results` (method = catpred / unikp) + `benchmark_calls`.

### S8 — Interpretation (v5: 12)
- SHAP per tier (binary, EC, Km global, Km per EC), aggregated per feature block.
- RuBisCO: SHAP of Form / Q149 / conservation in the RuBisCO Km model → Fig 5.

### S9 — Predict + release DB (v5: 10, 05, gmgc/05b)
- **Sequence sources:** main DB (2,380,446) = BRENDA + SwissProt + TrEMBL (+ negatives).
  Metagenomes are separate: **GMGC** (Paper 1, Fig 6) and **TARA** (Paper 2, frozen, not rerun).
- **GMGC candidates (~500 K, Pfam-prescreened with 19 carboxylase families)** need, like the main DB:
  - **full Pfam scan** with scores (same `04b_pfam_rescan.sh` settings; ~4–5 h at 64 threads);
  - **ESM-2** (not saved in v5; ~15 h A100) → `primary/embeddings/esm2_gmgc/`;
  - composition (minutes) and the other feature blocks of the final models.
  Queued after the main Pfam rescan so the 64-thread limit holds.
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
