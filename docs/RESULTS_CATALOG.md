# CarboDB results catalog — every analysis has one slot and one status

*Lives at `CarboDB/docs/RESULTS_CATALOG.md`. Built 2026-10-29 from the folder inventory
(scripts 09–24, `scripts/gmgc/01–16`, `figures/`, Aparna's motif work) and the paper scope
agreed on 24 June 2026. Change a slot only with a line in the roadmap decisions log.*

**Status codes**
- **RERUN**: depends on models/Km/predictions; regenerated from the v5.2 release DB, numbers change
- **KEEP**: independent of v5.2 (or only cosmetic); re-plot from stored data, numbers stay
- **P2**: moves to Paper 2; frozen as is, not rerun now
- **DROP**: superseded or not defensible; archived, not cited

**Rule for every figure:** one figure = one script = reads only `carbodb_v5.2.sqlite` (+ named
inputs) = writes `CarboDB_data/processed/figures/figN_*`. Then a rerun regenerates every figure.

---

## Paper 1 storyline (24 June): tiered recognition + resource + two use cases

| Slot | Question | Figure / table |
|---|---|---|
| **A** Resource | What is in CarboDB, how is it built? | Fig 1 pipeline + DB overview |
| **B** Tier 1 | Is it a carboxylase? | Fig 2a binary benchmark |
| **B** Tier 2 | Which EC class? | Fig 2b EC benchmark |
| **B** Tier 3 | What affinity (Km), globally and within EC? | Fig 2c Km benchmark |
| **B** Robustness | Low identity, fragments, kingdoms, no-Pfam | Fig 3 strata |
| **C** Interpretation | Which features drive each tier? | Fig 4 SHAP |
| **D** Use case 1 | RuBisCO: what determines Km within the class? | Fig 5 Gln149 |
| **E** Use case 2 | What does a marine metagenome look like? | Fig 6 GMGC scan |
| **F** Access | How do users query it? | Fig S webapp |
| **S** Supplement | TP/TN lists, per-EC tables, literature check, motif search | Tables S1–Sn |

---

## Catalog

### A — Resource and pipeline
| Analysis / script | Existing output | Slot | Status |
|---|---|---|---|
| 01–03 download, merge, CDB IDs | master.tsv, DB counts | Fig 1 numbers | KEEP (inputs frozen) — counts updated after Km cleaning |
| 04a–04e feature layers | feature tables | Fig 1 | KEEP; F2 Pfam rescan adds scores |
| 05 build DB | carbodb.sqlite v2.2 | Fig 1 | RERUN as one-pass v2.3 build |
| Pipeline overview figure | `2026-03-26_carbodb_pipeline_overview.svg` | Fig 1 | KEEP, redraw for v5.2 |

### B — Tiered benchmark
| Analysis / script | Existing output | Slot | Status |
|---|---|---|---|
| 09_benchmark | benchmark_report_v5.json | Fig 2a–c | RERUN |
| 13_publication_benchmark | `benchmark_figure.pdf` | Fig 2 | RERUN |
| 13_benchmark_identity | identity strata | Fig 3 | RERUN |
| 13_hard_tests (fragments, OOD) | fragment/OOD tables | Fig 3 | RERUN (+ regenerate missing BLAST 50/25 % queries) |
| BLAST-NN, Pfam-rule, PANTHER baselines | in 09/13 | Fig 2 rows | RERUN |
| UniKP comparison | in 13 | Fig 2c | RERUN |
| CatPred comparison | — (not done) | Fig 2c | NEW |
| RF vs MLP vs XGBoost (Simone) | partly | Fig 2 / S | NEW |
| Km global vs within-EC | per-EC R² table | Fig 2c | RERUN, extended |
| TP/TN hit lists per approach (Thomas) | — | Table S | NEW, from `benchmark_calls` |
| gmgc/16b literature benchmark v2 | 4/5 within 3× | Table S | RERUN (Km values change) |
| GreenAse external validation | `2026-06-01_greenase_A` (access failed) | Table S | NEW if access, else one sentence |

### C — Interpretation
| Analysis / script | Existing output | Slot | Status |
|---|---|---|---|
| 12_shap_analysis | shap_*.npz | Fig 4 | RERUN |
| 16_figure_feature_importance | `feature_importance_figure.pdf` | Fig 4 | RERUN |
| Local SHAP per sequence (Aman, app) | in app | Fig S | KEEP (app feature) |

### D — Use case 1: RuBisCO Km determinant
| Analysis / script | Existing output | Slot | Status |
|---|---|---|---|
| Form I/II/III classification (HAMAP) | motif bundles v3.1 | Fig 5 | KEEP |
| Gln149 within-clade test (Aparna) | Section 2.5, Figs 6–7 draft | Fig 5 | RERUN on cleaned Km; settle ρ=−0.34 vs −0.27 |
| Gln149 as feature in Km-within-RuBisCO model | — | Fig 5 / 2c | NEW |
| RF leave-genus-out AUC 0.77, SHAP flat, no epistasis | Section 2.5 | Fig 5 text | RERUN after Km cleaning |
| NRGEIA biotin-carboxylase marker | Aparna MEME | Fig 5 text / S | KEEP as family marker |
| fpocket channel classification | channel_analysis/ | — | P2 |
| gmgc/13–13c RuBisCO forms in GMGC | forms tables | Fig 6 | RERUN |

### E — Use case 2: GMGC marine scan
| Analysis / script | Existing output | Slot | Status |
|---|---|---|---|
| gmgc/01–04b translate, prescreen, filter | 500 K candidates | Fig 6 | KEEP (inputs unchanged) |
| gmgc/05b v5 cascade | 171,735 high-conf | Fig 6 | RERUN with v5.2 models |
| gmgc/06b metadata join, 08 novelty | tables | Fig 6 | RERUN |
| gmgc/09b per-EC/genus, 10 sample PCA | tables | Fig 6 / S | RERUN |
| gmgc/14 RuBisCO geography, 15b ocean basins | maps | Fig 6 | RERUN |
| gmgc/11 Km random forest, 12 bitscore RF | exploratory | — | DROP (superseded by v5.2 Km) |
| 22/23 TARA scan + figure | `tara_metagenome_figure.pdf` | — | P2 (GMGC is the Paper 1 metagenome) |

### F — Access
| Analysis / script | Existing output | Slot | Status |
|---|---|---|---|
| 11_annotate_sequence + webapp | CarboDB-App | Fig S | KEEP; app reads release DB |
| Webapp screenshots | `2026-05-08_webapp-screenshots.docx` | Fig S | RERUN (new screenshots on v5.2) |

### Paper 2 (frozen now)
| Analysis / script | Existing output | Status |
|---|---|---|
| 14 biological analysis, 17 biological figure | `biological_findings_figure.pdf` | P2 |
| 18 carnivorous plants (r=0.512) | `carnivorous_plant_figure.png` | P2 |
| 19 taxonomic deep dive | `taxonomy_figure.pdf` | P2 |
| 20 evolutionary analysis | `evolution_figure.pdf` | P2 |
| 21 ecological findings, explore_ecological | `ecological_findings_figure.pdf` | P2 |
| MEME 5-family discovery (65 motifs), controlled rerun | Aparna | P2 (or Table S if controlled run is ready by 20 Nov) |
| Gln149 MD / mutagenesis plan | doc3 | P2 |
| pLM motif arm (MotifAE / FiGS-MoD) | — | P2 |

### Dropped / superseded
| Item | Why |
|---|---|
| "pfam" and "composite" app modes | not real models (Aman); replaced by ablation models F6 |
| Motif bundles v1 (April) | contaminated by 18 mM HCO₃⁻ block and stuck Km values |
| Positions 150/229/216 for RuBisCO | off-by-one / not replicated; replaced by Gln149 |
| FIMO cross-family "3 motifs in all families" | refuted by BLAST-nr |
| "65 motifs" as a count | parameter ceiling, not a finding |
| v6 Ankh model, schema v3 long tables | v6 backlog, after Paper 1 |

---

## Order of work (follows from the catalog)
1. Release DB v5.2 (pipeline rerun) → everything marked RERUN depends on it.
2. Benchmark harness writes `benchmark_results` + `benchmark_calls` → Fig 2, 3, Table S (TP/TN).
3. SHAP → Fig 4.
4. Gln149 on cleaned Km (Aparna) + Km-within-RuBisCO model → Fig 5.
5. GMGC rerun with v5.2 models → Fig 6.
6. Fig 1 redraw + webapp screenshots last.

Each figure script lives in `CarboDB/scripts/figures/figN_<name>.py`; the old figure scripts
(15–18, 21, 23, gmgc/07) are the starting point for the ones marked RERUN.
