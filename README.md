# CarboxyDB

**A large-scale database and machine-learning pipeline for CO₂-fixing carboxylases.**

Predicts whether a protein sequence is a carboxylase, classifies its EC number,
and estimates its Michaelis constant (Km) for CO₂/HCO₃⁻ in mM.

---

## Quick start

```bash
# 1. Clone and enter
git clone <repo_url> CarboxyDB && cd CarboxyDB

# 2. Install Python dependencies
pip install requests tqdm pandas numpy biopython xgboost scikit-learn zeep

# 3. Set BRENDA credentials (needed for script 01)
export BRENDA_EMAIL="your@email.com"
export BRENDA_PASSWORD="sha256_of_your_password"

# 4. Create folders and check tools
python 00_setup_project.py

# 5. Run the pipeline (one script per step)
python scripts/01_brenda_download.py
python scripts/02_uniprot_download.py
python scripts/03_merge_all_sources.py
python scripts/04_annotate_features.py
python scripts/05_build_database.py
python scripts/06_cdhit_cluster.py
python scripts/07_prepare_ml_data.py
python scripts/08_train_models.py
python scripts/09_benchmark.py
python scripts/10_predict_database.py
```

---

## Project layout

```
CarboxyDB/
│
├── config.py                   ← Central config (paths, thresholds, params)
├── 00_setup_project.py         ← Run once: create dirs + check tools
│
├── scripts/
│   ├── 01_brenda_download.py   ← BRENDA SOAP → positives, negatives, Km gold standard
│   ├── 02_uniprot_download.py  ← UniProt REST → SwissProt + TrEMBL + negatives
│   ├── 03_merge_all_sources.py ← Merge, deduplicate, assign CDB_IDs
│   ├── 04_annotate_features.py ← Extract all feature layers (A–E)
│   ├── 05_build_database.py    ← Build carbodb.sqlite from master.tsv + features
│   ├── 06_cdhit_cluster.py     ← CD-HIT 90% clustering → cluster-aware splits
│   ├── 07_prepare_ml_data.py   ← Train/val/test split → data/ml/
│   ├── 08_train_models.py      ← XGBoost: binary + EC class + Km regression
│   ├── 09_benchmark.py         ← Benchmark vs BLAST and Pfam
│   └── 10_predict_database.py  ← Run models on all DB sequences + SHAP
│
├── data/
│   ├── raw/                    ← Downloaded files (not in git)
│   │   ├── brenda/
│   │   └── uniprot/
│   │       ├── swissprot/
│   │       ├── trembl/
│   │       └── negatives/
│   ├── interim/                ← Per-step intermediates
│   ├── primary/                ← master.tsv, master.fasta, id_map.tsv (source of truth)
│   ├── features/               ← Feature TSVs, one per layer
│   │   ├── composition/
│   │   ├── domains/
│   │   ├── motifs/
│   │   ├── blast/
│   │   ├── esm2/
│   │   └── meme/               ← PENDING (MEME subproject)
│   ├── ml/                     ← train/val/test splits
│   ├── benchmark/
│   ├── shap/
│   └── dbs/                    ← Pfam HMM, BLAST db, PROSITE (not in git)
│
├── database/
│   └── carbodb.sqlite
│
├── models/
│   ├── carboxy_binary_v3.pkl
│   ├── carboxy_ec_class_v3.pkl
│   ├── carboxy_km_v3.pkl
│   └── *.json                  ← Model metadata
│
└── logs/                       ← One log file per script run
```

---

## ID system

Every sequence gets **two IDs** assigned in script 03:

| Column | Example | Description |
|---|---|---|
| `cdb_id` | `CDB000001` | Internal primary key. Never changes. |
| `uniprot_id` | `P00187` | UniProt accession. Used as foreign key. May be absent for non-UniProt sequences. |

All feature tables and the SQLite database use `cdb_id` as the primary key.
`uniprot_id` has a unique index for fast lookup.

---

## Data sources and scale

| Source | Positives | Negatives | Notes |
|---|---|---|---|
| BRENDA (script 01) | ~695,000 | ~3,300,000 | Experimental Km for ~3,001 sequences |
| SwissProt (script 02) | ~50,000 | ~50,000 | Reviewed UniProt entries |
| TrEMBL (script 02) | ~600,000 | ~600,000 | Unreviewed, capped per EC |
| **Total (deduplicated)** | **~1,345,000** | **~3,975,000** | |
| ML training subset | ~100,000 | ~500,000 | CD-HIT 90% cluster split |

---

## Feature layers

| Layer | Name | Features | Tool | Status |
|---|---|---|---|---|
| A1 | AA composition | 20 | BioPython | ✓ |
| A2 | Dipeptide frequency | 400 | BioPython | ✓ |
| A3 | Pseudo-AAC | 30 | custom | ✓ |
| A4 | Physicochemical | ~20 | BioPython | ✓ |
| A5 | Catalytic core | ~17 | custom | ✓ |
| A6 | EC-specific motifs | 7 | regex | ✓ |
| B1 | Pfam domains | ~30 | HMMER3 | ✓ |
| B2 | PROSITE patterns | 14 | regex | ✓ |
| C | BLAST homology | 4 | BLAST+ | ✓ |
| D | MEME motifs | 65 | FIMO | ⏳ PENDING |
| E | ESM-2 embeddings | 1280 | GPU/HPC | optional |

**v3 model** = layers A+B+C (~523 features, R² Km = 0.91)  
**v5 model** = layers A+B+C+E (~1803 features, R² Km = 0.92+)

---

## Evidence tiers

| Tier | Label | Source |
|---|---|---|
| 1 | experimental | BRENDA Km measured |
| 2 | curated | SwissProt manually reviewed |
| 3 | predicted | TrEMBL / model output |
| 4 | inferred | BLAST best-hit / Pfam |

---

## External tools required

| Tool | Version | Install |
|---|---|---|
| HMMER3 | ≥ 3.3 | `conda install -c bioconda hmmer` |
| BLAST+ | ≥ 2.12 | `conda install -c bioconda blast` |
| CD-HIT | ≥ 4.8 | `conda install -c bioconda cd-hit` |
| FIMO (MEME) | ≥ 5.5 | `conda install -c bioconda meme` (optional) |

Check all at once: `python 00_setup_project.py`

---

## MEME motif subproject (pending)

When the MEME subproject delivers results, place the hit file at:

```
data/features/meme/meme_hits.tsv
```

Expected format:
```
cdb_id    meme_rubisco_1_GKST    meme_ca_2_HHC    ...
CDB000001    1    0    ...
```

Script 04 will detect this file automatically and merge it into the feature set.

---

## Km values

All Km values are stored in **millimolar (mM)**.  
The `km_log10_mM` column = log₁₀(Km_mM) and is used as the regression target.  
Gold-standard range from BRENDA: 0.0008 – 83.0 mM.

---

## Reproducibility

Each script writes a timestamped log to `logs/`.  
Raw data files are named with a timestamp: `brenda_positives_20250101_120000.tsv`.  
The `data/primary/master.tsv` is the single source of truth — regenerate it by
re-running scripts 01–03 in order.

## Model performance (v5, April 7 2026)
| Task | Metric | Value |
|------|--------|-------|
| Binary carboxylase detection | AUROC | 0.9999 |
| Binary carboxylase detection | AUPRC | 0.9997 |
| Binary carboxylase detection | F1 | 0.9965 |
| EC class prediction | pending | — |
| Km regression | pending | — |

## Model performance update (v5 final, April 8 2026)
| Task | Metric | Value | Notes |
|------|--------|-------|-------|
| Binary detection | AUROC | 0.9999 | cluster-based test set |
| Binary detection | AUPRC | 0.9997 | |
| Binary detection | F1 | 0.9965 | |
| EC class (26/40) | Top-1 acc | 0.9982 | BRENDA+SwissProt only |
| EC class (26/40) | Top-3 acc | 0.9996 | |
| EC class (26/40) | F1 macro | 0.9389 | |
| Km regression | R² | 0.9253 | random split, EC+kingdom features |
| Km regression | Pearson r | 0.9628 | log10 mM scale |
| Km regression | RMSE | 0.4072 | log10 mM |
