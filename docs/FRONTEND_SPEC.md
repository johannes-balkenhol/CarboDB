# CarboDB-App Frontend Specification

**Stack:** Flask (Python) backend + Vue.js frontend (existing structure in this repo)  
**Database:** `carbodb.sqlite` — 2,380,446 sequences, 510,454 predicted carboxylases  
**Annotation API:** `scripts/11_annotate_sequence.py` — see `docs/API.md`  

This document defines every page, its data requirements, UI components, and API fields. A frontend developer should be able to implement all pages from this spec alone.

---

## Architecture Overview

```
CarboDB-App/
├── backend/          ← Flask API (Python)
│   ├── routes.py     ← All API endpoints defined here
│   └── ...
├── frontend/         ← Vue.js SPA
│   └── ...
└── data/
    └── carbodb.sqlite  ← Read-only; mounted into backend container
```

The frontend calls backend API endpoints (JSON). The backend queries `carbodb.sqlite` and calls the annotation script.

---

## Pages

### 1. Home `/`

**Purpose:** Landing page — overview of the database and entry points.

**Components:**
- Hero with title, short description, search bar (UniProt ID or organism name)
- 4 stat cards: total sequences, predicted carboxylases, EC classes covered, sequences with Km

**Stat card data — API endpoint:** `GET /api/stats`

```json
{
  "total_sequences": 2380446,
  "predicted_carboxylases": 510454,
  "high_confidence": 485081,
  "ec_classes_covered": 26,
  "sequences_with_km_pred": 352888,
  "sequences_with_km_experimental": 3152,
  "km_range_mM": [0.0017, 114.0]
}
```

- EC class distribution bar chart (top 10 EC classes by count)
- Quick links: Browse database | Annotate new sequence | Experimental data

---

### 2. Browse `/browse`

**Purpose:** Filterable, paginated table of all predicted carboxylases.

**Filters (left sidebar or top bar):**
- EC class (dropdown, multi-select)
- Confidence: high / medium / low (checkboxes)
- Kingdom: bacteria / plant / archaea / fungi (checkboxes)  
- Km range: min/max slider (mM, log scale)
- Has experimental Km: yes/no toggle
- Reviewed (SwissProt): yes/no toggle

**Table columns:**
| Column | Source field | Notes |
|---|---|---|
| CDB ID | `sequences.cdb_id` | Link to detail page |
| UniProt | `sequences.uniprot_id` | Link to UniProt |
| Organism | `sequences.organism` | |
| EC (experimental) | `sequences.ec_number` | Grey if null |
| EC (predicted) | `predictions.ec_pred` + `ec_prob` | Show prob as badge |
| Km experimental | `sequences.km_best_mM` | Grey if null |
| Km predicted | `predictions.km_pred_mM` | Grey if null |
| Confidence | `confidence_scores.confidence_label` | Color-coded badge |
| Reviewed | `sequences.reviewed` | ✓/— |

**API endpoint:** `GET /api/browse`

Query parameters:
```
ec_pred=4.1.1.39
confidence=high,medium
km_min=0.01
km_max=10
has_km_exp=true
reviewed=true
kingdom=bacteria
page=1
per_page=50
sort=km_pred_mM
order=asc
```

Response:
```json
{
  "total": 12483,
  "page": 1,
  "per_page": 50,
  "results": [
    {
      "cdb_id": "CDB000042",
      "uniprot_id": "P00875",
      "organism": "Spinacia oleracea",
      "ec_experimental": "4.1.1.39",
      "ec_predicted": "4.1.1.39",
      "ec_prob": 1.0,
      "km_experimental_mM": 0.01,
      "km_predicted_mM": 0.0101,
      "confidence_label": "high",
      "reviewed": true,
      "length": 475
    }
  ]
}
```

---

### 3. Sequence Detail `/sequence/:cdb_id`

**Purpose:** Full detail view for a single sequence. This is the most important page.

**Sections:**

#### 3a. Header
- UniProt ID (linked), organism, sequence length, source (SwissProt/TrEMBL/BRENDA)
- Confidence badge (high/medium/low)
- Links: UniProt page, AlphaFold structure (if available)

#### 3b. Carboxylase Prediction
- Large probability gauge or bar: `co2_prob` 
- `is_co2_pred` label
- `confidence_label`

#### 3c. EC Class Prediction
- Predicted EC class + name (large, prominent)
- Probability bar chart — top 5 EC classes from `ec_probabilities`
- Experimental EC class if available (`sequences.ec_number`)
- Agreement indicator: ✓ predicted matches experimental / ✗ mismatch / — no experimental

**Top SHAP features for EC prediction:**  
Display top 5 feature groups as a horizontal bar chart:
- Data source: pre-computed from `data/shap/shap_ec_per_class.json` for the predicted EC class
- Show: Pfam domains, ESM-2 (%), Dipeptide (%), etc.
- Label each bar with group name and % importance
- Note: these are class-level importances, not per-sequence SHAP

#### 3d. Km Prediction
- Predicted Km: large number (mM) with log₁₀ value
- Experimental Km from BRENDA if available — show both and the difference
- Context: where does this Km fall within the EC class distribution?  
  Show mini histogram of all Km values for this EC class, with a marker for this sequence

**Top SHAP features for Km within this EC class:**  
- Data source: pre-computed from `data/shap/shap_km_per_ec.json` for the predicted EC class
- Show top 5 features that drive Km prediction within this EC class
- Show direction: + pushes Km higher, − pushes Km lower

#### 3e. Domain Architecture
- List of Pfam hits from `features_domains` table
- For each hit: Pfam ID, domain name, description, E-value
- Visual domain diagram (horizontal bar showing domain positions along sequence length)

#### 3f. AlphaFold Structure
- Embed AlphaFold viewer using UniProt ID:  
  `https://alphafold.ebi.ac.uk/entry/{uniprot_id}`  
  Use the public AlphaFold iframe or link to the page
- Show pLDDT score distribution if available from AlphaFold API

#### 3g. Experimental Data (if available)
- Table of all Km values from `km_evidence` where `evidence_tier = 1`
- Columns: km_value_mM, substrate, source (BRENDA), commentary
- Table of EC annotations from `ec_evidence`

#### 3h. Sequence
- Collapsible section showing the full amino acid sequence
- Copy-to-clipboard button

**API endpoint:** `GET /api/sequence/:cdb_id`

```json
{
  "cdb_id": "CDB000042",
  "uniprot_id": "P00875",
  "organism": "Spinacia oleracea",
  "length": 475,
  "source": "brenda",
  "reviewed": true,
  "sequence": "MSPQTETK...",
  "ec_experimental": "4.1.1.39",
  "km_experimental_mM": 0.01,

  "prediction": {
    "is_co2_pred": true,
    "co2_prob": 1.0,
    "confidence_label": "high",
    "ec_pred": "4.1.1.39",
    "ec_name": "ribulose-bisphosphate carboxylase (RuBisCO)",
    "ec_prob": 1.0,
    "ec_top5": {"4.1.1.39": 1.0, "4.2.1.1": 0.0, ...},
    "km_pred_mM": 0.0101,
    "km_pred_log10": -1.996
  },

  "pfam_hits": [
    {"pfam_id": "PF00016", "name": "RuBisCO large subunit", "evalue": 1e-120},
    {"pfam_id": "PF02788", "name": "Biotin carboxylase N-terminal", "evalue": 1e-45}
  ],

  "shap_ec_class": {
    "ec": "4.1.1.39",
    "group_importance": {"Pfam domains": 62.0, "ESM-2 embedding": 18.0, "Dipeptide": 13.0},
    "top_features": [
      {"feature": "pfam_PF00016", "group": "Pfam domains", "mean_abs_shap": 0.42}
    ]
  },

  "shap_km_within_ec": {
    "ec": "4.1.1.39",
    "top_features": [
      {"feature": "dp_DQ", "group": "Dipeptide", "diff_shap": 0.104, "direction": "high_km"},
      {"feature": "dp_VA", "group": "Dipeptide", "diff_shap": -0.015, "direction": "low_km"}
    ]
  },

  "km_evidence": [
    {"km_value_mM": 0.01, "substrate": "CO2", "source": "brenda_experimental", "evidence_tier": 1}
  ]
}
```

---

### 4. Annotate New Sequence `/annotate`

**Purpose:** Submit a new protein sequence for carboxylase annotation.

**Form:**
- Text area: paste FASTA sequence (single sequence)
- Kingdom selector: bacteria / plant / archaea / fungi (default: bacteria)
- Fast mode toggle: "Skip ESM-2 (faster, ~3s, less accurate)"
- Submit button

**Flow:**
1. User submits → POST to `/api/annotate/submit` → returns `task_id`
2. Frontend polls `GET /api/annotate/result/:task_id` every 2s
3. Show progress spinner with estimated time (30s full / 3s fast)
4. On completion: show result in the same format as the sequence detail page (sections 3b–3e)
5. Note: results are NOT saved to the database — this is a live prediction only

**API endpoints:** See `docs/API.md` for full async job queue spec.

**Result display:**  
Same components as sequence detail page (EC prediction, Km prediction, Pfam hits, SHAP), but without the experimental data section and without AlphaFold (no UniProt ID for novel sequences).

---

### 5. Batch Annotation `/batch`

**Purpose:** Upload a multi-FASTA file and download annotation results as TSV/JSON.

**Form:**
- File upload: FASTA file (multi-sequence, max 100 sequences per job)
- Kingdom selector
- Fast mode toggle (recommended for batch)
- Output format: TSV / JSON

**Flow:**
1. Upload file → POST `/api/batch/submit`
2. Backend splits into individual sequences, runs annotation pipeline
3. Poll for completion
4. Download results file

**TSV output columns:**
`query_id, sequence_length, is_carboxylase, co2_prob, confidence, ec_predicted, ec_prob, km_predicted_mM, pfam_hits, warnings`

---

### 6. Experimental Data `/experimental`

**Purpose:** Browse only sequences with experimentally measured CO₂ Km values from BRENDA.

**Table columns:** UniProt, organism, EC class, Km (mM), substrate, BRENDA source, sequence length, predicted Km, error (|log₁₀ predicted − log₁₀ experimental|)

**Filters:** EC class, Km range, organism, reviewed

**API endpoint:** `GET /api/experimental`

```json
{
  "total": 3152,
  "results": [
    {
      "uniprot_id": "P00875",
      "organism": "Spinacia oleracea",
      "ec_number": "4.1.1.39",
      "km_experimental_mM": 0.01,
      "km_predicted_mM": 0.0101,
      "log10_error": 0.004,
      "substrate": "CO2"
    }
  ]
}
```

---

### 7. Statistics `/stats`

**Purpose:** Database overview and model performance.

**Charts:**
- EC class distribution (horizontal bar chart, top 15 classes)
- Km distribution (histogram, log scale, all predicted + experimental)
- Confidence distribution (pie: high/medium/low/review)
- Organism kingdom distribution (pie)
- Scatter plot: experimental vs. predicted Km (log scale, colored by EC class)

**Model performance table:**

| Metric | Value |
|---|---|
| Binary AUROC | 0.9999 |
| EC Top-1 accuracy | 99.82% |
| EC Top-3 accuracy | 99.96% |
| Km R² | 0.9503 |
| Km Pearson r | 0.975 |
| vs. PANTHER EC accuracy | 99.8% vs 94% |
| vs. Pfam rule-based | 99.8% vs 31% |

**API endpoint:** `GET /api/stats` — returns all numbers for the page.

---

## API Endpoint Summary

| Method | Endpoint | Description |
|---|---|---|
| GET | `/api/stats` | Database statistics |
| GET | `/api/browse` | Paginated sequence table with filters |
| GET | `/api/sequence/:cdb_id` | Full detail for one sequence |
| GET | `/api/search?q=` | Search by UniProt ID or organism name |
| GET | `/api/experimental` | BRENDA experimental Km data |
| POST | `/api/annotate/submit` | Submit annotation job |
| GET | `/api/annotate/result/:task_id` | Poll annotation job result |
| POST | `/api/batch/submit` | Submit batch annotation job |
| GET | `/api/batch/result/:task_id` | Poll batch job result |
| GET | `/api/ec_classes` | List of all 26 EC classes with names and counts |

---

## Data Files for SHAP Display

The SHAP feature importance data is pre-computed — the frontend does NOT need to call the ML model for this. Load these JSON files at backend startup and serve them via API:

| File | Used on page | Content |
|---|---|---|
| `data/shap/shap_ec_per_class.json` | Sequence detail §3c | Per-EC top features for EC prediction |
| `data/shap/shap_km_per_ec.json` | Sequence detail §3d | Per-EC within-class Km drivers |
| `data/shap/shap_km_within_ec.json` | Sequence detail §3d | High vs low Km feature comparison |
| `data/shap/shap_binary_global.json` | Stats page | Global binary model feature importance |

---

## Notes for the Frontend Developer

- The database is **read-only** — no writes from the API
- All Km values are in **mM**; show both raw mM and log₁₀ where space allows
- `ec_pred` can be null for non-carboxylase predictions — handle gracefully
- `km_pred_mM` is null for EC classes outside the 10 trainable classes — show "N/A"
- The SHAP data (`shap_ec_per_class.json` etc.) contains class-level importances, not per-sequence — this is expected and should be explained with a tooltip: "Feature importance shown is the average for all [EC class] sequences"
- AlphaFold viewer: use `https://alphafold.ebi.ac.uk/entry/{uniprot_id}` — only available for UniProt sequences; not available for novel sequences submitted via `/annotate`
- Annotation runtime: full pipeline with ESM-2 takes ~30s — always use async with a progress indicator
- ESM-2 model (~2.5 GB) must be cached on the server — do not re-download on each request
