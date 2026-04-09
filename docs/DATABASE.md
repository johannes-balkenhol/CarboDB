# CarboDB Database Documentation

**File:** `data/primary/carbodb.sqlite`  
**Size:** ~30 GB  
**Sequences:** 2,380,446  
**Predicted carboxylases:** 510,454 (21.4%)  
**Experimental Km values:** 352,888 sequences  

---

## Table Overview

| Table | Rows | Description |
|---|---|---|
| `sequences` | 2,380,446 | Master sequence table — one row per protein |
| `predictions` | 2,380,446 | ML predictions for every sequence (v5 models) |
| `confidence_scores` | 2,380,446 | Aggregated confidence per sequence |
| `km_evidence` | ~352,888 | All Km values (experimental + predicted) |
| `ec_evidence` | ~2.38M | All EC annotations (experimental + predicted) |
| `best_evidence` | 2,380,446 | Best available EC + Km per sequence (denormalized) |
| `features_composition` | 2,380,446 | 489-dim composition feature vectors |
| `features_domains` | variable | Pfam/InterPro domain hits |
| `features_esm2` | 2,380,446 | 1,280-dim ESM-2 embeddings |
| `features_ankh` | 2,380,446 | 1,280-dim Ankh embeddings (in progress) |
| `features_interpro` | variable | InterPro annotations |
| `features_blast` | variable | BLAST hits against training set |
| `features_expert_motifs` | variable | Manual carboxylase motif hits |
| `features_fimo` | variable | MEME/FIMO motif scan results |
| `id_map` | 2,380,446 | CDB_ID ↔ UniProt ID mapping |
| `db_metadata` | — | Schema version, build timestamps |

---

## Core Tables

### `sequences`
The master sequence table. Every protein in the database has exactly one row.

| Column | Type | Description |
|---|---|---|
| `id` | INTEGER PK | Internal row ID |
| `cdb_id` | TEXT | CarboDB ID (e.g. `CDB000001`) — stable identifier |
| `uniprot_id` | TEXT | UniProt accession (e.g. `P00875`) |
| `ec_number` | TEXT | Experimentally assigned EC class (NULL if unknown) |
| `label` | INTEGER | Training label: 1=carboxylase, 0=negative, 2=ancestral-related |
| `source` | TEXT | Data source: `brenda`, `swissprot`, `trembl`, `brenda_neg`, `uniprot_neg` |
| `sequence` | TEXT | Full amino acid sequence |
| `length` | INTEGER | Sequence length in amino acids |
| `organism` | TEXT | Organism name from UniProt |
| `reviewed` | INTEGER | 1=SwissProt reviewed, 0=TrEMBL unreviewed |
| `km_best_mM` | REAL | Best experimental Km in mM (NULL if no BRENDA data) |
| `km_log10_mM` | REAL | log₁₀(km_best_mM) |
| `seq_valid` | INTEGER | 1=passes validation filters, 0=rejected |
| `created_at` | TEXT | Timestamp |

**Sources:**
- `brenda`: positive carboxylases with experimental Km from BRENDA
- `swissprot`: SwissProt-reviewed carboxylases without Km data
- `trembl`: TrEMBL predicted carboxylases
- `brenda_neg` / `uniprot_neg`: negative training sequences

---

### `predictions`
ML model predictions for every sequence. The primary table for webapp queries.

| Column | Type | Description |
|---|---|---|
| `id` | INTEGER PK | Internal row ID |
| `sequence_id` | INTEGER | FK → sequences.id |
| `uniprot_id` | TEXT | UniProt accession (denormalized for fast lookup) |
| `model_version` | TEXT | Model version used (e.g. `v5`) |
| `is_co2_pred` | INTEGER | 1=predicted carboxylase, 0=non-carboxylase |
| `co2_prob` | REAL | Carboxylase probability [0–1] |
| `ec_pred` | TEXT | Predicted EC class (e.g. `4.1.1.39`) |
| `ec_prob` | REAL | Confidence in EC prediction [0–1] |
| `km_pred_mM` | REAL | Predicted CO₂ Km in mM (NULL if EC not in trainable set) |
| `km_pred_log10` | REAL | log₁₀(km_pred_mM) |
| `created_at` | TEXT | Timestamp |

**Trainable EC classes for Km prediction:**
`4.2.1.1`, `4.1.1.39`, `4.1.1.31`, `4.1.1.49`, `6.3.4.14`, `4.1.1.32`, `6.4.1.1`, `6.4.1.2`, `6.4.1.3`, `6.4.1.4`

---

### `confidence_scores`
Aggregated confidence assessment per sequence.

| Column | Type | Description |
|---|---|---|
| `sequence_id` | INTEGER | FK → sequences.id |
| `uniprot_id` | TEXT | UniProt accession |
| `method_agreement` | INTEGER | 1=Pfam and ML agree, 0=disagreement |
| `ec_confidence` | REAL | EC prediction confidence score |
| `km_confidence` | REAL | Km prediction confidence score |
| `overall_score` | REAL | Combined confidence [0–1] |
| `confidence_label` | TEXT | `high` (≥0.90), `medium` (0.70–0.90), `low` (0.50–0.70), `review` (<0.50) |

**Distribution in current database:**
- high: 485,081
- medium: 6,977
- low: 18,396
- review: 1,869,992

---

### `km_evidence`
All Km values — both experimental (BRENDA) and predicted (ML model).

| Column | Type | Description |
|---|---|---|
| `sequence_id` | INTEGER | FK → sequences.id |
| `uniprot_id` | TEXT | UniProt accession |
| `ec_number` | TEXT | EC class this Km applies to |
| `km_value_mM` | REAL | Km value in mM |
| `km_log10_mM` | REAL | log₁₀(km_value_mM) |
| `km_unit` | TEXT | Unit (always `mM`) |
| `substrate` | TEXT | Substrate description from BRENDA |
| `source` | TEXT | `brenda_experimental`, `model_v5`, etc. |
| `evidence_tier` | INTEGER | 1=experimental, 2=curated, 3=predicted, 4=inferred |
| `commentary` | TEXT | BRENDA assay notes |
| `model_version` | TEXT | Model version (for predicted values) |

---

### `ec_evidence`
All EC annotations — experimental and predicted.

| Column | Type | Description |
|---|---|---|
| `sequence_id` | INTEGER | FK → sequences.id |
| `uniprot_id` | TEXT | UniProt accession |
| `ec_number` | TEXT | EC class |
| `source` | TEXT | `uniprot`, `brenda`, `model_v5`, `pfam_rule` |
| `evidence_tier` | INTEGER | 1=experimental, 2=curated, 3=predicted |
| `confidence` | REAL | Confidence score [0–1] |
| `model_version` | TEXT | Model version (for predicted values) |

---

## Key Queries for Webapp

### Browse predicted carboxylases with filters
```sql
SELECT 
    s.cdb_id,
    s.uniprot_id,
    s.organism,
    s.ec_number        AS ec_experimental,
    s.km_best_mM       AS km_experimental_mM,
    p.ec_pred,
    p.ec_prob,
    p.km_pred_mM,
    c.confidence_label
FROM sequences s
JOIN predictions p ON p.sequence_id = s.id
JOIN confidence_scores c ON c.sequence_id = s.id
WHERE p.is_co2_pred = 1
  AND c.confidence_label = 'high'
  AND p.ec_pred = '4.1.1.39'         -- filter by EC
  AND p.km_pred_mM < 1.0             -- filter by Km
ORDER BY p.km_pred_mM ASC
LIMIT 100 OFFSET 0;
```

### Lookup by UniProt ID
```sql
SELECT 
    s.*,
    p.is_co2_pred, p.co2_prob, p.ec_pred, p.ec_prob,
    p.km_pred_mM, p.km_pred_log10,
    c.confidence_label, c.overall_score
FROM sequences s
JOIN predictions p ON p.sequence_id = s.id
JOIN confidence_scores c ON c.sequence_id = s.id
WHERE s.uniprot_id = 'P00875';
```

### Top predicted carboxylases by confidence
```sql
SELECT s.uniprot_id, s.organism, p.ec_pred, p.co2_prob, p.km_pred_mM
FROM sequences s
JOIN predictions p ON p.sequence_id = s.id
WHERE p.is_co2_pred = 1
ORDER BY p.co2_prob DESC
LIMIT 50;
```

### EC class distribution
```sql
SELECT ec_pred, COUNT(*) AS n
FROM predictions
WHERE is_co2_pred = 1
GROUP BY ec_pred
ORDER BY n DESC;
```

### Sequences with both experimental and predicted Km
```sql
SELECT 
    s.uniprot_id, s.organism, s.ec_number,
    s.km_best_mM AS km_experimental,
    p.km_pred_mM AS km_predicted,
    ABS(s.km_log10_mM - p.km_pred_log10) AS log10_error
FROM sequences s
JOIN predictions p ON p.sequence_id = s.id
WHERE s.km_best_mM IS NOT NULL
  AND p.km_pred_mM IS NOT NULL
ORDER BY log10_error ASC;
```

### Kingdom breakdown
```sql
SELECT 
    CASE 
        WHEN s.organism LIKE '%sapiens%' OR s.organism LIKE '%musculus%' THEN 'Eukaryota'
        WHEN s.reviewed = 1 THEN 'SwissProt'
        ELSE 'TrEMBL'
    END AS group_label,
    COUNT(*) AS n,
    AVG(p.km_pred_mM) AS mean_km_pred
FROM sequences s
JOIN predictions p ON p.sequence_id = s.id
WHERE p.is_co2_pred = 1
GROUP BY group_label;
```

---

## Recommended Indexes for Webapp Performance

Run once after database is complete:

```sql
CREATE INDEX IF NOT EXISTS idx_pred_is_co2 ON predictions(is_co2_pred);
CREATE INDEX IF NOT EXISTS idx_pred_ec ON predictions(ec_pred);
CREATE INDEX IF NOT EXISTS idx_pred_km ON predictions(km_pred_mM);
CREATE INDEX IF NOT EXISTS idx_pred_prob ON predictions(co2_prob);
CREATE INDEX IF NOT EXISTS idx_seq_uniprot ON sequences(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_seq_organism ON sequences(organism);
CREATE INDEX IF NOT EXISTS idx_conf_label ON confidence_scores(confidence_label);
CREATE INDEX IF NOT EXISTS idx_pred_seq_id ON predictions(sequence_id);
```

---

## Notes

- The database is read-only for the webapp. Never write to it from the API.
- For the webapp, expose only `sequences`, `predictions`, `confidence_scores`, `km_evidence`, `ec_evidence`. The `features_*` tables are large (14 GB total) and should not be queried in the webapp.
- `features_ankh` insertion is in progress as of April 2026 — do not depend on it being complete.
- All Km values are in **mM** and stored as both raw and log₁₀.
- `co2_prob` is the raw XGBoost output [0–1]. `confidence_label` is derived from it via thresholds: high ≥ 0.90, medium ≥ 0.70, low ≥ 0.50.
