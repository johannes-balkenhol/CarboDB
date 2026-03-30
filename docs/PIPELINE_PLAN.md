#!/usr/bin/env python3
"""
CarboxyDB Clean Pipeline — Master Plan
=======================================
This file is the definitive reference. Read before running anything.

PIPELINE OVERVIEW
─────────────────
  Script 01 — BRENDA download (positive EC list + Km gold standard + BRENDA negatives)
  Script 02 — UniProt download (SwissProt + TrEMBL positives + large negative pool)
  Script 03 — Merge & deduplicate all sources → master TSV + FASTA
  Script 04 — Feature annotation (multi-layer, batched, modular)
  Script 05 — SQLite database build (evidence-tiered schema)
  [Script 06 — ML training — separate, runs after features are ready]

DATA SCALE TARGETS (confirmed from working runs + new goal)
────────────────────────────────────────────────────────────

  SOURCE           POSITIVES        NEGATIVES        NOTES
  BRENDA           695,384          3,324,901        From working Dec 2025 run
  SwissProt        ~50,000          ~50,000          reviewed:true AND ec:{co2_ec}
  TrEMBL           ~600,000         ~600,000         reviewed:false AND ec:{co2_ec}
  ─────────────────────────────────────────────────────────────────────────────
  TOTAL            ~1,345,000       ~3,974,901       Deduplicated by UniProt ID
  ML TRAINING      ~100,000 pos     ~500,000 neg     Balanced subset, CD-HIT 90%

  Gold standard Km (BRENDA):  ~3,001 unique sequences
  Km range:                   0.0008 – 83.0 mM

CONFIRMED CO2 EC CLASSES (found by scanning ALL BRENDA Km, Dec 2025)
─────────────────────────────────────────────────────────────────────
  Top by count:
    4.2.1.1   — Carbonic anhydrase (999 Km entries)
    4.1.1.39  — RuBisCO            (930)
    6.3.5.5   — Carbamoyl-P synth  (281)
    6.3.3.3   — PhosphoribosylGAR  (164)  ← unexpected!
    4.1.1.49  — PEPCK-GTP          (123)
    4.1.1.31  — PEPC               (122)
    6.4.1.2   — Acetyl-CoA carbox  (103)
    6.3.4.14  — Biotin carboxylase  (57)
    6.4.1.1   — Pyruvate carbox     (53)
    4.1.1.32  — PEPCK-ATP           (42)
    + 29 more EC classes (39 total)
  
  IMPORTANT: Do NOT hardcode this list. Script 01 discovers it dynamically
  by scanning all BRENDA Km entries. New BRENDA releases may add more.
  Use the ec_list saved by Script 01 as input to Script 02.

FEATURE LAYERS (confirmed working from CarboxylaseDatabase, Dec 2025)
───────────────────────────────────────────────────────────────────────
  Layer   Name              Features  Tool         Notes
  ─────────────────────────────────────────────────────────────────────
  A1      AA composition       20     BioPython    Standard
  A2      Dipeptide freq       400    BioPython    dp_XY format, confirmed top features
  A3      Pseudo-AAC (PseAAC)  30     Custom       pse_ prefix
  A4      Physicochemical      ~20    BioPython    MW, pI, GRAVY, aromaticity, instability
  A5      Catalytic core       ~17    Custom       inv_cat_* (middle 50% residue stats)
  A6      Interface features   ~12    Custom       interface_* (binding site proxies)
  ─────────────────────────────────────────────────────────────────────
  B1      Pfam domains         ~30    HMMER3       pfam_PFxxxxx binary + n_domains
  B2      PROSITE patterns      14    regex        prosite_PS* count + present
  ─────────────────────────────────────────────────────────────────────
  C       EC-specific motifs     7    regex        motif_rubisco_kk, motif_ca_hh, etc.
  ─────────────────────────────────────────────────────────────────────
  D       MEME de novo motifs   65    FIMO batch   binary hit/no-hit (PENDING subproject)
  ─────────────────────────────────────────────────────────────────────
  E       BLAST homology         4    blastp batch pident, evalue, best_ec, has_hit
  ─────────────────────────────────────────────────────────────────────
  F       ESM-2 embeddings    1280    HPC GPU      esm2_ prefix, mean-pooled
  ─────────────────────────────────────────────────────────────────────

  v3 feature set (no ESM-2): layers A+B+C+E = ~523 features confirmed
  v5 feature set (+ ESM-2):  layers A+B+C+E+F = ~1803 features

  MEME motifs (Layer D): injected once motif discovery pipeline is complete.
  Format: binary TSV with columns meme_{family}_{motif_id}, merged on uniprot_id.

CONFIRMED FEATURE PERFORMANCE (from Dec 2025 runs)
────────────────────────────────────────────────────
  Feature set          R² Km   RMSE   Top features
  474 features          0.89   0.46   dp_GF (42%), pse_aa_Y (5.5%), pse_aa_G (5%)
  545 features (+new)   0.91   0.43   dp_GF (42%), catalytic_K_percent (4.3%)
  + ESM-2               0.92+  0.40   esm2 dims dominate, then dipeptides

DATABASE SCHEMA (evidence-tiered, matches working v2 DB)
──────────────────────────────────────────────────────────
  Tables:
    sequences          — core (uniprot_id PK, sequence, length, organism, label)
    ec_evidence        — multiple sources per UID (source, evidence_tier, ec_number)
    km_evidence        — multiple sources per UID (km_value_mM, evidence_tier)
    features_composition  — AA comp + dipeptide JSON + PseAAC + physicochemical
    features_domains      — Pfam binary columns + JSON blob of all hits
    features_motifs       — PROSITE binary + MEME JSON blob
    features_blast        — pident, evalue, best_ec
    features_esm2         — 1280-dim float32 blob
    predictions           — v3/v5 binary + EC + Km predictions
    confidence_scores     — composite score per UID

  Evidence tiers:
    1 = experimental  (BRENDA Km measured)
    2 = curated       (SwissProt reviewed)
    3 = predicted     (TrEMBL or model output)
    4 = inferred      (BLAST/Pfam)

MOTIF SUBPROJECT (pending — Layer D)
─────────────────────────────────────
  When MEME motifs are ready, they are added as follows:
  1. MEME outputs per-family .meme files → saved to data/features/meme/
  2. Script: fimo_batch.py runs FIMO once on full FASTA per motif file
  3. Output: meme_hits.tsv with columns: uniprot_id, meme_{name}_hit (0/1)
  4. Script 04 detects data/features/meme/meme_hits.tsv and merges automatically
  5. Database: stored in features_motifs.meme_hits_json

  This file is the specification for the motif team:
    - Output format: TSV with header row
    - Column names: uniprot_id, meme_{family}_{motif_number}_{consensus_abbrev}
    - Values: 0 or 1 (FIMO p < 1e-4)
    - File location: data/features/meme/meme_hits.tsv
"""

# This is a documentation-only file. No executable code.
# See scripts 01–05 for the actual implementation.
print(__doc__)

## Km Gold Standard Note (March 30, 2026)
master.tsv contains 2,971 sequences with experimental Km:
- 2,077 from current BRENDA download (matched via genus+species join)
- 894 from CarboxyPred v3 legacy file (source='brenda_km_v3')
  These are valid full sequences (mean 385 aa) with experimental Km.
  Script 01 third pass will recover these automatically on future runs
  by fetching from UniProt REST API for any Km UID not in BRENDA sequences.
