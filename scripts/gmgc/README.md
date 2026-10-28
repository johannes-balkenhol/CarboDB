# GMGC10 marine scan pipeline

Scans the GMGC10 (Coelho et al. 2022) marine sub-catalog for carboxylase
candidates using the CarboDB v5 prediction cascade.

Source: https://gmgc.embl.de/
Paper: 10.1038/s41586-021-04233-4

## Scope
Marine sub-catalog only (130 samples). Future expansions may add
freshwater (104 samples) and soil (312 samples) as comparison habitats.

## Metadata
Available: GPS coordinates, habitat, country, sequencing depth metrics.
NOT available in GMGC10: water depth, temperature, salinity, nutrients.
(For depth/temperature analyses, would need Tara OM-RGC.v2 instead.)

## Pipeline order
01_download_marine_faa.sh        - Download marine .faa.gz from GMGC
02_extract_marine_unigenes.py    - Filter unigene-environment table to marine
03_filter_fasta.py               - Keep only marine sequences
04_chunk_fasta.py                - Split into chunks for parallel processing
05_pfam_prefilter.sh             - HMMER scan vs carboxylase Pfams
06_extract_candidates.py         - FASTA of Pfam-positive sequences only
07_predict_cascade.py            - CarboDB v5: binary → EC → Km
08_join_metadata.py              - Merge predictions with sample meta + taxonomy
09_analysis_geo.py               - GPS / latitude patterns
10_analysis_taxonomy.py          - Genus/family/order vs Km
11_analysis_habitat.py           - Cross-habitat (after Phase 2 expansion)

## Carboxylase Pfam IDs used as prefilter
PF00016  RuBisCO_large
PF02788  RuBisCO_large_N
PF00210  Ferritin (control, NOT carboxylase — sanity check)
PF00682  HMG_CoA_synt
PF00289  Biotin_carb_N
PF02786  Biotin_carb_C
PF02785  Biotin_carb_M
PF01039  Carboxyl_trans
[full list TBD - check feature_extraction_v5.py for what v5 actually uses]
