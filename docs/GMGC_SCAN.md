# GMGC10 metagenome scan

Replaces the planned direct Tara Oceans scan with GMGC10 (Coelho et al.
2022), as suggested by Simone.

GMGC10 includes Tara Oceans samples alongside other marine studies,
clustered at 95% nucleotide / 90% amino acid identity. Less redundant
than direct Tara scanning, with consistent pipeline across all 14
habitats.

## Status
- [x] Phase 0: metadata exploration complete (2026-05-21)
- [ ] Phase 1: marine sub-catalog scan (in progress)
- [ ] Phase 2: freshwater + soil expansion (planned)
- [ ] Phase 3: cross-habitat comparative analyses (planned)

## Source
- Catalog: https://gmgc.embl.de/
- Paper:   10.1038/s41586-021-04233-4

## Available sample metadata (limited)
- GPS coordinates (latitude, longitude)
- Habitat (14 categories)
- Country, income group
- Sequencing metrics (GC fraction, basepairs, ORFs)

## NOT available in GMGC10
- Water depth (surface / DCM / mesopelagic)
- Temperature, salinity, oxygen
- Nutrients (nitrate, phosphate, silicate)
- These exist in Tara native metadata (PANGAEA) but were dropped
  when GMGC was assembled. Out of scope for this scan.

## Pipeline
See scripts/gmgc/ for the numbered pipeline scripts.
