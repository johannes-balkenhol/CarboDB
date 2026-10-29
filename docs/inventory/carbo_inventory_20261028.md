# Carbo inventory Wed Oct 28 09:21:33 PM CET 2026

## CarboDB_v3
```
size: 330G   owner: job37yv   modified: 2026-05-31
--- git
origin	git@github.com:johannes-balkenhol/CarboDB.git (fetch)
origin	git@github.com:johannes-balkenhol/CarboDB.git (push)
5b81111 2026-05-08 chore(motifs): commit v2 generator for traceability
f0440d0 2026-05-08 fix(browse): include uniprot_id field in /db/seq response
751d308 2026-05-08 feat(motifs): v3.1 — six per-Form balanced Km contrasts
c86bd24 2026-05-08 fix(motifs): correct Form classification rules
5a0f5fa 2026-05-07 docs: handoff bundle for shared-owner colleague + v3 motif bundle generator
A  docs/GMGC_SCAN.md
A  scripts/gmgc/README.md
 M webapp/app/pipeline/annotate.py
 M webapp/app/pipeline/predict.py
?? CarboDB-App-v2/
?? carbodb.sqlite
?? models/
?? scripts/gmgc/01_translate_marine.sh
?? scripts/gmgc/02_subsample_pilot.py
?? scripts/gmgc/03_run_hmmscan_parallel.sh
--- tree (depth 2)
CarboDB_v3
CarboDB_v3/00_setup_project.py
CarboDB_v3/CarboDB-App-v2
CarboDB_v3/CarboDB-App-v2/app
CarboDB_v3/CarboDB-App-v2/docker
CarboDB_v3/CarboDB-App-v2/.env
CarboDB_v3/CarboDB-App-v2/frontend
CarboDB_v3/CarboDB-App-v2/jobs
CarboDB_v3/CarboDB-App-v2/logs
CarboDB_v3/CarboDB-App-v2/models
CarboDB_v3/CarboDB-App-v2/README.md
CarboDB_v3/CarboDB-App-v2/README.md.OLD.before_rewrite
CarboDB_v3/CarboDB-App-v2/requirements.txt
CarboDB_v3/CarboDB-App-v2/resources
CarboDB_v3/CarboDB-App-v2/start_app.sh
CarboDB_v3/CarboDB-App-v2/static
CarboDB_v3/carbodb.sqlite
CarboDB_v3/config.py
CarboDB_v3/data
CarboDB_v3/database
CarboDB_v3/data/benchmark
CarboDB_v3/data/biological
CarboDB_v3/data/blast_ec_dbs
CarboDB_v3/data/blast_ec_dbs_exp
CarboDB_v3/data/dbs
CarboDB_v3/data/evolution
CarboDB_v3/data/external_annotations.sqlite
CarboDB_v3/data/features
CarboDB_v3/data/formI_clade_control.py
CarboDB_v3/data/gmgc
CarboDB_v3/data/greenase
CarboDB_v3/data/interim
CarboDB_v3/data/metagenome
CarboDB_v3/data/ml
CarboDB_v3/data/models
CarboDB_v3/data/motifs
CarboDB_v3/data/motifs_v2_clean
CarboDB_v3/data/motifs_v2_clean.tar.gz
CarboDB_v3/data/motifs_v3_1_per_form_balanced
CarboDB_v3/data/motifs_v3_1_per_form_balanced.tar.gz
CarboDB_v3/data/motifs_v3_form_split
CarboDB_v3/data/motifs_v3_form_split.tar.gz
CarboDB_v3/data/primary
CarboDB_v3/data/raw
CarboDB_v3/data/shap
CarboDB_v3/data/splits
CarboDB_v3/data/taxonomy
CarboDB_v3/docs
CarboDB_v3/docs/2026-03-31_carbodb_schema_v2.html
CarboDB_v3/docs/API.md
CarboDB_v3/docs/ARCHITECTURE.md
CarboDB_v3/docs/DATABASE.md
CarboDB_v3/docs/DATA_INGESTION.md
CarboDB_v3/docs/DEPLOYMENT.md
CarboDB_v3/docs/EXTERNAL_INTEGRATION.md
CarboDB_v3/docs/FRONTEND_SPEC.md
CarboDB_v3/docs/GMGC_SCAN.md
CarboDB_v3/docs/HANDOFF.md
CarboDB_v3/docs/MOTIF_ANALYSIS_v2.md
CarboDB_v3/docs/PIPELINE_PLAN.md
CarboDB_v3/docs/PLAN_2026-04-26.md
CarboDB_v3/docs/ROADMAP.md
CarboDB_v3/docs/WEBAPP.md
CarboDB_v3/environment.yml
CarboDB_v3/figures
CarboDB_v3/figures/benchmark_figure.pdf
CarboDB_v3/figures/benchmark_figure.png
CarboDB_v3/figures/biological_findings_figure.pdf
CarboDB_v3/figures/biological_findings_figure.png
CarboDB_v3/figures/carnivorous_plant_figure.png
CarboDB_v3/figures/ecological_findings_figure.pdf
CarboDB_v3/figures/ecological_findings_figure.png
CarboDB_v3/figures/evolution_figure.pdf
CarboDB_v3/figures/evolution_figure.png
CarboDB_v3/figures/feature_importance_figure.pdf
CarboDB_v3/figures/feature_importance_figure.png
CarboDB_v3/figures/tara_metagenome_figure.pdf
CarboDB_v3/figures/tara_metagenome_figure.png
CarboDB_v3/figures/taxonomy_figure.pdf
CarboDB_v3/figures/taxonomy_figure.png
--- scripts
CarboDB_v3/00_setup_project.py
CarboDB_v3/CarboDB-App-v2/app/__init__.py
CarboDB_v3/CarboDB-App-v2/app/main.py
CarboDB_v3/CarboDB-App-v2/app/startup.py
CarboDB_v3/CarboDB-App-v2/start_app.sh
CarboDB_v3/config.py
CarboDB_v3/data/formI_clade_control.py
CarboDB_v3/scripts/01_brenda_download.py
CarboDB_v3/scripts/02_uniprot_download.py
CarboDB_v3/scripts/03_merge_all_sources.py
CarboDB_v3/scripts/04a_composition.py
CarboDB_v3/scripts/04_annotate_features.py
CarboDB_v3/scripts/04b_hmmer.sh
CarboDB_v3/scripts/04c_interproscan.sh
CarboDB_v3/scripts/04d_ankh.sh
CarboDB_v3/scripts/04e_esm2.py
CarboDB_v3/scripts/04f_blast_benchmark.py
CarboDB_v3/scripts/05_build_database.py
CarboDB_v3/scripts/06_cluster_and_split.py
CarboDB_v3/scripts/07_build_feature_matrix.py
CarboDB_v3/scripts/08_train_models.py
CarboDB_v3/scripts/09_benchmark.py
CarboDB_v3/scripts/10_predict_all.py
CarboDB_v3/scripts/11_annotate_sequence.py
CarboDB_v3/scripts/12_shap_analysis.py
CarboDB_v3/scripts/13_benchmark_identity.py
CarboDB_v3/scripts/13_hard_tests.py
CarboDB_v3/scripts/13_publication_benchmark.py
CarboDB_v3/scripts/14_biological_analysis.py
CarboDB_v3/scripts/15_publication_figure.py
CarboDB_v3/scripts/16_figure_feature_importance.py
CarboDB_v3/scripts/17_figure_biological.py
CarboDB_v3/scripts/18_carnivorous_figure.py
CarboDB_v3/scripts/19_taxonomic_deepdive.py
CarboDB_v3/scripts/20_evolutionary_analysis.py
CarboDB_v3/scripts/21_ecological_findings_figure.py
CarboDB_v3/scripts/22_metagenome_scan.py
CarboDB_v3/scripts/23_tara_figure.py
CarboDB_v3/scripts/build_motif_bundles_v2.py
CarboDB_v3/scripts/build_motif_bundles_v3_1_per_form_balanced.py
CarboDB_v3/scripts/build_motif_bundles_v3_form_split.py
CarboDB_v3/scripts/explore_ecological.py
CarboDB_v3/scripts/explore.py
CarboDB_v3/scripts/gmgc/01_translate_marine.sh
CarboDB_v3/scripts/gmgc/02_subsample_pilot.py
CarboDB_v3/scripts/gmgc/03b_build_carboxy_hmm_db.sh
CarboDB_v3/scripts/gmgc/03b_run_hmmsearch_full.sh
CarboDB_v3/scripts/gmgc/03_run_hmmscan_parallel.sh
CarboDB_v3/scripts/gmgc/04b_filter_carboxylase_candidates_full.py
CarboDB_v3/scripts/gmgc/04_filter_carboxylase_candidates.py
CarboDB_v3/scripts/gmgc/05b_run_v5_cascade_full.py
CarboDB_v3/scripts/gmgc/05_run_v5_cascade.py
CarboDB_v3/scripts/gmgc/06b_join_metadata_full.py
CarboDB_v3/scripts/gmgc/06_join_metadata.py
CarboDB_v3/scripts/gmgc/07_pilot_figures.py
CarboDB_v3/scripts/gmgc/08_novelty_check.py
CarboDB_v3/scripts/gmgc/09b_per_ec_genus_full.py
CarboDB_v3/scripts/gmgc/09_per_ec_genus_stats.py
CarboDB_v3/scripts/gmgc/10_sample_ec_pca.py
CarboDB_v3/scripts/gmgc/11_km_random_forest.py
CarboDB_v3/scripts/gmgc/12_bitscore_rf.py
CarboDB_v3/scripts/gmgc/13b_rubisco_forms_rescue.py
CarboDB_v3/scripts/gmgc/13c_rubisco_form2_4_disambiguation.py
CarboDB_v3/scripts/gmgc/13_rubisco_forms.py
CarboDB_v3/scripts/gmgc/14_rubisco_geography.py
CarboDB_v3/scripts/gmgc/15b_ocean_basin_fixed.py
CarboDB_v3/scripts/gmgc/15_ocean_basin_analysis.py
CarboDB_v3/scripts/gmgc/16b_literature_benchmark_v2.py
CarboDB_v3/scripts/gmgc/16_literature_benchmark.py
CarboDB_v3/scripts/run_ankh_wbbi206.sh
CarboDB_v3/scripts/run_hmmer_wbbi206.sh
CarboDB_v3/scripts/run_interproscan_wbbi206.sh
CarboDB_v3/scripts/run_msfragger_swieringa.sh
CarboDB_v3/scripts/schema.sql
CarboDB_v3/webapp/app/__init__.py
CarboDB_v3/webapp/app/main.py
CarboDB_v3/webapp/app/startup.py
CarboDB_v3/webapp/scripts/build_ec_blast_dbs_experimental.py
CarboDB_v3/webapp/scripts/build_ec_blast_dbs.py
--- big files (>500 MB)
49.0 GB	CarboDB_v3/data/primary/carbodb.sqlite
17.3 GB	CarboDB_v3/data/gmgc/interim_full/marine.95nr.faa
13.6 GB	CarboDB_v3/data/gmgc/downloads/GMGC10.marine.95nr.fna.gz
9.8 GB	CarboDB_v3/data/ml/X_binary_train.npz
9.2 GB	CarboDB_v3/data/gmgc/interim/marine.95nr.faa.gz
5.9 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0-64-bit.tar.gz
5.2 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/gene3d/4.3.0/gene3d_main.hmm
3.4 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/panther/19.0/famhmm/binHmm
3.0 GB	CarboDB_v3/data/gmgc/downloads/GMGC10.meta.tsv.gz
2.5 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/gene3d/4.3.0/gene3d_main.hmm.h3p
2.3 GB	CarboDB_v3/data/gmgc/interim/marine_unigene_ids.txt
2.2 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/gene3d/4.3.0/gene3d_main.hmm.h3m
2.1 GB	CarboDB_v3/data/gmgc/downloads/GMGC10.taxonomy.tsv.gz
1.8 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/pfam/37.1/pfam_a.hmm
1.8 GB	CarboDB_v3/data/ml/X_ec_train_fixed.npz
1.7 GB	CarboDB_v3/data/ml/X_ec_train.npz
1.7 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/panther/19.0/famhmm/binHmm.h3p
1.5 GB	CarboDB_v3/data/primary/master.tsv
1.5 GB	CarboDB_v3/data/gmgc/downloads/GMGC10.unigene-environment.tsv.gz
1.4 GB	CarboDB_v3/data/dbs/interpro/interproscan-5.72-103.0/data/panther/19.0/famhmm/binHmm.h3m
--- databases / models
CarboDB_v3/data/primary/carbodb.sqlite
CarboDB_v3/data/benchmark/benchmark_report_v5.json
CarboDB_v3/data/models/binary_v5.json
CarboDB_v3/data/models/training_report_v5.json
CarboDB_v3/data/models/ec_v5.json
CarboDB_v3/data/models/km_v5.json
CarboDB_v3/data/models/km_v5_final.json
CarboDB_v3/data/models/km_v5_weighted.json
CarboDB_v3/data/external_annotations.sqlite
CarboDB_v3/webapp/models/binary_v5.json
CarboDB_v3/webapp/models/ec_v5.json
CarboDB_v3/webapp/models/km_v5_weighted.json
CarboDB_v3/CarboDB-App-v2/models/binary_v5.json
CarboDB_v3/CarboDB-App-v2/models/ec_v5.json
CarboDB_v3/CarboDB-App-v2/models/km_v5_weighted.json
CarboDB_v3/carbodb.sqlite
--- download provenance
CarboDB_v3/scripts/01_brenda_download.py
CarboDB_v3/scripts/02_uniprot_download.py
CarboDB_v3/scripts/04_annotate_features.py
CarboDB_v3/scripts/03_merge_all_sources.py
CarboDB_v3/scripts/05_build_database.py
CarboDB_v3/scripts/06_cluster_and_split.py
CarboDB_v3/scripts/build_motif_bundles_v2.py
CarboDB_v3/scripts/gmgc/01_translate_marine.sh
CarboDB_v3/scripts/gmgc/02_subsample_pilot.py
CarboDB_v3/scripts/gmgc/03_run_hmmscan_parallel.sh
CarboDB_v3/scripts/gmgc/04_filter_carboxylase_candidates.py
CarboDB_v3/scripts/gmgc/05_run_v5_cascade.py
CarboDB_v3/scripts/gmgc/06_join_metadata.py
CarboDB_v3/scripts/gmgc/07_pilot_figures.py
CarboDB_v3/scripts/gmgc/08_novelty_check.py
```

## CarboDB-App-v2
```
size: 192M   owner: job37yv   modified: 2026-06-24
--- git
origin	git@github.com:johannes-balkenhol/CarboDB-App.git (fetch)
origin	git@github.com:johannes-balkenhol/CarboDB-App.git (push)
a48a19f 2026-06-24 docs: runbook for operating the app (Aman handoff)
4343ef5 2026-05-24 chore: stop tracking __pycache__ bytecode
56e4844 2026-05-24 chore: gitignore Python __pycache__
ac32e8e 2026-05-24 feat(frontend): make vite ports configurable via env vars
6b1dc8a 2026-05-24 feat(frontend): make vite port and proxy target configurable via env
?? docs/.ipynb_checkpoints/
--- tree (depth 2)
CarboDB-App-v2
CarboDB-App-v2/app
CarboDB-App-v2/app/db
CarboDB-App-v2/app/__init__.py
CarboDB-App-v2/app/main.py
CarboDB-App-v2/app/pipeline
CarboDB-App-v2/app/routes
CarboDB-App-v2/app/startup.py
CarboDB-App-v2/backend
CarboDB-App-v2/backend/app
CarboDB-App-v2/backend/carboxylase_search
CarboDB-App-v2/backend/config.py
CarboDB-App-v2/backend/data_acquisition
CarboDB-App-v2/backend/Dockerfile
CarboDB-App-v2/backend/domain
CarboDB-App-v2/backend/__init__.py
CarboDB-App-v2/backend/main.py
CarboDB-App-v2/backend/repository
CarboDB-App-v2/backend/requirements.txt
CarboDB-App-v2/backend/tests
CarboDB-App-v2/docker
CarboDB-App-v2/docker-compose.yml
CarboDB-App-v2/docker/docker-compose.yml
CarboDB-App-v2/docker/Dockerfile
CarboDB-App-v2/docs
CarboDB-App-v2/docs/2026-06-24_AMAN_REPORT.md
CarboDB-App-v2/docs/.ipynb_checkpoints
CarboDB-App-v2/docs/SETUP_AMAN.md
CarboDB-App-v2/.env
CarboDB-App-v2/frontend
CarboDB-App-v2/frontend/Dockerfile
CarboDB-App-v2/frontend/index.html
CarboDB-App-v2/frontend/jsconfig.json
CarboDB-App-v2/frontend/package.json
CarboDB-App-v2/frontend/package-lock.json
CarboDB-App-v2/frontend/public
CarboDB-App-v2/frontend/README.md
CarboDB-App-v2/frontend/src
CarboDB-App-v2/frontend/utils
CarboDB-App-v2/frontend/vite.config.js
CarboDB-App-v2/frontend/vitest.config.js
CarboDB-App-v2/frontend/.vscode
CarboDB-App-v2/jobs
CarboDB-App-v2/logs
CarboDB-App-v2/logs/webapp.log
CarboDB-App-v2/models
CarboDB-App-v2/models/binary_v5.json
CarboDB-App-v2/models/ec_label_map.json
CarboDB-App-v2/models/ec_v5.json
CarboDB-App-v2/models/feature_names_binary.json
CarboDB-App-v2/models/feature_names_km.json
CarboDB-App-v2/models/km_v5_weighted.json
CarboDB-App-v2/README.md
CarboDB-App-v2/README.md.OLD.before_rewrite
CarboDB-App-v2/requirements.txt
CarboDB-App-v2/resources
CarboDB-App-v2/resources/carboxylases
CarboDB-App-v2/start_app.sh
CarboDB-App-v2/static
CarboDB-App-v2/static/css
CarboDB-App-v2/static/index.html
CarboDB-App-v2/static/js
--- scripts
CarboDB-App-v2/app/db/__init__.py
CarboDB-App-v2/app/db/models_old.py
CarboDB-App-v2/app/__init__.py
CarboDB-App-v2/app/main.py
CarboDB-App-v2/app/pipeline/annotate.py
CarboDB-App-v2/app/pipeline/carbodb_config.py
CarboDB-App-v2/app/pipeline/carbodb_features_ref.py
CarboDB-App-v2/app/pipeline/composition.py
CarboDB-App-v2/app/pipeline/feature_extraction_old.py
CarboDB-App-v2/app/pipeline/feature_extraction.py
CarboDB-App-v2/app/pipeline/feature_extraction_v5.py
CarboDB-App-v2/app/pipeline/__init__.py
CarboDB-App-v2/app/pipeline/predict.py
CarboDB-App-v2/app/routes/batch.py
CarboDB-App-v2/app/routes/browse.py
CarboDB-App-v2/app/routes/__init__.py
CarboDB-App-v2/app/routes/predict_old.py
CarboDB-App-v2/app/routes/predict.py
CarboDB-App-v2/app/startup.py
CarboDB-App-v2/backend/app/__init__.py
CarboDB-App-v2/backend/app/routes.py
CarboDB-App-v2/backend/app/tasks.py
CarboDB-App-v2/backend/app/utils.py
CarboDB-App-v2/backend/carboxylase_search/export_as_pdf.py
CarboDB-App-v2/backend/carboxylase_search/__init__.py
CarboDB-App-v2/backend/carboxylase_search/run_all_searches_task.py
CarboDB-App-v2/backend/carboxylase_search/search_utils.py
CarboDB-App-v2/backend/config.py
CarboDB-App-v2/backend/data_acquisition/data_acquisition_ebi.py
CarboDB-App-v2/backend/data_acquisition/__init__.py
CarboDB-App-v2/backend/domain/HmmerSearchResult.py
CarboDB-App-v2/backend/domain/HmmProfile.py
CarboDB-App-v2/backend/domain/PrositeSearchResult.py
CarboDB-App-v2/backend/__init__.py
CarboDB-App-v2/backend/main.py
CarboDB-App-v2/backend/repository/HmmProfileRepository.py
CarboDB-App-v2/backend/repository/PrositePatternRepository.py
CarboDB-App-v2/backend/tests/__init__.py
CarboDB-App-v2/start_app.sh
--- big files (>500 MB)
--- databases / models
CarboDB-App-v2/models/binary_v5.json
CarboDB-App-v2/models/ec_v5.json
CarboDB-App-v2/models/km_v5_weighted.json
--- download provenance
CarboDB-App-v2/backend/data_acquisition/data_acquisition_ebi.py
CarboDB-App-v2/app/pipeline/carbodb_config.py
```

## Carboxylase_Aparna
```
size: 330M   owner: lumpi   modified: 2026-09-23
--- tree (depth 2)
Carboxylase_Aparna
Carboxylase_Aparna/00_SETUP_AND_RUN.md
Carboxylase_Aparna/06_MEME_controlled_rerun
Carboxylase_Aparna/06_MEME_controlled_rerun/controls
Carboxylase_Aparna/06_MEME_controlled_rerun/dedup
Carboxylase_Aparna/06_MEME_controlled_rerun/features
Carboxylase_Aparna/06_MEME_controlled_rerun/gold
Carboxylase_Aparna/06_MEME_controlled_rerun/meme
Carboxylase_Aparna/06_MEME_controlled_rerun/replicates
Carboxylase_Aparna/06_MEME_controlled_rerun/scan
Carboxylase_Aparna/06_MEME_controlled_rerun/tmp
Carboxylase_Aparna/all_families_combined.fasta
Carboxylase_Aparna/aln_ACC.fasta
Carboxylase_Aparna/aln_PCC.fasta
Carboxylase_Aparna/aln_PEPC.fasta
Carboxylase_Aparna/aln_RuBisCO.fasta
Carboxylase_Aparna/aln_VKGC.fasta
Carboxylase_Aparna/annotate_log.txt
Carboxylase_Aparna/annotate_motifs_databases.py
Carboxylase_Aparna/answer_key.csv
Carboxylase_Aparna/carboxylase-motifs
Carboxylase_Aparna/carboxylase-motifs/carboxylase-motifs-repo
Carboxylase_Aparna/carboxylase-motifs.d0c1205c-2452-4f04-970d-dfadca2699f5
Carboxylase_Aparna/carboxylase-motifs/.ipynb_checkpoints
Carboxylase_Aparna/carboxylase-motifs/tests
Carboxylase_Aparna/carboxylase-motifs_v2
Carboxylase_Aparna/carboxylase-motifs_v2/carboxylase-motifs-repo
Carboxylase_Aparna/carboxylase-motifs_v2/config.example.yaml
Carboxylase_Aparna/carboxylase-motifs_v2/data
Carboxylase_Aparna/carboxylase-motifs_v2/demo.log
Carboxylase_Aparna/carboxylase-motifs_v2/docs
Carboxylase_Aparna/carboxylase-motifs_v2/.ipynb_checkpoints
Carboxylase_Aparna/carboxylase-motifs_v2/LICENSE
Carboxylase_Aparna/carboxylase-motifs_v2/README.md
Carboxylase_Aparna/carboxylase-motifs_v2/requirements.txt
Carboxylase_Aparna/carboxylase-motifs_v2/results
Carboxylase_Aparna/carboxylase-motifs_v2/run.sh
Carboxylase_Aparna/carboxylase-motifs_v2/scripts
Carboxylase_Aparna/carboxylase-motifs_v2/tests
Carboxylase_Aparna/carboxylase-motifs_v2/WORKFLOW.md
Carboxylase_Aparna/carboxylases_mixed.fasta
Carboxylase_Aparna/carboxylases_mixed.tsv
Carboxylase_Aparna/carboxylase_survey_v3.ipynb
Carboxylase_Aparna/carboxylase_survey_v4.ipynb
Carboxylase_Aparna/channel_analysis
Carboxylase_Aparna/channel_analysis/chimerax_scripts
Carboxylase_Aparna/channel_analysis/figures
Carboxylase_Aparna/channel_analysis/fpocket_runs
Carboxylase_Aparna/channel_analysis/motif_classification.csv
Carboxylase_Aparna/channel_analysis/pymol_scripts
Carboxylase_Aparna/channel_analysis/structures
Carboxylase_Aparna/fpocket_channel_analysis.py
Carboxylase_Aparna/fpocket_channel_analysis_v1.py
Carboxylase_Aparna/fpocket_log.txt
Carboxylase_Aparna/ground_truth.csv
Carboxylase_Aparna/.ipynb_checkpoints
Carboxylase_Aparna/.ipynb_checkpoints/00_SETUP_AND_RUN-checkpoint.md
Carboxylase_Aparna/.ipynb_checkpoints/aln_ACC-checkpoint.fasta
Carboxylase_Aparna/.ipynb_checkpoints/annotate_motifs_databases-checkpoint.py
Carboxylase_Aparna/.ipynb_checkpoints/carboxylase_survey_v3-checkpoint.ipynb
Carboxylase_Aparna/.ipynb_checkpoints/carboxylase_survey_v4-checkpoint.ipynb
Carboxylase_Aparna/.ipynb_checkpoints/meme_log-checkpoint.txt
Carboxylase_Aparna/.ipynb_checkpoints/meme_motif_discovery-checkpoint.py
Carboxylase_Aparna/.ipynb_checkpoints/tmp_ACC-checkpoint.fasta
Carboxylase_Aparna/Johannes_data
Carboxylase_Aparna/Johannes_data/high_km_top50_check.txt
Carboxylase_Aparna/Johannes_data/.ipynb_checkpoints
Carboxylase_Aparna/Johannes_data/km_stratified_meme_analysis.ipynb
Carboxylase_Aparna/Johannes_data/meme_km_results
Carboxylase_Aparna/Johannes_data/motifs
Carboxylase_Aparna/Johannes_data/motifs_v2_clean
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced
Carboxylase_Aparna/Johannes_data/Motifs_v3_1_per_form_balanced.ipynb
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced.tar.gz
Carboxylase_Aparna/make_chimerax_scripts.py
Carboxylase_Aparna/make_paper_figures.py
Carboxylase_Aparna/meme_experiment.py
Carboxylase_Aparna/meme_highkm_log.txt
Carboxylase_Aparna/meme_log.txt
Carboxylase_Aparna/meme_lowkm_log.txt
--- scripts
Carboxylase_Aparna/annotate_motifs_databases.py
Carboxylase_Aparna/carboxylase-motifs/carboxylase-motifs-repo/run.sh
Carboxylase_Aparna/carboxylase-motifs/tests/test_pipeline.py
Carboxylase_Aparna/carboxylase-motifs_v2/carboxylase-motifs-repo/run.sh
Carboxylase_Aparna/carboxylase-motifs_v2/run.sh
Carboxylase_Aparna/carboxylase-motifs_v2/scripts/discover_motifs.py
Carboxylase_Aparna/carboxylase-motifs_v2/scripts/make_figures.py
Carboxylase_Aparna/carboxylase-motifs_v2/scripts/map_to_structure.py
Carboxylase_Aparna/carboxylase-motifs_v2/tests/test_pipeline.py
Carboxylase_Aparna/carboxylase_survey_v3.ipynb
Carboxylase_Aparna/carboxylase_survey_v4.ipynb
Carboxylase_Aparna/fpocket_channel_analysis.py
Carboxylase_Aparna/fpocket_channel_analysis_v1.py
Carboxylase_Aparna/.ipynb_checkpoints/annotate_motifs_databases-checkpoint.py
Carboxylase_Aparna/.ipynb_checkpoints/carboxylase_survey_v3-checkpoint.ipynb
Carboxylase_Aparna/.ipynb_checkpoints/carboxylase_survey_v4-checkpoint.ipynb
Carboxylase_Aparna/.ipynb_checkpoints/meme_motif_discovery-checkpoint.py
Carboxylase_Aparna/Johannes_data/.ipynb_checkpoints/km_stratified_meme_analysis-checkpoint.ipynb
Carboxylase_Aparna/Johannes_data/.ipynb_checkpoints/Motifs_v3_1_per_form_balanced-checkpoint.ipynb
Carboxylase_Aparna/Johannes_data/km_stratified_meme_analysis.ipynb
Carboxylase_Aparna/Johannes_data/motifs_v2_clean/carbodb_motif_analysis_v2.ipynb
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/fetch_lineage.py
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/formI_clade_control.py
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/formI_clade_control_v2.py
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/form_I_within_km_analysis.ipynb
Carboxylase_Aparna/Johannes_data/Motifs_v3_1_per_form_balanced.ipynb
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/johannes_workflow_v2.ipynb
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/residue149_trajectory.py
Carboxylase_Aparna/make_chimerax_scripts.py
Carboxylase_Aparna/make_paper_figures.py
Carboxylase_Aparna/meme_experiment.py
Carboxylase_Aparna/meme_motif_discovery.py
--- big files (>500 MB)
--- databases / models
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/rf_results.pkl
--- download provenance
Carboxylase_Aparna/.ipynb_checkpoints/annotate_motifs_databases-checkpoint.py
Carboxylase_Aparna/fpocket_channel_analysis.py
Carboxylase_Aparna/fpocket_channel_analysis_v1.py
Carboxylase_Aparna/annotate_motifs_databases.py
Carboxylase_Aparna/Johannes_data/motifs_v3_1_per_form_balanced/fetch_lineage.py
```

## CarboDB_motifs_collab
```
size: 4.0K   owner: job37yv   modified: 2026-04-14
--- tree (depth 2)
CarboDB_motifs_collab
CarboDB_motifs_collab/sequences_from_johannes
--- scripts
--- big files (>500 MB)
--- databases / models
--- download provenance
```

## CarboxyPred_v2
```
size: 52G   owner: job37yv   modified: 2026-01-13
--- tree (depth 2)
CarboxyPred_v2
CarboxyPred_v2/binary_classifier_v3.pkl
CarboxyPred_v2/binary_classifier_v5.pkl
CarboxyPred_v2/carboxylase_webapp_db_v2.zip
CarboxyPred_v2/carboxylase_webapp_db.zip
CarboxyPred_v2/CHANGELOG.md
CarboxyPred_v2/CITATION.cff
CarboxyPred_v2/configs
CarboxyPred_v2/configs/default.yaml
CarboxyPred_v2/data
CarboxyPred_v2/data/consensus_positives
CarboxyPred_v2/data/downloaded
CarboxyPred_v2/data/external
CarboxyPred_v2/data/negative_controls
CarboxyPred_v2/data/processed
CarboxyPred_v2/data/raw
CarboxyPred_v2/data/README.md
CarboxyPred_v2/docs
CarboxyPred_v2/docs/FAIR_RDM
CarboxyPred_v2/docs/features.md
CarboxyPred_v2/docs/figures
CarboxyPred_v2/docs/index.md
CarboxyPred_v2/docs/installation.md
CarboxyPred_v2/docs/PIPELINE_FLOWCHART.md
CarboxyPred_v2/docs/usage.md
CarboxyPred_v2/ec_classifier_v3.pkl
CarboxyPred_v2/environment.yml
CarboxyPred_v2/FP_FN_SUMMARY_TABLE_20260113.md
CarboxyPred_v2/km_predictor_v3.pkl
CarboxyPred_v2/LICENSE
CarboxyPred_v2/logs
CarboxyPred_v2/logs/full_695k_20260107_190146.log
CarboxyPred_v2/Makefile
CarboxyPred_v2/migrate_for_webapp.py
CarboxyPred_v2/MIGRATION_REPORT.md
CarboxyPred_v2/mkdocs.yml
CarboxyPred_v2/models
CarboxyPred_v2/models/binary_classifier_v3.pkl
CarboxyPred_v2/models/binary_classifier_v5.pkl
CarboxyPred_v2/models_compat
CarboxyPred_v2/models_compat/binary_classifier_v3.pkl
CarboxyPred_v2/models_compat/binary_classifier_v5.pkl
CarboxyPred_v2/models_compat/ec_classifier_v3.pkl
CarboxyPred_v2/models_compat/km_predictor_v3.pkl
CarboxyPred_v2/models_compat.zip
CarboxyPred_v2/models/ec_classifier_v3.pkl
CarboxyPred_v2/models/ec_multiclass_classifier.pkl
CarboxyPred_v2/models/km_predictor_v3.pkl
CarboxyPred_v2/models/km_regressor.pkl
CarboxyPred_v2/models/README.md
CarboxyPred_v2/models/taxonomy_classifier.pkl
CarboxyPred_v2/models/v3_binary_classifier.pkl
CarboxyPred_v2/models/v5_binary_classifier.pkl
CarboxyPred_v2/notebooks
CarboxyPred_v2/PROJECT_FINAL_SUMMARY.md
CarboxyPred_v2/PROJECT_SUMMARY_20260107.md
CarboxyPred_v2/PROJECT_SUMMARY_20260108.md
CarboxyPred_v2/pyproject.toml
CarboxyPred_v2/README.md
CarboxyPred_v2/requirements.txt
CarboxyPred_v2/results
CarboxyPred_v2/results/classifications
CarboxyPred_v2/results/database_evaluation
CarboxyPred_v2/results/figures
CarboxyPred_v2/results/full_695k_scan
CarboxyPred_v2/results/km_predictions
CarboxyPred_v2/results/legacy_results
CarboxyPred_v2/results/meeting_reports
CarboxyPred_v2/results/phylogenetics
CarboxyPred_v2/results/reclassification_v2
CarboxyPred_v2/results/tables
CarboxyPred_v2/results/taxonomy_classifier
CarboxyPred_v2/results/webapp_export
CarboxyPred_v2/scripts
CarboxyPred_v2/scripts/01_full_discovery_pipeline.py
CarboxyPred_v2/scripts/01_multicore_corrected.py
CarboxyPred_v2/scripts/02_evolutionary_visualization.py
CarboxyPred_v2/scripts/04_classify_sequences.py
CarboxyPred_v2/scripts/2026-02-27_carboxylase_task.md
CarboxyPred_v2/scripts/75_full_uniprot_scan.py
--- scripts
CarboxyPred_v2/migrate_for_webapp.py
CarboxyPred_v2/scripts/01_full_discovery_pipeline.py
CarboxyPred_v2/scripts/01_multicore_corrected.py
CarboxyPred_v2/scripts/02_evolutionary_visualization.py
CarboxyPred_v2/scripts/04_classify_sequences.py
CarboxyPred_v2/scripts/75_full_uniprot_scan.py
CarboxyPred_v2/scripts/90_full_reclassification_pipeline.py
CarboxyPred_v2/scripts/91_reclassify_fixed.py
CarboxyPred_v2/scripts/92_full_reclassify.py
CarboxyPred_v2/scripts/93_process_695k_fasta.py
CarboxyPred_v2/scripts/94_post_process_results.py
CarboxyPred_v2/scripts/analysis/11_deep_analysis.py
CarboxyPred_v2/scripts/analysis/76_add_km_predictions.py
CarboxyPred_v2/scripts/analysis/77_km_with_ec.py
CarboxyPred_v2/scripts/analyze_km_by_taxonomy.py
CarboxyPred_v2/scripts/benchmarks/52_comprehensive_benchmark.py
CarboxyPred_v2/scripts/benchmarks/62_comprehensive_v5_benchmark.py
CarboxyPred_v2/scripts/benchmarks/66_final_benchmark.py
CarboxyPred_v2/scripts/build_phylogenetic_tree.py
CarboxyPred_v2/scripts/create_summary_docx.py
CarboxyPred_v2/scripts/evaluate_webapp_db_final.py
CarboxyPred_v2/scripts/evaluate_webapp_db_full.py
CarboxyPred_v2/scripts/evaluate_webapp_db.py
CarboxyPred_v2/scripts/extract_top_km_lists.py
CarboxyPred_v2/scripts/legacy/01_full_discovery_pipeline.py
CarboxyPred_v2/scripts/legacy/01_multicore_corrected.py
CarboxyPred_v2/scripts/legacy/01_multicore_discovery_fixed.py
CarboxyPred_v2/scripts/legacy/01_multicore_discovery.py
CarboxyPred_v2/scripts/legacy/02_evolutionary_visualization.py
CarboxyPred_v2/scripts/legacy/03_multi_level_validation.py
CarboxyPred_v2/scripts/legacy/04_quality_control.py
CarboxyPred_v2/scripts/legacy/05_improvement_pipeline.py
CarboxyPred_v2/scripts/legacy/06_discovery_v3.py
CarboxyPred_v2/scripts/legacy/07_qc_dashboard.py
CarboxyPred_v2/scripts/legacy/08_proper_evaluation.py
CarboxyPred_v2/scripts/legacy/09_large_scale_evaluation.py
CarboxyPred_v2/scripts/legacy/10_local_large_evaluation.py
CarboxyPred_v2/scripts/legacy/11_deep_analysis.py
CarboxyPred_v2/scripts/legacy/12_comparison_novel.py
CarboxyPred_v2/scripts/legacy/15_qc_benchmark.py
CarboxyPred_v2/scripts/legacy/16_full_benchmark.py
CarboxyPred_v2/scripts/legacy/18_overnight_benchmark.py
CarboxyPred_v2/scripts/legacy/19_parse_fixed.py
CarboxyPred_v2/scripts/legacy/20_investigation.py
CarboxyPred_v2/scripts/legacy/21_deep_fn_fp_analysis.py
CarboxyPred_v2/scripts/legacy/22_fair_benchmark.py
CarboxyPred_v2/scripts/legacy/24_fair_same_dataset.py
CarboxyPred_v2/scripts/legacy/25_comprehensive_benchmark.py
CarboxyPred_v2/scripts/legacy/26_standardized_benchmark.py
CarboxyPred_v2/scripts/legacy/27_simple_fair_comparison.py
CarboxyPred_v2/scripts/legacy/29_clean_benchmark_pipeline.py
CarboxyPred_v2/scripts/legacy/30_proper_methodology.py
CarboxyPred_v2/scripts/legacy/31_final_benchmark.py
CarboxyPred_v2/scripts/legacy/32_complete_benchmark.py
CarboxyPred_v2/scripts/legacy/33_stratified_benchmark.py
CarboxyPred_v2/scripts/legacy/34_stratified_benchmark_fixed.py
CarboxyPred_v2/scripts/legacy/36_circularity_analysis.py
CarboxyPred_v2/scripts/legacy/40_comprehensive_pipeline_plan.py
CarboxyPred_v2/scripts/legacy/44_optimize_full_pipeline.py
CarboxyPred_v2/scripts/legacy/45_benchmark_similarity_stratified.py
CarboxyPred_v2/scripts/legacy/46b_stage2_random.py
CarboxyPred_v2/scripts/legacy/46_full_uniprot_analysis.py
CarboxyPred_v2/scripts/legacy/47_full_uniprot_overnight.py
CarboxyPred_v2/scripts/legacy/48_complete_analysis.py
CarboxyPred_v2/scripts/legacy/49_benchmark_vs_traditional.py
CarboxyPred_v2/scripts/legacy/50_create_co2_database.py
CarboxyPred_v2/scripts/legacy/51_continue_prediction.py
CarboxyPred_v2/scripts/legacy/52_comprehensive_benchmark.py
CarboxyPred_v2/scripts/legacy/53_uniprot20k_benchmark.py
CarboxyPred_v2/scripts/legacy/54_multimodel_benchmark.py
CarboxyPred_v2/scripts/legacy/55_multimodel_benchmark_fixed.py
CarboxyPred_v2/scripts/legacy/56_definitive_benchmark.py
CarboxyPred_v2/scripts/legacy/57_analyze_training_data.py
CarboxyPred_v2/scripts/legacy/59_prepare_expanded_training.py
CarboxyPred_v2/scripts/legacy/60_expand_negatives.py
CarboxyPred_v2/scripts/legacy/61_train_v5_model.py
CarboxyPred_v2/scripts/legacy/62_comprehensive_v5_benchmark.py
CarboxyPred_v2/scripts/legacy/64_benchmark_existing_testsets.py
CarboxyPred_v2/scripts/legacy/65_benchmark_with_metadata.py
CarboxyPred_v2/scripts/legacy/66_final_benchmark.py
--- big files (>500 MB)
38.4 GB	CarboxyPred_v2/logs/full_695k_20260107_190146.log
8.8 GB	CarboxyPred_v2/results/full_695k_scan/full_predictions_20260107_190147.csv
1.9 GB	CarboxyPred_v2/data/external/pfam/Pfam-A.hmm
0.9 GB	CarboxyPred_v2/data/external/pfam/Pfam-A.hmm.h3p
0.8 GB	CarboxyPred_v2/data/external/pfam/Pfam-A.hmm.h3m
0.6 GB	CarboxyPred_v2/results/webapp_export/carboxypred_v2_20260109_100916.db
0.5 GB	CarboxyPred_v2/webapp_database/carboxylase_webapp_v3.db
--- databases / models
CarboxyPred_v2/models/binary_classifier_v3.pkl
CarboxyPred_v2/models/ec_classifier_v3.pkl
CarboxyPred_v2/models/km_predictor_v3.pkl
CarboxyPred_v2/models/binary_classifier_v5.pkl
CarboxyPred_v2/models/v3_binary_classifier.pkl
CarboxyPred_v2/models/v5_binary_classifier.pkl
CarboxyPred_v2/models/ec_multiclass_classifier.pkl
CarboxyPred_v2/models/km_regressor.pkl
CarboxyPred_v2/models/taxonomy_classifier.pkl
CarboxyPred_v2/results/legacy_results/complete_benchmark/domain_ec_mappings.pkl
CarboxyPred_v2/results/legacy_results/final_benchmark/pfam_ec_mappings.pkl
CarboxyPred_v2/results/webapp_export/carboxypred_v2_20260109_100916.db
CarboxyPred_v2/models_compat/binary_classifier_v3.pkl
CarboxyPred_v2/models_compat/ec_classifier_v3.pkl
CarboxyPred_v2/models_compat/km_predictor_v3.pkl
CarboxyPred_v2/models_compat/binary_classifier_v5.pkl
CarboxyPred_v2/webapp_database/carboxylase.db
CarboxyPred_v2/webapp_database/carboxylase_v2.db
CarboxyPred_v2/webapp_database/carboxylase_webapp_v3.db
CarboxyPred_v2/ec_classifier_v3.pkl
--- download provenance
CarboxyPred_v2/setup/download_databases.sh
CarboxyPred_v2/scripts/legacy/01_full_discovery_pipeline.py
CarboxyPred_v2/scripts/legacy/05_improvement_pipeline.py
CarboxyPred_v2/scripts/legacy/07_qc_dashboard.py
CarboxyPred_v2/scripts/legacy/08_proper_evaluation.py
CarboxyPred_v2/scripts/legacy/09_large_scale_evaluation.py
CarboxyPred_v2/scripts/legacy/11_deep_analysis.py
CarboxyPred_v2/scripts/legacy/12_comparison_novel.py
CarboxyPred_v2/scripts/legacy/15_qc_benchmark.py
CarboxyPred_v2/scripts/legacy/16_full_benchmark.py
CarboxyPred_v2/scripts/legacy/18_overnight_benchmark.py
CarboxyPred_v2/scripts/legacy/25_comprehensive_benchmark.py
CarboxyPred_v2/scripts/legacy/26_standardized_benchmark.py
CarboxyPred_v2/scripts/legacy/29_clean_benchmark_pipeline.py
CarboxyPred_v2/scripts/legacy/30_proper_methodology.py
```

## CarboxyPred_Final
```
size: 27G   owner: job37yv   modified: 2025-12-22
--- tree (depth 2)
CarboxyPred_Final
CarboxyPred_Final/analysis
CarboxyPred_Final/analysis/improvements
CarboxyPred_Final/analyze_channel3_distribution.py
CarboxyPred_Final/BENCHMARK_SUMMARY.md
CarboxyPred_Final/data
CarboxyPred_Final/databases
CarboxyPred_Final/databases/Pfam-A.hmm
CarboxyPred_Final/databases/Pfam-A.hmm.h3f
CarboxyPred_Final/databases/Pfam-A.hmm.h3i
CarboxyPred_Final/databases/Pfam-A.hmm.h3m
CarboxyPred_Final/databases/Pfam-A.hmm.h3p
CarboxyPred_Final/databases/swissprot_db.pdb
CarboxyPred_Final/databases/swissprot_db.phr
CarboxyPred_Final/databases/swissprot_db.pin
CarboxyPred_Final/databases/swissprot_db.pog
CarboxyPred_Final/databases/swissprot_db.pos
CarboxyPred_Final/databases/swissprot_db.pot
CarboxyPred_Final/databases/swissprot_db.psq
CarboxyPred_Final/databases/swissprot_db.ptf
CarboxyPred_Final/databases/swissprot_db.pto
CarboxyPred_Final/databases/uniprot_sprot.fasta
CarboxyPred_Final/data/domains
CarboxyPred_Final/data/processed
CarboxyPred_Final/data/raw
CarboxyPred_Final/data/v5_training
CarboxyPred_Final/docs
CarboxyPred_Final/docs/BENCHMARK_FINAL.md
CarboxyPred_Final/docs/BENCHMARK_REPRODUCIBILITY.md
CarboxyPred_Final/docs/BLAST_vs_PFAM_ANALYSIS.md
CarboxyPred_Final/docs/CORRECTED_RESULTS.md
CarboxyPred_Final/docs/FAQ.md
CarboxyPred_Final/docs/FINAL_CORRECTED_SUMMARY.md
CarboxyPred_Final/docs/FINAL_RESULTS.md
CarboxyPred_Final/docs/KM_BIOLOGICAL_VALIDATION.md
CarboxyPred_Final/docs/KM_VALIDATION_REPORT.md
CarboxyPred_Final/docs/PAPER_RESULTS.md
CarboxyPred_Final/docs/PIPELINE.md
CarboxyPred_Final/docs/QUICK_STATUS.md
CarboxyPred_Final/docs/RESEARCH_PLAN.md
CarboxyPred_Final/docs/SUMMARY.md
CarboxyPred_Final/docs/VERSION_HISTORY.md
CarboxyPred_Final/figures
CarboxyPred_Final/final_package
CarboxyPred_Final/final_package/data
CarboxyPred_Final/final_package/FAQ.md
CarboxyPred_Final/final_package/FINAL_RESULTS.md
CarboxyPred_Final/final_package/models
CarboxyPred_Final/final_package/PIPELINE.md
CarboxyPred_Final/final_package/README.md
CarboxyPred_Final/final_package/scripts
CarboxyPred_Final/final_package/SUMMARY.md
CarboxyPred_Final/logs
CarboxyPred_Final/logs/add_km_predictions.log
CarboxyPred_Final/logs/benchmark_comparison.log
CarboxyPred_Final/logs/benchmark_existing.log
CarboxyPred_Final/logs/benchmark_similarity.log
CarboxyPred_Final/logs/benchmark_with_labels.log
CarboxyPred_Final/logs/clean_benchmark.log
CarboxyPred_Final/logs_clean_benchmark.txt
CarboxyPred_Final/logs/clean_benchmark_v2.log
CarboxyPred_Final/logs_comparison.txt
CarboxyPred_Final/logs_complete_benchmark.txt
CarboxyPred_Final/logs/comprehensive_benchmark.log
CarboxyPred_Final/logs_comprehensive.txt
CarboxyPred_Final/logs/continue_prediction.log
CarboxyPred_Final/logs_corrected.txt
CarboxyPred_Final/logs_deep_analysis.txt
CarboxyPred_Final/logs/definitive_benchmark.log
CarboxyPred_Final/logs_discovery.txt
CarboxyPred_Final/logs/discovery_v3_20251223_022040.log
CarboxyPred_Final/logs_evaluation.txt
CarboxyPred_Final/logs/expand_negatives.log
CarboxyPred_Final/logs_fair_benchmark.txt
CarboxyPred_Final/logs_fair_same.txt
CarboxyPred_Final/logs/final_benchmark.log
CarboxyPred_Final/logs_final_benchmark.txt
CarboxyPred_Final/logs_fixed.txt
CarboxyPred_Final/logs/full_analysis.log
CarboxyPred_Final/logs_full_benchmark.txt
--- scripts
CarboxyPred_Final/analyze_channel3_distribution.py
CarboxyPred_Final/final_package/scripts/benchmark.py
CarboxyPred_Final/final_package/scripts/download_data.py
CarboxyPred_Final/final_package/scripts/predict.py
CarboxyPred_Final/final_package/scripts/train.py
CarboxyPred_Final/scripts/01_full_discovery_pipeline.py
CarboxyPred_Final/scripts/01_multicore_corrected.py
CarboxyPred_Final/scripts/01_multicore_discovery_fixed.py
CarboxyPred_Final/scripts/01_multicore_discovery.py
CarboxyPred_Final/scripts/02_evolutionary_visualization.py
CarboxyPred_Final/scripts/03_multi_level_validation.py
CarboxyPred_Final/scripts/04_quality_control.py
CarboxyPred_Final/scripts/05_improvement_pipeline.py
CarboxyPred_Final/scripts/06_discovery_v3.py
CarboxyPred_Final/scripts/07_qc_dashboard.py
CarboxyPred_Final/scripts/08_proper_evaluation.py
CarboxyPred_Final/scripts/09_large_scale_evaluation.py
CarboxyPred_Final/scripts/10_local_large_evaluation.py
CarboxyPred_Final/scripts/11_deep_analysis.py
CarboxyPred_Final/scripts/12_comparison_novel.py
CarboxyPred_Final/scripts/15_qc_benchmark.py
CarboxyPred_Final/scripts/16_full_benchmark.py
CarboxyPred_Final/scripts/18_overnight_benchmark.py
CarboxyPred_Final/scripts/19_parse_fixed.py
CarboxyPred_Final/scripts/20_investigation.py
CarboxyPred_Final/scripts/21_deep_fn_fp_analysis.py
CarboxyPred_Final/scripts/22_fair_benchmark.py
CarboxyPred_Final/scripts/24_fair_same_dataset.py
CarboxyPred_Final/scripts/25_comprehensive_benchmark.py
CarboxyPred_Final/scripts/26_standardized_benchmark.py
CarboxyPred_Final/scripts/27_simple_fair_comparison.py
CarboxyPred_Final/scripts/29_clean_benchmark_pipeline.py
CarboxyPred_Final/scripts/30_proper_methodology.py
CarboxyPred_Final/scripts/31_final_benchmark.py
CarboxyPred_Final/scripts/32_complete_benchmark.py
CarboxyPred_Final/scripts/33_stratified_benchmark.py
CarboxyPred_Final/scripts/34_stratified_benchmark_fixed.py
CarboxyPred_Final/scripts/36_circularity_analysis.py
CarboxyPred_Final/scripts/40_comprehensive_pipeline_plan.py
CarboxyPred_Final/scripts/44_optimize_full_pipeline.py
CarboxyPred_Final/scripts/45_benchmark_similarity_stratified.py
CarboxyPred_Final/scripts/46b_stage2_random.py
CarboxyPred_Final/scripts/46_full_uniprot_analysis.py
CarboxyPred_Final/scripts/47_full_uniprot_overnight.py
CarboxyPred_Final/scripts/48_complete_analysis.py
CarboxyPred_Final/scripts/49_benchmark_vs_traditional.py
CarboxyPred_Final/scripts/50_create_co2_database.py
CarboxyPred_Final/scripts/51_continue_prediction.py
CarboxyPred_Final/scripts/52_comprehensive_benchmark.py
CarboxyPred_Final/scripts/53_uniprot20k_benchmark.py
CarboxyPred_Final/scripts/54_multimodel_benchmark.py
CarboxyPred_Final/scripts/55_multimodel_benchmark_fixed.py
CarboxyPred_Final/scripts/56_definitive_benchmark.py
CarboxyPred_Final/scripts/57_analyze_training_data.py
CarboxyPred_Final/scripts/59_prepare_expanded_training.py
CarboxyPred_Final/scripts/60_expand_negatives.py
CarboxyPred_Final/scripts/61_train_v5_model.py
CarboxyPred_Final/scripts/62_comprehensive_v5_benchmark.py
CarboxyPred_Final/scripts/64_benchmark_existing_testsets.py
CarboxyPred_Final/scripts/65_benchmark_with_metadata.py
CarboxyPred_Final/scripts/66_final_benchmark.py
CarboxyPred_Final/scripts/67_analyze_fp.py
CarboxyPred_Final/scripts/68_reproducibility_test.py
CarboxyPred_Final/scripts/69_reproducibility_fixed.py
CarboxyPred_Final/scripts/70_investigate_fp.py
CarboxyPred_Final/scripts/71_investigate_unknown.py
CarboxyPred_Final/scripts/72_clean_benchmark.py
CarboxyPred_Final/scripts/73_clean_benchmark_v2.py
CarboxyPred_Final/scripts/74_full_swissprot_benchmark.py
CarboxyPred_Final/scripts/75_full_uniprot_scan.py
CarboxyPred_Final/scripts/76_add_km_predictions.py
CarboxyPred_Final/scripts/77_km_with_ec.py
CarboxyPred_Final/scripts/migrate_carboxypred.sh
CarboxyPred_Final/scripts/run_improvements.sh
CarboxyPred_Final/scripts/setup_carboxypred_project.sh
--- big files (>500 MB)
17.3 GB	CarboxyPred_Final/logs_discovery.txt
4.7 GB	CarboxyPred_Final/logs_local_eval.txt
1.9 GB	CarboxyPred_Final/databases/Pfam-A.hmm
0.9 GB	CarboxyPred_Final/databases/Pfam-A.hmm.h3p
0.8 GB	CarboxyPred_Final/databases/Pfam-A.hmm.h3m
--- databases / models
CarboxyPred_Final/results/final_benchmark/pfam_ec_mappings.pkl
CarboxyPred_Final/results/complete_benchmark/domain_ec_mappings.pkl
CarboxyPred_Final/models_v3_improved/binary_classifier_v3.pkl
CarboxyPred_Final/models_v3_improved/ec_classifier_v3.pkl
CarboxyPred_Final/models_v3_improved/km_predictor_v3.pkl
CarboxyPred_Final/models_v4_optimized/binary_classifier.pkl
CarboxyPred_Final/models_v4_optimized/km_predictor.pkl
CarboxyPred_Final/models_v4_optimized/scaler.pkl
CarboxyPred_Final/models_v4_optimized/feature_cols.pkl
CarboxyPred_Final/models_v4_optimized/top_domains.pkl
CarboxyPred_Final/models_v4_optimized/domain_features.pkl
CarboxyPred_Final/models_v4_fixed/binary_classifier.pkl
CarboxyPred_Final/models_v4_fixed/km_predictor.pkl
CarboxyPred_Final/models_v4_fixed/scaler.pkl
CarboxyPred_Final/models_v4_fixed/feature_cols.pkl
CarboxyPred_Final/models_v4_fixed/top_domains.pkl
CarboxyPred_Final/models_v4_fixed/domain_features.pkl
CarboxyPred_Final/models_archive/v3_improved/binary_classifier_v3.pkl
CarboxyPred_Final/models_archive/v3_improved/ec_classifier_v3.pkl
CarboxyPred_Final/models_archive/v3_improved/km_predictor_v3.pkl
--- download provenance
CarboxyPred_Final/scripts/07_qc_dashboard.py
CarboxyPred_Final/scripts/05_improvement_pipeline.py
CarboxyPred_Final/scripts/08_proper_evaluation.py
CarboxyPred_Final/scripts/09_large_scale_evaluation.py
CarboxyPred_Final/scripts/11_deep_analysis.py
CarboxyPred_Final/scripts/12_comparison_novel.py
CarboxyPred_Final/scripts/15_qc_benchmark.py
CarboxyPred_Final/scripts/16_full_benchmark.py
CarboxyPred_Final/scripts/18_overnight_benchmark.py
CarboxyPred_Final/scripts/25_comprehensive_benchmark.py
CarboxyPred_Final/scripts/26_standardized_benchmark.py
CarboxyPred_Final/scripts/29_clean_benchmark_pipeline.py
CarboxyPred_Final/scripts/30_proper_methodology.py
CarboxyPred_Final/scripts/31_final_benchmark.py
CarboxyPred_Final/scripts/32_complete_benchmark.py
```

## CarboxyPred_Final_baseline_v1
```
size: 955M   owner: job37yv   modified: 2025-12-23
--- tree (depth 2)
CarboxyPred_Final_baseline_v1
CarboxyPred_Final_baseline_v1/analysis
CarboxyPred_Final_baseline_v1/analysis/improvements
CarboxyPred_Final_baseline_v1/data
CarboxyPred_Final_baseline_v1/data/domains
CarboxyPred_Final_baseline_v1/data/processed
CarboxyPred_Final_baseline_v1/data/raw
CarboxyPred_Final_baseline_v1/figures
CarboxyPred_Final_baseline_v1/logs_corrected.txt
CarboxyPred_Final_baseline_v1/logs_discovery.txt
CarboxyPred_Final_baseline_v1/logs_fixed.txt
CarboxyPred_Final_baseline_v1/logs_multicore.txt
CarboxyPred_Final_baseline_v1/models
CarboxyPred_Final_baseline_v1/paper
CarboxyPred_Final_baseline_v1/README.md
CarboxyPred_Final_baseline_v1/results
CarboxyPred_Final_baseline_v1/results/full_discovery
CarboxyPred_Final_baseline_v1/results/full_discovery_corrected
CarboxyPred_Final_baseline_v1/results/validation
CarboxyPred_Final_baseline_v1/scripts
CarboxyPred_Final_baseline_v1/scripts/01_full_discovery_pipeline.py
CarboxyPred_Final_baseline_v1/scripts/01_multicore_corrected.py
CarboxyPred_Final_baseline_v1/scripts/01_multicore_discovery_fixed.py
CarboxyPred_Final_baseline_v1/scripts/01_multicore_discovery.py
CarboxyPred_Final_baseline_v1/scripts/02_evolutionary_visualization.py
CarboxyPred_Final_baseline_v1/scripts/03_multi_level_validation.py
CarboxyPred_Final_baseline_v1/scripts/04_quality_control.py
--- scripts
CarboxyPred_Final_baseline_v1/scripts/01_full_discovery_pipeline.py
CarboxyPred_Final_baseline_v1/scripts/01_multicore_corrected.py
CarboxyPred_Final_baseline_v1/scripts/01_multicore_discovery_fixed.py
CarboxyPred_Final_baseline_v1/scripts/01_multicore_discovery.py
CarboxyPred_Final_baseline_v1/scripts/02_evolutionary_visualization.py
CarboxyPred_Final_baseline_v1/scripts/03_multi_level_validation.py
CarboxyPred_Final_baseline_v1/scripts/04_quality_control.py
--- big files (>500 MB)
--- databases / models
--- download provenance
CarboxyPred_Final_baseline_v1/scripts/01_full_discovery_pipeline.py
```

## CarboxyPred_Clean
```
size: 28M   owner: job37yv   modified: 2026-01-05
--- tree (depth 2)
CarboxyPred_Clean
CarboxyPred_Clean/00_config.py
CarboxyPred_Clean/data
CarboxyPred_Clean/data/domains
CarboxyPred_Clean/data/processed
CarboxyPred_Clean/data/raw
CarboxyPred_Clean/figures
CarboxyPred_Clean/logs_evolution.txt
CarboxyPred_Clean/logs_validation.txt
CarboxyPred_Clean/models
CarboxyPred_Clean/paper
CarboxyPred_Clean/README.md
CarboxyPred_Clean/results
CarboxyPred_Clean/results/discovery
CarboxyPred_Clean/results/evolution
CarboxyPred_Clean/results/full_discovery
CarboxyPred_Clean/results/validation
CarboxyPred_Clean/run_pipeline.sh
CarboxyPred_Clean/scripts
CarboxyPred_Clean/scripts/08_multi_level_validation.py
CarboxyPred_Clean/scripts/10_evolutionary_tree.py
--- scripts
CarboxyPred_Clean/00_config.py
CarboxyPred_Clean/run_pipeline.sh
CarboxyPred_Clean/scripts/08_multi_level_validation.py
CarboxyPred_Clean/scripts/10_evolutionary_tree.py
--- big files (>500 MB)
--- databases / models
CarboxyPred_Clean/results/full_discovery/checkpoint.pkl
--- download provenance
CarboxyPred_Clean/00_config.py
```

## CarboxylaseDatabase
```
size: 9.5G   owner: job37yv   modified: 2025-12-22
--- tree (depth 2)
CarboxylaseDatabase
CarboxylaseDatabase/10,000
CarboxylaseDatabase/analysis
CarboxylaseDatabase/data
CarboxylaseDatabase/data/alphafold_structures
CarboxylaseDatabase/data/blast_db
CarboxylaseDatabase/data/domain_annotations_complete.json
CarboxylaseDatabase/data/esmfold_structures
CarboxylaseDatabase/data/negative_sequences
CarboxylaseDatabase/data/pfam_hmms
CarboxylaseDatabase/data/pfam_test
CarboxylaseDatabase/data/processed_data
CarboxylaseDatabase/data/raw_data
CarboxylaseDatabase/data/reference_rubisco_sequences.json
CarboxylaseDatabase/data/rubisco_controls
CarboxylaseDatabase/data/rubisco_controls_pdb
CarboxylaseDatabase/data/rubisco_controls_v2
CarboxylaseDatabase/data/rubisco_mutants
CarboxylaseDatabase/data/rubisco_mutants_pdb
CarboxylaseDatabase/data/rubisco_optimized
CarboxylaseDatabase/data/rubisco_pdb_fixed
CarboxylaseDatabase/data/rubisco_pdb_with_mg
CarboxylaseDatabase/data/uniprot_domain_cache.json
CarboxylaseDatabase/docs
CarboxylaseDatabase/FINAL_PROJECT_SUMMARY.md
CarboxylaseDatabase/logs
CarboxylaseDatabase/logs/discovery.log
CarboxylaseDatabase/logs/esmfold_api.log
CarboxylaseDatabase/logs/esmfold_gpu.log
CarboxylaseDatabase/logs/esmfold.log
CarboxylaseDatabase/logs/full_pipeline.log
CarboxylaseDatabase/logs/train_improved.log
CarboxylaseDatabase/models
CarboxylaseDatabase/models/binary_carboxylase_classifier.pkl
CarboxylaseDatabase/models/combined_domain_dipeptide_classifier.pkl
CarboxylaseDatabase/models/domain_classifier.pkl
CarboxylaseDatabase/models/ec_classifier.pkl
CarboxylaseDatabase/models/hierarchical_km_model.pkl
CarboxylaseDatabase/models/km_final_predictor.pkl
CarboxylaseDatabase/models/km_prediction_pipeline.pkl
CarboxylaseDatabase/models/km_thorough_best.pkl
CarboxylaseDatabase/models/km_with_domains.pkl
CarboxylaseDatabase/models/mlp_final_classifier.pkl
CarboxylaseDatabase/models_v2
CarboxylaseDatabase/models_v2/binary_classifier.pkl
CarboxylaseDatabase/models_v2/ec_classifier.pkl
CarboxylaseDatabase/models_v2/km_predictor.pkl
CarboxylaseDatabase/PROJECT_FINAL_SUMMARY.md
CarboxylaseDatabase/README.md
CarboxylaseDatabase/results
CarboxylaseDatabase/results_classification
CarboxylaseDatabase/results_combined
CarboxylaseDatabase/results_combined/3d_correlations.csv
CarboxylaseDatabase/results_combined/3d_feature_correlations.csv
CarboxylaseDatabase/results_combined/3d_within_ec_correlations.csv
CarboxylaseDatabase/results_combined/ANALYSIS_SUMMARY.md
CarboxylaseDatabase/results_combined/binary_classifier_summary.json
CarboxylaseDatabase/results_combined/carboxylase_domain_analysis.json
CarboxylaseDatabase/results_combined/catalytic_feature_analysis.csv
CarboxylaseDatabase/results_combined/classifier_optimization_results.json
CarboxylaseDatabase/results_combined/complete_benchmark_summary.json
CarboxylaseDatabase/results_combined/complete_feature_table.json
CarboxylaseDatabase/results_combined/complete_feature_table.md
CarboxylaseDatabase/results_combined/comprehensive_feature_analysis.png
CarboxylaseDatabase/results_combined/comprehensive_learning_strategy.png
CarboxylaseDatabase/results_combined/domain_annotation_results.json
CarboxylaseDatabase/results_combined/domain_invariant_all_features.csv
CarboxylaseDatabase/results_combined/domain_invariant_features.csv
CarboxylaseDatabase/results_combined/domain_ml_classifier_simulated.json
CarboxylaseDatabase/results_combined/ec_results_20251215_171328.json
CarboxylaseDatabase/results_combined/ec_specific_features.pdf
CarboxylaseDatabase/results_combined/ec_specific_features.png
CarboxylaseDatabase/results_combined/ec_specific_performance.csv
CarboxylaseDatabase/results_combined/esm2_dipeptides_comparison.csv
CarboxylaseDatabase/results_combined/esm2_vs_dipeptides_comparison.png
CarboxylaseDatabase/results_combined/feature_importance_20251215_171328.csv
CarboxylaseDatabase/results_combined/feature_importance_20251215_171328.png
CarboxylaseDatabase/results_combined/feature_importance_4_1_1_31.csv
CarboxylaseDatabase/results_combined/feature_importance_4_1_1_39.csv
CarboxylaseDatabase/results_combined/feature_importance_4_2_1_1.csv
--- scripts
CarboxylaseDatabase/data/rubisco_pdb_fixed/20_investigation.py
CarboxylaseDatabase/data/rubisco_pdb_fixed/21_deep_fn_fp_analysis.py
CarboxylaseDatabase/data/rubisco_pdb_fixed/benchmark_v4_fast.py
CarboxylaseDatabase/data/rubisco_pdb_fixed/lda_cifar_benchmark.py
CarboxylaseDatabase/run_discovery.sh
CarboxylaseDatabase/run_full_pipeline.sh
CarboxylaseDatabase/scripts/python/01_brenda_download.py
CarboxylaseDatabase/scripts/python/02_clean_and_cluster.py
CarboxylaseDatabase/scripts/python/03_create_km_dataset.py
CarboxylaseDatabase/scripts/python/04_create_curated_km_dataset.py
CarboxylaseDatabase/scripts/python/06_advanced_features.py
CarboxylaseDatabase/scripts/python/07_rubisco_deep_analysis.py
CarboxylaseDatabase/scripts/python/08_ml_framework.py
CarboxylaseDatabase/scripts/python/08_ml_framework_v2.py
CarboxylaseDatabase/scripts/python/100_final_domain_summary.py
CarboxylaseDatabase/scripts/python/101_full_benchmark_mlp.py
CarboxylaseDatabase/scripts/python/102_full_benchmark_with_km.py
CarboxylaseDatabase/scripts/python/103b_km_mlp_fixed.py
CarboxylaseDatabase/scripts/python/103_km_mlp_thorough.py
CarboxylaseDatabase/scripts/python/104_feature_analysis_and_improvement.py
CarboxylaseDatabase/scripts/python/105_final_documentation.py
CarboxylaseDatabase/scripts/python/106_prediction_pipeline.py
CarboxylaseDatabase/scripts/python/107_improved_pipeline.py
CarboxylaseDatabase/scripts/python/108_train_ec_classifier.py
CarboxylaseDatabase/scripts/python/109_discover_new_carboxylases.py
CarboxylaseDatabase/scripts/python/10_summary_visualization.py
CarboxylaseDatabase/scripts/python/110_analyze_top_km.py
CarboxylaseDatabase/scripts/python/111_improved_km_training_analysis.py
CarboxylaseDatabase/scripts/python/113_full_retrain_and_discovery.py
CarboxylaseDatabase/scripts/python/114_paper_preparation.py
CarboxylaseDatabase/scripts/python/11_comprehensive_features.py
CarboxylaseDatabase/scripts/python/12_classification_model.py
CarboxylaseDatabase/scripts/python/13_predict_improvements.py
CarboxylaseDatabase/scripts/python/15_within_ec_analysis.py
CarboxylaseDatabase/scripts/python/16_correct_summary.py
CarboxylaseDatabase/scripts/python/17_improved_features_analysis.py
CarboxylaseDatabase/scripts/python/20_esm2_extraction.py
CarboxylaseDatabase/scripts/python/21_domain_analysis.py
CarboxylaseDatabase/scripts/python/22_combined_esm2_dipeptides.py
CarboxylaseDatabase/scripts/python/24_extended_analysis.py
CarboxylaseDatabase/scripts/python/24_extended_analysis_v2.py
CarboxylaseDatabase/scripts/python/25_alphafold_download.py
CarboxylaseDatabase/scripts/python/25_alphafold_features.py
CarboxylaseDatabase/scripts/python/26_final_report.py
CarboxylaseDatabase/scripts/python/27_alphafold_complete_pipeline.py
CarboxylaseDatabase/scripts/python/31_comprehensive_km_learning.py
CarboxylaseDatabase/scripts/python/31_unbiased_learning_strategy.py
CarboxylaseDatabase/scripts/python/32_esmfold_api_pipeline.py
CarboxylaseDatabase/scripts/python/32_esmfold_background_pipeline.py
CarboxylaseDatabase/scripts/python/33_analyze_3d_features.py
CarboxylaseDatabase/scripts/python/33_esmfold_gpu_pipeline.py
CarboxylaseDatabase/scripts/python/34_analyze_all_features.py
CarboxylaseDatabase/scripts/python/35_extract_3d_features_now.py
CarboxylaseDatabase/scripts/python/36_3d_within_ec_analysis.py
CarboxylaseDatabase/scripts/python/37_advanced_feature_strategy.py
CarboxylaseDatabase/scripts/python/38_active_site_features.py
CarboxylaseDatabase/scripts/python/40_advanced_km_prediction_framework.py
CarboxylaseDatabase/scripts/python/41_final_conclusions.py
CarboxylaseDatabase/scripts/python/42_hierarchical_km_model.py
CarboxylaseDatabase/scripts/python/43_detailed_feature_explanation.py
CarboxylaseDatabase/scripts/python/44_design_optimized_rubisco.py
CarboxylaseDatabase/scripts/python/45_fetch_best_rubisco.py
CarboxylaseDatabase/scripts/python/46_create_rubisco_mutants.py
CarboxylaseDatabase/scripts/python/47_predict_mutant_km.py
CarboxylaseDatabase/scripts/python/48_optimize_zoysia_template.py
CarboxylaseDatabase/scripts/python/49_fold_mutants_esmfold.py
CarboxylaseDatabase/scripts/python/50_fix_pdb_for_docking.py
CarboxylaseDatabase/scripts/python/51_add_mg_to_pdb.py
CarboxylaseDatabase/scripts/python/52_create_control_sequences.py
CarboxylaseDatabase/scripts/python/53_fetch_control_sequences.py
CarboxylaseDatabase/scripts/python/54_fold_controls.py
CarboxylaseDatabase/scripts/python/55_find_good_controls.py
CarboxylaseDatabase/scripts/python/56_fold_controls_v2.py
CarboxylaseDatabase/scripts/python/57_debug_and_fix_controls.py
CarboxylaseDatabase/scripts/python/58_find_rubisco_controls.py
CarboxylaseDatabase/scripts/python/59_fold_rubisco_controls.py
CarboxylaseDatabase/scripts/python/60_find_medium_high_km_rubisco.py
CarboxylaseDatabase/scripts/python/61_find_moderate_km_rubisco.py
CarboxylaseDatabase/scripts/python/63_prepare_complete_controls.py
CarboxylaseDatabase/scripts/python/70_test_classifier_negatives.py
--- big files (>500 MB)
2.2 GB	CarboxylaseDatabase/data/raw_data/all_sequences_20251214_211952.fasta
1.8 GB	CarboxylaseDatabase/data/raw_data/negative_non_carboxylases_20251214_211952.fasta
1.7 GB	CarboxylaseDatabase/data/processed_data/negative_cleaned_20251214_220554.fasta
0.7 GB	CarboxylaseDatabase/data/processed_data/negative_clustered_20251214_220554.fasta
0.6 GB	CarboxylaseDatabase/data/raw_data/brenda_raw_20251214_204820.json
--- databases / models
CarboxylaseDatabase/models/km_prediction_pipeline.pkl
CarboxylaseDatabase/models/hierarchical_km_model.pkl
CarboxylaseDatabase/models/binary_carboxylase_classifier.pkl
CarboxylaseDatabase/models/combined_domain_dipeptide_classifier.pkl
CarboxylaseDatabase/models/domain_classifier.pkl
CarboxylaseDatabase/models/mlp_final_classifier.pkl
CarboxylaseDatabase/models/km_final_predictor.pkl
CarboxylaseDatabase/models/km_thorough_best.pkl
CarboxylaseDatabase/models/km_with_domains.pkl
CarboxylaseDatabase/models/ec_classifier.pkl
CarboxylaseDatabase/models_v2/binary_classifier.pkl
CarboxylaseDatabase/models_v2/ec_classifier.pkl
CarboxylaseDatabase/models_v2/km_predictor.pkl
--- download provenance
CarboxylaseDatabase/scripts/python/01_brenda_download.py
CarboxylaseDatabase/scripts/python/06_advanced_features.py
CarboxylaseDatabase/scripts/python/11_comprehensive_features.py
CarboxylaseDatabase/scripts/python/25_alphafold_download.py
CarboxylaseDatabase/scripts/python/25_alphafold_features.py
CarboxylaseDatabase/scripts/python/27_alphafold_complete_pipeline.py
CarboxylaseDatabase/scripts/python/44_design_optimized_rubisco.py
CarboxylaseDatabase/scripts/python/45_fetch_best_rubisco.py
CarboxylaseDatabase/scripts/python/48_optimize_zoysia_template.py
CarboxylaseDatabase/scripts/python/53_fetch_control_sequences.py
CarboxylaseDatabase/scripts/python/76_download_diverse_negatives.py
CarboxylaseDatabase/scripts/python/84_test_real_domain_annotations.py
CarboxylaseDatabase/scripts/python/87_find_new_carboxylases.py
CarboxylaseDatabase/scripts/python/88_validate_km_predictions.py
CarboxylaseDatabase/scripts/python/92_real_pfam_interpro_test.py
```

## carboxylase_project
```
size: 20K   owner: job37yv   modified: 2025-10-08
--- tree (depth 2)
carboxylase_project
carboxylase_project/data
carboxylase_project/data/models
carboxylase_project/data/processed
carboxylase_project/data/raw
carboxylase_project/documents
carboxylase_project/documents/workflow_overview.md
carboxylase_project/environment.yml
carboxylase_project/main.nf
carboxylase_project/modules
carboxylase_project/nextflow.config
carboxylase_project/notebooks
carboxylase_project/README.md
carboxylase_project/results
carboxylase_project/results/features
carboxylase_project/results/logs
carboxylase_project/results/motifs
carboxylase_project/results/scores
carboxylase_project/scripts
carboxylase_project/scripts/carboxylase_pipeline
carboxylase_project/scripts/create_env.sh
carboxylase_project/scripts/environment_ai.yml
carboxylase_project/scripts/environment.yml
carboxylase_project/scripts/setup_carboxylase_project.sh
carboxylase_project/src
carboxylase_project/src/annotate
carboxylase_project/src/harvest
carboxylase_project/src/ml
carboxylase_project/src/persist
carboxylase_project/src/utils
carboxylase_project/webapp
--- scripts
carboxylase_project/scripts/create_env.sh
carboxylase_project/scripts/setup_carboxylase_project.sh
--- big files (>500 MB)
--- databases / models
--- download provenance
carboxylase_project/scripts/setup_carboxylase_project.sh
```

## Carboxylase_Server.git
```
size: 8.9M   owner: job37yv   modified: 2026-03-31
--- tree (depth 2)
Carboxylase_Server.git
Carboxylase_Server.git/branches
Carboxylase_Server.git/config
Carboxylase_Server.git/description
Carboxylase_Server.git/HEAD
Carboxylase_Server.git/hooks
Carboxylase_Server.git/hooks/applypatch-msg.sample
Carboxylase_Server.git/hooks/commit-msg.sample
Carboxylase_Server.git/hooks/fsmonitor-watchman.sample
Carboxylase_Server.git/hooks/post-update.sample
Carboxylase_Server.git/hooks/pre-applypatch.sample
Carboxylase_Server.git/hooks/pre-commit.sample
Carboxylase_Server.git/hooks/pre-merge-commit.sample
Carboxylase_Server.git/hooks/prepare-commit-msg.sample
Carboxylase_Server.git/hooks/pre-push.sample
Carboxylase_Server.git/hooks/pre-rebase.sample
Carboxylase_Server.git/hooks/pre-receive.sample
Carboxylase_Server.git/hooks/push-to-checkout.sample
Carboxylase_Server.git/hooks/update.sample
Carboxylase_Server.git/info
Carboxylase_Server.git/info/exclude
Carboxylase_Server.git/objects
Carboxylase_Server.git/objects/info
Carboxylase_Server.git/objects/pack
Carboxylase_Server.git/packed-refs
Carboxylase_Server.git/refs
Carboxylase_Server.git/refs/heads
Carboxylase_Server.git/refs/tags
--- scripts
--- big files (>500 MB)
--- databases / models
--- download provenance
```

## carbodb.sqlite contents (read-only)
db: CarboDB_v3/data/primary/carbodb.sqlite
```
id_map: 2380446
sqlite_sequence: 11
sequences: 2380446
ec_evidence: 2380446
km_evidence: 2971
features_composition: 2379552
features_domains: 2379446
features_interpro: 2380446
features_expert_motifs: 0
features_esm2: 2379552
features_ankh: 2380446
features_blast: 0
predictions: 2380446
confidence_scores: 2380446
db_metadata: 21
features_fimo: 0
external_annotations_cache: 3
predictions_v5_count|2380446
predictions_v5_date|20260408_162114
schema_version|2.1
km_unit|mM
created|2026-04-09 10:41:25
label_system|0=negative|1=carboxylase|2=ancestral_CO2
model_current|v6
esm2_dim|1280
ankh_dim|1536
motif_status|pending_expert_review
n_sequences|2380446
n_label1|503841
n_label2|121512
n_label0|1755093
n_km_gold|2971
n_features_comp|2379552
n_features_domains|2379446
n_features_interpro|2380446
n_features_esm2|2379552
n_features_ankh|2380446
```

## FULL SQL SCHEMA — /storage/users/job37yv/Projects/CarboDB_v3/data/primary/carbodb.sqlite
```sql
CREATE TABLE id_map (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    cdb_id     TEXT UNIQUE NOT NULL,   -- CDB000001 format
    uniprot_id TEXT UNIQUE NOT NULL
);
CREATE TABLE sqlite_sequence(name,seq);
CREATE INDEX idx_idmap_cdb ON id_map(cdb_id);
CREATE INDEX idx_idmap_uid ON id_map(uniprot_id);
CREATE TABLE sequences (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    cdb_id       TEXT    UNIQUE NOT NULL,
    uniprot_id   TEXT    UNIQUE NOT NULL,
    ec_number    TEXT    NOT NULL,
    label        INTEGER NOT NULL,
    source       TEXT    NOT NULL,
    sequence     TEXT    NOT NULL,
    length       INTEGER NOT NULL,
    organism     TEXT,
    reviewed     INTEGER DEFAULT 0,      -- 1=SwissProt curated
    km_best_mM   REAL,                   -- best experimental Km (mM); NULL if none
    km_log10_mM  REAL,
    seq_valid    INTEGER DEFAULT 1,      -- 0=failed AA validation
    created_at   TEXT DEFAULT (datetime('now'))
);
CREATE INDEX idx_seq_uid    ON sequences(uniprot_id);
CREATE INDEX idx_seq_cdb    ON sequences(cdb_id);
CREATE INDEX idx_seq_label  ON sequences(label);
CREATE INDEX idx_seq_ec     ON sequences(ec_number);
CREATE INDEX idx_seq_source ON sequences(source);
CREATE TABLE ec_evidence (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id    TEXT    NOT NULL,
    ec_number     TEXT    NOT NULL,
    source        TEXT    NOT NULL,
    evidence_tier INTEGER NOT NULL,
    confidence    REAL,               -- model probability if predicted
    model_version TEXT,               -- v6 for ML predictions
    created_at    TEXT DEFAULT (datetime('now'))
);
CREATE INDEX idx_ec_uid  ON ec_evidence(uniprot_id);
CREATE INDEX idx_ec_num  ON ec_evidence(ec_number);
CREATE INDEX idx_ec_tier ON ec_evidence(evidence_tier);
CREATE TABLE km_evidence (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id    TEXT    NOT NULL,
    ec_number     TEXT,
    km_value_mM   REAL    NOT NULL,
    km_log10_mM   REAL,
    km_unit       TEXT    DEFAULT 'mM',
    substrate     TEXT,               -- CO2|HCO3-|bicarbonate
    source        TEXT    NOT NULL,   -- brenda|model_v6
    evidence_tier INTEGER NOT NULL,
    commentary    TEXT,
    model_version TEXT,
    created_at    TEXT DEFAULT (datetime('now'))
);
CREATE INDEX idx_km_uid  ON km_evidence(uniprot_id);
CREATE INDEX idx_km_tier ON km_evidence(evidence_tier);
CREATE INDEX idx_km_val  ON km_evidence(km_value_mM);
CREATE VIEW best_evidence AS
SELECT
    s.cdb_id, s.uniprot_id, s.ec_number, s.organism, s.length,
    s.label, s.reviewed, s.source,
    ec.ec_number    AS ec_best,
    ec.evidence_tier AS ec_tier,
    ec.source        AS ec_source,
    km.km_value_mM,
    km.km_log10_mM,
    km.evidence_tier AS km_tier,
    km.source        AS km_source
FROM sequences s
LEFT JOIN ec_evidence ec ON ec.uniprot_id = s.uniprot_id
    AND ec.evidence_tier = (SELECT MIN(e2.evidence_tier) FROM ec_evidence e2 WHERE e2.uniprot_id = s.uniprot_id)
LEFT JOIN km_evidence km ON km.uniprot_id = s.uniprot_id
    AND km.evidence_tier = (SELECT MIN(k2.evidence_tier) FROM km_evidence k2 WHERE k2.uniprot_id = s.uniprot_id)
/* best_evidence(cdb_id,uniprot_id,ec_number,organism,length,label,reviewed,source,ec_best,ec_tier,ec_source,km_value_mM,km_log10_mM,km_tier,km_source) */;
CREATE TABLE features_composition (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    -- A1: AA composition (20)
    aac_A REAL, aac_C REAL, aac_D REAL, aac_E REAL, aac_F REAL,
    aac_G REAL, aac_H REAL, aac_I REAL, aac_K REAL, aac_L REAL,
    aac_M REAL, aac_N REAL, aac_P REAL, aac_Q REAL, aac_R REAL,
    aac_S REAL, aac_T REAL, aac_V REAL, aac_W REAL, aac_Y REAL,
    -- A4: Physicochemical (15)
    phys_length REAL, phys_length_log REAL, phys_mw REAL, phys_pi REAL,
    phys_charge_ph7 REAL, phys_gravy REAL, phys_aromaticity REAL,
    phys_instability REAL, phys_frac_glycine REAL, phys_frac_proline REAL,
    phys_frac_charged REAL, phys_frac_aromatic REAL, phys_frac_polar REAL,
    phys_frac_nonpolar REAL, phys_frac_small REAL,
    -- A5: Catalytic core (17)
    inv_cat_D REAL, inv_cat_E REAL, inv_cat_H REAL, inv_cat_K REAL,
    inv_cat_C REAL, inv_cat_S REAL, inv_cat_T REAL,
    inv_cat_mean_dist REAL, inv_cat_std_dist REAL,
    inv_cat_min_dist REAL, inv_cat_max_dist REAL, inv_cat_clustering REAL,
    inv_hydrophobic REAL, inv_charged REAL, inv_polar REAL,
    inv_aromatic REAL, inv_net_charge REAL,
    -- A6: EC-specific regex motifs (7)
    motif_rubisco_kk REAL, motif_rubisco_gk REAL,
    motif_ca_hh REAL, motif_ca_his_cluster REAL,
    motif_pepc_rr REAL, motif_biotin_mk REAL, motif_biotin_amk REAL,
    -- A2+A3: JSON blobs
    dp_AA REAL DEFAULT 0.0,
    dp_AC REAL DEFAULT 0.0,
    dp_AD REAL DEFAULT 0.0,
    dp_AE REAL DEFAULT 0.0,
    dp_AF REAL DEFAULT 0.0,
    dp_AG REAL DEFAULT 0.0,
    dp_AH REAL DEFAULT 0.0,
    dp_AI REAL DEFAULT 0.0,
    dp_AK REAL DEFAULT 0.0,
    dp_AL REAL DEFAULT 0.0,
    dp_AM REAL DEFAULT 0.0,
    dp_AN REAL DEFAULT 0.0,
    dp_AP REAL DEFAULT 0.0,
    dp_AQ REAL DEFAULT 0.0,
    dp_AR REAL DEFAULT 0.0,
    dp_AS REAL DEFAULT 0.0,
    dp_AT REAL DEFAULT 0.0,
    dp_AV REAL DEFAULT 0.0,
    dp_AW REAL DEFAULT 0.0,
    dp_AY REAL DEFAULT 0.0,
    dp_CA REAL DEFAULT 0.0,
    dp_CC REAL DEFAULT 0.0,
    dp_CD REAL DEFAULT 0.0,
    dp_CE REAL DEFAULT 0.0,
    dp_CF REAL DEFAULT 0.0,
    dp_CG REAL DEFAULT 0.0,
    dp_CH REAL DEFAULT 0.0,
    dp_CI REAL DEFAULT 0.0,
    dp_CK REAL DEFAULT 0.0,
    dp_CL REAL DEFAULT 0.0,
    dp_CM REAL DEFAULT 0.0,
    dp_CN REAL DEFAULT 0.0,
    dp_CP REAL DEFAULT 0.0,
    dp_CQ REAL DEFAULT 0.0,
    dp_CR REAL DEFAULT 0.0,
    dp_CS REAL DEFAULT 0.0,
    dp_CT REAL DEFAULT 0.0,
    dp_CV REAL DEFAULT 0.0,
    dp_CW REAL DEFAULT 0.0,
    dp_CY REAL DEFAULT 0.0,
    dp_DA REAL DEFAULT 0.0,
    dp_DC REAL DEFAULT 0.0,
    dp_DD REAL DEFAULT 0.0,
    dp_DE REAL DEFAULT 0.0,
    dp_DF REAL DEFAULT 0.0,
    dp_DG REAL DEFAULT 0.0,
    dp_DH REAL DEFAULT 0.0,
    dp_DI REAL DEFAULT 0.0,
    dp_DK REAL DEFAULT 0.0,
    dp_DL REAL DEFAULT 0.0,
    dp_DM REAL DEFAULT 0.0,
    dp_DN REAL DEFAULT 0.0,
    dp_DP REAL DEFAULT 0.0,
    dp_DQ REAL DEFAULT 0.0,
    dp_DR REAL DEFAULT 0.0,
    dp_DS REAL DEFAULT 0.0,
    dp_DT REAL DEFAULT 0.0,
    dp_DV REAL DEFAULT 0.0,
    dp_DW REAL DEFAULT 0.0,
    dp_DY REAL DEFAULT 0.0,
    dp_EA REAL DEFAULT 0.0,
    dp_EC REAL DEFAULT 0.0,
    dp_ED REAL DEFAULT 0.0,
    dp_EE REAL DEFAULT 0.0,
    dp_EF REAL DEFAULT 0.0,
    dp_EG REAL DEFAULT 0.0,
    dp_EH REAL DEFAULT 0.0,
    dp_EI REAL DEFAULT 0.0,
    dp_EK REAL DEFAULT 0.0,
    dp_EL REAL DEFAULT 0.0,
    dp_EM REAL DEFAULT 0.0,
    dp_EN REAL DEFAULT 0.0,
    dp_EP REAL DEFAULT 0.0,
    dp_EQ REAL DEFAULT 0.0,
    dp_ER REAL DEFAULT 0.0,
    dp_ES REAL DEFAULT 0.0,
    dp_ET REAL DEFAULT 0.0,
    dp_EV REAL DEFAULT 0.0,
    dp_EW REAL DEFAULT 0.0,
    dp_EY REAL DEFAULT 0.0,
    dp_FA REAL DEFAULT 0.0,
    dp_FC REAL DEFAULT 0.0,
    dp_FD REAL DEFAULT 0.0,
    dp_FE REAL DEFAULT 0.0,
    dp_FF REAL DEFAULT 0.0,
    dp_FG REAL DEFAULT 0.0,
    dp_FH REAL DEFAULT 0.0,
    dp_FI REAL DEFAULT 0.0,
    dp_FK REAL DEFAULT 0.0,
    dp_FL REAL DEFAULT 0.0,
    dp_FM REAL DEFAULT 0.0,
    dp_FN REAL DEFAULT 0.0,
    dp_FP REAL DEFAULT 0.0,
    dp_FQ REAL DEFAULT 0.0,
    dp_FR REAL DEFAULT 0.0,
    dp_FS REAL DEFAULT 0.0,
    dp_FT REAL DEFAULT 0.0,
    dp_FV REAL DEFAULT 0.0,
    dp_FW REAL DEFAULT 0.0,
    dp_FY REAL DEFAULT 0.0,
    dp_GA REAL DEFAULT 0.0,
    dp_GC REAL DEFAULT 0.0,
    dp_GD REAL DEFAULT 0.0,
    dp_GE REAL DEFAULT 0.0,
    dp_GF REAL DEFAULT 0.0,
    dp_GG REAL DEFAULT 0.0,
    dp_GH REAL DEFAULT 0.0,
    dp_GI REAL DEFAULT 0.0,
    dp_GK REAL DEFAULT 0.0,
    dp_GL REAL DEFAULT 0.0,
    dp_GM REAL DEFAULT 0.0,
    dp_GN REAL DEFAULT 0.0,
    dp_GP REAL DEFAULT 0.0,
    dp_GQ REAL DEFAULT 0.0,
    dp_GR REAL DEFAULT 0.0,
    dp_GS REAL DEFAULT 0.0,
    dp_GT REAL DEFAULT 0.0,
    dp_GV REAL DEFAULT 0.0,
    dp_GW REAL DEFAULT 0.0,
    dp_GY REAL DEFAULT 0.0,
    dp_HA REAL DEFAULT 0.0,
    dp_HC REAL DEFAULT 0.0,
    dp_HD REAL DEFAULT 0.0,
    dp_HE REAL DEFAULT 0.0,
    dp_HF REAL DEFAULT 0.0,
    dp_HG REAL DEFAULT 0.0,
    dp_HH REAL DEFAULT 0.0,
    dp_HI REAL DEFAULT 0.0,
    dp_HK REAL DEFAULT 0.0,
    dp_HL REAL DEFAULT 0.0,
    dp_HM REAL DEFAULT 0.0,
    dp_HN REAL DEFAULT 0.0,
    dp_HP REAL DEFAULT 0.0,
    dp_HQ REAL DEFAULT 0.0,
    dp_HR REAL DEFAULT 0.0,
    dp_HS REAL DEFAULT 0.0,
    dp_HT REAL DEFAULT 0.0,
    dp_HV REAL DEFAULT 0.0,
    dp_HW REAL DEFAULT 0.0,
    dp_HY REAL DEFAULT 0.0,
    dp_IA REAL DEFAULT 0.0,
    dp_IC REAL DEFAULT 0.0,
    dp_ID REAL DEFAULT 0.0,
    dp_IE REAL DEFAULT 0.0,
    dp_IF REAL DEFAULT 0.0,
    dp_IG REAL DEFAULT 0.0,
    dp_IH REAL DEFAULT 0.0,
    dp_II REAL DEFAULT 0.0,
    dp_IK REAL DEFAULT 0.0,
    dp_IL REAL DEFAULT 0.0,
    dp_IM REAL DEFAULT 0.0,
    dp_IN REAL DEFAULT 0.0,
    dp_IP REAL DEFAULT 0.0,
    dp_IQ REAL DEFAULT 0.0,
    dp_IR REAL DEFAULT 0.0,
    dp_IS REAL DEFAULT 0.0,
    dp_IT REAL DEFAULT 0.0,
    dp_IV REAL DEFAULT 0.0,
    dp_IW REAL DEFAULT 0.0,
    dp_IY REAL DEFAULT 0.0,
    dp_KA REAL DEFAULT 0.0,
    dp_KC REAL DEFAULT 0.0,
    dp_KD REAL DEFAULT 0.0,
    dp_KE REAL DEFAULT 0.0,
    dp_KF REAL DEFAULT 0.0,
    dp_KG REAL DEFAULT 0.0,
    dp_KH REAL DEFAULT 0.0,
    dp_KI REAL DEFAULT 0.0,
    dp_KK REAL DEFAULT 0.0,
    dp_KL REAL DEFAULT 0.0,
    dp_KM REAL DEFAULT 0.0,
    dp_KN REAL DEFAULT 0.0,
    dp_KP REAL DEFAULT 0.0,
    dp_KQ REAL DEFAULT 0.0,
    dp_KR REAL DEFAULT 0.0,
    dp_KS REAL DEFAULT 0.0,
    dp_KT REAL DEFAULT 0.0,
    dp_KV REAL DEFAULT 0.0,
    dp_KW REAL DEFAULT 0.0,
    dp_KY REAL DEFAULT 0.0,
    dp_LA REAL DEFAULT 0.0,
    dp_LC REAL DEFAULT 0.0,
    dp_LD REAL DEFAULT 0.0,
    dp_LE REAL DEFAULT 0.0,
    dp_LF REAL DEFAULT 0.0,
    dp_LG REAL DEFAULT 0.0,
    dp_LH REAL DEFAULT 0.0,
    dp_LI REAL DEFAULT 0.0,
    dp_LK REAL DEFAULT 0.0,
    dp_LL REAL DEFAULT 0.0,
    dp_LM REAL DEFAULT 0.0,
    dp_LN REAL DEFAULT 0.0,
    dp_LP REAL DEFAULT 0.0,
    dp_LQ REAL DEFAULT 0.0,
    dp_LR REAL DEFAULT 0.0,
    dp_LS REAL DEFAULT 0.0,
    dp_LT REAL DEFAULT 0.0,
    dp_LV REAL DEFAULT 0.0,
    dp_LW REAL DEFAULT 0.0,
    dp_LY REAL DEFAULT 0.0,
    dp_MA REAL DEFAULT 0.0,
    dp_MC REAL DEFAULT 0.0,
    dp_MD REAL DEFAULT 0.0,
    dp_ME REAL DEFAULT 0.0,
    dp_MF REAL DEFAULT 0.0,
    dp_MG REAL DEFAULT 0.0,
    dp_MH REAL DEFAULT 0.0,
    dp_MI REAL DEFAULT 0.0,
    dp_MK REAL DEFAULT 0.0,
    dp_ML REAL DEFAULT 0.0,
    dp_MM REAL DEFAULT 0.0,
    dp_MN REAL DEFAULT 0.0,
    dp_MP REAL DEFAULT 0.0,
    dp_MQ REAL DEFAULT 0.0,
    dp_MR REAL DEFAULT 0.0,
    dp_MS REAL DEFAULT 0.0,
    dp_MT REAL DEFAULT 0.0,
    dp_MV REAL DEFAULT 0.0,
    dp_MW REAL DEFAULT 0.0,
    dp_MY REAL DEFAULT 0.0,
    dp_NA REAL DEFAULT 0.0,
    dp_NC REAL DEFAULT 0.0,
    dp_ND REAL DEFAULT 0.0,
    dp_NE REAL DEFAULT 0.0,
    dp_NF REAL DEFAULT 0.0,
    dp_NG REAL DEFAULT 0.0,
    dp_NH REAL DEFAULT 0.0,
    dp_NI REAL DEFAULT 0.0,
    dp_NK REAL DEFAULT 0.0,
    dp_NL REAL DEFAULT 0.0,
    dp_NM REAL DEFAULT 0.0,
    dp_NN REAL DEFAULT 0.0,
    dp_NP REAL DEFAULT 0.0,
    dp_NQ REAL DEFAULT 0.0,
    dp_NR REAL DEFAULT 0.0,
    dp_NS REAL DEFAULT 0.0,
    dp_NT REAL DEFAULT 0.0,
    dp_NV REAL DEFAULT 0.0,
    dp_NW REAL DEFAULT 0.0,
    dp_NY REAL DEFAULT 0.0,
    dp_PA REAL DEFAULT 0.0,
    dp_PC REAL DEFAULT 0.0,
    dp_PD REAL DEFAULT 0.0,
    dp_PE REAL DEFAULT 0.0,
    dp_PF REAL DEFAULT 0.0,
    dp_PG REAL DEFAULT 0.0,
    dp_PH REAL DEFAULT 0.0,
    dp_PI REAL DEFAULT 0.0,
    dp_PK REAL DEFAULT 0.0,
    dp_PL REAL DEFAULT 0.0,
    dp_PM REAL DEFAULT 0.0,
    dp_PN REAL DEFAULT 0.0,
    dp_PP REAL DEFAULT 0.0,
    dp_PQ REAL DEFAULT 0.0,
    dp_PR REAL DEFAULT 0.0,
    dp_PS REAL DEFAULT 0.0,
    dp_PT REAL DEFAULT 0.0,
    dp_PV REAL DEFAULT 0.0,
    dp_PW REAL DEFAULT 0.0,
    dp_PY REAL DEFAULT 0.0,
    dp_QA REAL DEFAULT 0.0,
    dp_QC REAL DEFAULT 0.0,
    dp_QD REAL DEFAULT 0.0,
    dp_QE REAL DEFAULT 0.0,
    dp_QF REAL DEFAULT 0.0,
    dp_QG REAL DEFAULT 0.0,
    dp_QH REAL DEFAULT 0.0,
    dp_QI REAL DEFAULT 0.0,
    dp_QK REAL DEFAULT 0.0,
    dp_QL REAL DEFAULT 0.0,
    dp_QM REAL DEFAULT 0.0,
    dp_QN REAL DEFAULT 0.0,
    dp_QP REAL DEFAULT 0.0,
    dp_QQ REAL DEFAULT 0.0,
    dp_QR REAL DEFAULT 0.0,
    dp_QS REAL DEFAULT 0.0,
    dp_QT REAL DEFAULT 0.0,
    dp_QV REAL DEFAULT 0.0,
    dp_QW REAL DEFAULT 0.0,
    dp_QY REAL DEFAULT 0.0,
    dp_RA REAL DEFAULT 0.0,
    dp_RC REAL DEFAULT 0.0,
    dp_RD REAL DEFAULT 0.0,
    dp_RE REAL DEFAULT 0.0,
    dp_RF REAL DEFAULT 0.0,
    dp_RG REAL DEFAULT 0.0,
    dp_RH REAL DEFAULT 0.0,
    dp_RI REAL DEFAULT 0.0,
    dp_RK REAL DEFAULT 0.0,
    dp_RL REAL DEFAULT 0.0,
    dp_RM REAL DEFAULT 0.0,
    dp_RN REAL DEFAULT 0.0,
    dp_RP REAL DEFAULT 0.0,
    dp_RQ REAL DEFAULT 0.0,
    dp_RR REAL DEFAULT 0.0,
    dp_RS REAL DEFAULT 0.0,
    dp_RT REAL DEFAULT 0.0,
    dp_RV REAL DEFAULT 0.0,
    dp_RW REAL DEFAULT 0.0,
    dp_RY REAL DEFAULT 0.0,
    dp_SA REAL DEFAULT 0.0,
    dp_SC REAL DEFAULT 0.0,
    dp_SD REAL DEFAULT 0.0,
    dp_SE REAL DEFAULT 0.0,
    dp_SF REAL DEFAULT 0.0,
    dp_SG REAL DEFAULT 0.0,
    dp_SH REAL DEFAULT 0.0,
    dp_SI REAL DEFAULT 0.0,
    dp_SK REAL DEFAULT 0.0,
    dp_SL REAL DEFAULT 0.0,
    dp_SM REAL DEFAULT 0.0,
    dp_SN REAL DEFAULT 0.0,
    dp_SP REAL DEFAULT 0.0,
    dp_SQ REAL DEFAULT 0.0,
    dp_SR REAL DEFAULT 0.0,
    dp_SS REAL DEFAULT 0.0,
    dp_ST REAL DEFAULT 0.0,
    dp_SV REAL DEFAULT 0.0,
    dp_SW REAL DEFAULT 0.0,
    dp_SY REAL DEFAULT 0.0,
    dp_TA REAL DEFAULT 0.0,
    dp_TC REAL DEFAULT 0.0,
    dp_TD REAL DEFAULT 0.0,
    dp_TE REAL DEFAULT 0.0,
    dp_TF REAL DEFAULT 0.0,
    dp_TG REAL DEFAULT 0.0,
    dp_TH REAL DEFAULT 0.0,
    dp_TI REAL DEFAULT 0.0,
    dp_TK REAL DEFAULT 0.0,
    dp_TL REAL DEFAULT 0.0,
    dp_TM REAL DEFAULT 0.0,
    dp_TN REAL DEFAULT 0.0,
    dp_TP REAL DEFAULT 0.0,
    dp_TQ REAL DEFAULT 0.0,
    dp_TR REAL DEFAULT 0.0,
    dp_TS REAL DEFAULT 0.0,
    dp_TT REAL DEFAULT 0.0,
    dp_TV REAL DEFAULT 0.0,
    dp_TW REAL DEFAULT 0.0,
    dp_TY REAL DEFAULT 0.0,
    dp_VA REAL DEFAULT 0.0,
    dp_VC REAL DEFAULT 0.0,
    dp_VD REAL DEFAULT 0.0,
    dp_VE REAL DEFAULT 0.0,
    dp_VF REAL DEFAULT 0.0,
    dp_VG REAL DEFAULT 0.0,
    dp_VH REAL DEFAULT 0.0,
    dp_VI REAL DEFAULT 0.0,
    dp_VK REAL DEFAULT 0.0,
    dp_VL REAL DEFAULT 0.0,
    dp_VM REAL DEFAULT 0.0,
    dp_VN REAL DEFAULT 0.0,
    dp_VP REAL DEFAULT 0.0,
    dp_VQ REAL DEFAULT 0.0,
    dp_VR REAL DEFAULT 0.0,
    dp_VS REAL DEFAULT 0.0,
    dp_VT REAL DEFAULT 0.0,
    dp_VV REAL DEFAULT 0.0,
    dp_VW REAL DEFAULT 0.0,
    dp_VY REAL DEFAULT 0.0,
    dp_WA REAL DEFAULT 0.0,
    dp_WC REAL DEFAULT 0.0,
    dp_WD REAL DEFAULT 0.0,
    dp_WE REAL DEFAULT 0.0,
    dp_WF REAL DEFAULT 0.0,
    dp_WG REAL DEFAULT 0.0,
    dp_WH REAL DEFAULT 0.0,
    dp_WI REAL DEFAULT 0.0,
    dp_WK REAL DEFAULT 0.0,
    dp_WL REAL DEFAULT 0.0,
    dp_WM REAL DEFAULT 0.0,
    dp_WN REAL DEFAULT 0.0,
    dp_WP REAL DEFAULT 0.0,
    dp_WQ REAL DEFAULT 0.0,
    dp_WR REAL DEFAULT 0.0,
    dp_WS REAL DEFAULT 0.0,
    dp_WT REAL DEFAULT 0.0,
    dp_WV REAL DEFAULT 0.0,
    dp_WW REAL DEFAULT 0.0,
    dp_WY REAL DEFAULT 0.0,
    dp_YA REAL DEFAULT 0.0,
    dp_YC REAL DEFAULT 0.0,
    dp_YD REAL DEFAULT 0.0,
    dp_YE REAL DEFAULT 0.0,
    dp_YF REAL DEFAULT 0.0,
    dp_YG REAL DEFAULT 0.0,
    dp_YH REAL DEFAULT 0.0,
    dp_YI REAL DEFAULT 0.0,
    dp_YK REAL DEFAULT 0.0,
    dp_YL REAL DEFAULT 0.0,
    dp_YM REAL DEFAULT 0.0,
    dp_YN REAL DEFAULT 0.0,
    dp_YP REAL DEFAULT 0.0,
    dp_YQ REAL DEFAULT 0.0,
    dp_YR REAL DEFAULT 0.0,
    dp_YS REAL DEFAULT 0.0,
    dp_YT REAL DEFAULT 0.0,
    dp_YV REAL DEFAULT 0.0,
    dp_YW REAL DEFAULT 0.0,
    dp_YY REAL DEFAULT 0.0,
    pse_A REAL DEFAULT 0.0,
    pse_C REAL DEFAULT 0.0,
    pse_D REAL DEFAULT 0.0,
    pse_E REAL DEFAULT 0.0,
    pse_F REAL DEFAULT 0.0,
    pse_G REAL DEFAULT 0.0,
    pse_H REAL DEFAULT 0.0,
    pse_I REAL DEFAULT 0.0,
    pse_K REAL DEFAULT 0.0,
    pse_L REAL DEFAULT 0.0,
    pse_M REAL DEFAULT 0.0,
    pse_N REAL DEFAULT 0.0,
    pse_P REAL DEFAULT 0.0,
    pse_Q REAL DEFAULT 0.0,
    pse_R REAL DEFAULT 0.0,
    pse_S REAL DEFAULT 0.0,
    pse_T REAL DEFAULT 0.0,
    pse_V REAL DEFAULT 0.0,
    pse_W REAL DEFAULT 0.0,
    pse_Y REAL DEFAULT 0.0,
    pse_corr_1 REAL DEFAULT 0.0,
    pse_corr_2 REAL DEFAULT 0.0,
    pse_corr_3 REAL DEFAULT 0.0,
    pse_corr_4 REAL DEFAULT 0.0,
    pse_corr_5 REAL DEFAULT 0.0,
    pse_corr_6 REAL DEFAULT 0.0,
    pse_corr_7 REAL DEFAULT 0.0,
    pse_corr_8 REAL DEFAULT 0.0,
    pse_corr_9 REAL DEFAULT 0.0,
    pse_corr_10 REAL DEFAULT 0.0,
    feature_version TEXT DEFAULT 'v1'
);
CREATE INDEX idx_comp_uid ON features_composition(uniprot_id);
CREATE TABLE features_domains (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    pfam_hits_json  TEXT,               -- full JSON list of all Pfam hits
    pfam_PF00016 INTEGER DEFAULT 0,     -- RuBisCO large subunit
    pfam_PF02788 INTEGER DEFAULT 0,     -- RuBisCO small subunit
    pfam_PF00101 INTEGER DEFAULT 0,     -- RuBisCO-like
    pfam_PF00194 INTEGER DEFAULT 0,     -- Carbonic anhydrase
    pfam_PF03119 INTEGER DEFAULT 0,     -- Biotin carboxylase N-term
    pfam_PF00311 INTEGER DEFAULT 0,     -- PEPC
    pfam_PF00821 INTEGER DEFAULT 0,     -- Pyruvate carboxylase
    pfam_PF02785 INTEGER DEFAULT 0,     -- Biotin carboxylase C-term
    pfam_PF00364 INTEGER DEFAULT 0,     -- Biotin/lipoyl attachment
    pfam_PF01039 INTEGER DEFAULT 0,     -- Carbamoyl-phosphate synthase
    pfam_PF02786 INTEGER DEFAULT 0,     -- CPSase large N
    pfam_PF02787 INTEGER DEFAULT 0,     -- CPSase large C
    pfam_PF00289 INTEGER DEFAULT 0,     -- CPSase small
    pfam_PF01309 INTEGER DEFAULT 0,     -- PEPC C-terminal
    pfam_PF03599 INTEGER DEFAULT 0,     -- Acetyl-CoA carboxylase
    pfam_PF03590 INTEGER DEFAULT 0,     -- ACC central domain
    pfam_PF00384 INTEGER DEFAULT 0,     -- Molybdopterin oxidoreductase
    pfam_PF00682 INTEGER DEFAULT 0,     -- HMGL-like
    pfam_n_hits  INTEGER DEFAULT 0
);
CREATE INDEX idx_dom_uid  ON features_domains(uniprot_id);
CREATE INDEX idx_dom_pf16 ON features_domains(pfam_PF00016);
CREATE TABLE features_interpro (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    panther_family    TEXT,             -- PTHR accession
    panther_subfamily TEXT,
    n_panther         INTEGER DEFAULT 0,
    cath_superfamily  TEXT,             -- Gene3D primary hit
    gene3d_domains_json TEXT,
    n_gene3d          INTEGER DEFAULT 0,
    tigrfam_hits_json TEXT,
    n_tigrfam         INTEGER DEFAULT 0,
    superfamily_json  TEXT,
    scop_primary      TEXT,               -- primary SCOP superfamily (SSF accession)
    cdd_hits_json     TEXT,
    hamap_hits_json   TEXT,
    prosite_profiles_json TEXT,
    n_prosite_prof    INTEGER DEFAULT 0,
    prosite_patterns_json TEXT,
    n_prosite_pat     INTEGER DEFAULT 0,
    raw_ipr_json      TEXT              -- full raw output for extensibility
);
CREATE INDEX idx_ipr_uid     ON features_interpro(uniprot_id);
CREATE INDEX idx_ipr_panther ON features_interpro(panther_family);
CREATE INDEX idx_ipr_cath    ON features_interpro(cath_superfamily);
CREATE TABLE features_expert_motifs (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id           INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id            TEXT    UNIQUE NOT NULL,
    motif_hits_json       TEXT,
    n_motif_hits          INTEGER DEFAULT 0,
    motif_library_version TEXT    DEFAULT 'pending'
);
CREATE INDEX idx_motif_uid ON features_expert_motifs(uniprot_id);
CREATE TABLE features_esm2 (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    embedding_blob  BLOB,
    model_version   TEXT DEFAULT 'esm2_t33_650M_UR50D',
    computed_at     TEXT
);
CREATE INDEX idx_esm_uid ON features_esm2(uniprot_id);
CREATE TABLE features_ankh (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    embedding_blob  BLOB,
    model_version   TEXT DEFAULT 'ankh-large',
    computed_at     TEXT
);
CREATE INDEX idx_ankh_uid ON features_ankh(uniprot_id);
CREATE TABLE features_blast (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id       INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id        TEXT    UNIQUE NOT NULL,
    blast_best_pident REAL    DEFAULT 0.0,
    blast_best_evalue REAL    DEFAULT 999.0,
    blast_best_ec     TEXT    DEFAULT '',
    blast_has_hit     INTEGER DEFAULT 0
);
CREATE INDEX idx_blast_uid ON features_blast(uniprot_id);
CREATE TABLE predictions (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id    TEXT    NOT NULL,
    model_version TEXT    NOT NULL,   -- v6 | v7 | ...
    is_co2_pred   INTEGER,            -- 1=predicted carboxylase
    co2_prob      REAL,               -- binary classification probability
    ec_pred       TEXT,               -- predicted EC class
    ec_prob       REAL,               -- EC class confidence
    km_pred_mM    REAL,               -- predicted Km in mM
    km_pred_log10 REAL,
    created_at    TEXT DEFAULT (datetime('now'))
);
CREATE INDEX idx_pred_uid   ON predictions(uniprot_id);
CREATE INDEX idx_pred_model ON predictions(model_version);
CREATE INDEX idx_pred_ec    ON predictions(ec_pred);
CREATE INDEX idx_pred_km    ON predictions(km_pred_mM);
CREATE TABLE confidence_scores (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id      INTEGER UNIQUE REFERENCES sequences(id),
    uniprot_id       TEXT    UNIQUE NOT NULL,
    method_agreement INTEGER,   -- n independent methods agreeing on EC (0-4)
    ec_confidence    REAL,
    km_confidence    REAL,
    overall_score    REAL,
    confidence_label TEXT
);
CREATE INDEX idx_conf_uid   ON confidence_scores(uniprot_id);
CREATE INDEX idx_conf_label ON confidence_scores(confidence_label);
CREATE TABLE db_metadata (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE features_fimo (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id       INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id        TEXT    UNIQUE NOT NULL,
    motif_hits_json   TEXT,               -- {motif_id: q_value} all hits
    motif_binary_json TEXT,               -- {motif_id: 0/1} binary presence
    n_motif_hits      INTEGER DEFAULT 0,
    motif_db_version  TEXT    DEFAULT 'pending'  -- MEME motif DB version used
);
CREATE INDEX idx_fimo_uid ON features_fimo(uniprot_id);
CREATE INDEX idx_seq_organism ON sequences(organism);
CREATE INDEX idx_pred_seqid ON predictions(sequence_id);
CREATE INDEX idx_seq_label_ec ON sequences(label, ec_number);
CREATE INDEX idx_pred_seqid_modelver ON predictions(sequence_id, model_version);
CREATE TABLE external_annotations_cache (
                uniprot_id  TEXT PRIMARY KEY,
                data_json   TEXT NOT NULL,
                cached_at   INTEGER NOT NULL,
                fetch_status TEXT NOT NULL  -- 'ok' | 'partial' | 'error'
            );
CREATE INDEX idx_extcache_cached_at
                ON external_annotations_cache(cached_at)
        ;
```

## schema.sql files in repo (intended schema)
### /storage/users/job37yv/Projects/CarboDB_v3/scripts/schema.sql
```sql
-- CarboxyDB SQLite Schema v2.1
-- Single-file SQLite database. All Km in mM. Evidence-tiered design.
-- Model naming: v6 = current XGBoost (composition+domains+InterPro+embeddings+expert_motifs).
-- Expert motifs: provided by structural biology collaborator (pending).
-- ═══════════════════════════════════════════════════════════════════════════

PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

-- ── CDB_ID registry ──────────────────────────────────────────────────────────
-- Permanent append-only map. CDB_ID never changes or gets reassigned.
CREATE TABLE IF NOT EXISTS id_map (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    cdb_id     TEXT UNIQUE NOT NULL,   -- CDB000001 format
    uniprot_id TEXT UNIQUE NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_idmap_cdb ON id_map(cdb_id);
CREATE INDEX IF NOT EXISTS idx_idmap_uid ON id_map(uniprot_id);

-- ── Core sequence table ──────────────────────────────────────────────────────
-- One row per unique sequence. Central join key for all feature/evidence tables.
-- label: 1=carboxylase  2=ancestral_CO2  0=negative
-- source: brenda|trembl|swissprot|uniprot_neg|brenda_neg
CREATE TABLE IF NOT EXISTS sequences (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    cdb_id       TEXT    UNIQUE NOT NULL,
    uniprot_id   TEXT    UNIQUE NOT NULL,
    ec_number    TEXT    NOT NULL,
    label        INTEGER NOT NULL,
    source       TEXT    NOT NULL,
    sequence     TEXT    NOT NULL,
    length       INTEGER NOT NULL,
    organism     TEXT,
    reviewed     INTEGER DEFAULT 0,      -- 1=SwissProt curated
    km_best_mM   REAL,                   -- best experimental Km (mM); NULL if none
    km_log10_mM  REAL,
    seq_valid    INTEGER DEFAULT 1,      -- 0=failed AA validation
    created_at   TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_seq_uid    ON sequences(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_seq_cdb    ON sequences(cdb_id);
CREATE INDEX IF NOT EXISTS idx_seq_label  ON sequences(label);
CREATE INDEX IF NOT EXISTS idx_seq_ec     ON sequences(ec_number);
CREATE INDEX IF NOT EXISTS idx_seq_source ON sequences(source);

-- ── EC evidence ──────────────────────────────────────────────────────────────
-- Multiple EC assignments per sequence from different sources.
-- Tiers: 1=experimental  2=curated  3=predicted  4=inferred
-- source: brenda|swissprot|trembl|model_v6|pfam|blast
CREATE TABLE IF NOT EXISTS ec_evidence (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id    TEXT    NOT NULL,
    ec_number     TEXT    NOT NULL,
    source        TEXT    NOT NULL,
    evidence_tier INTEGER NOT NULL,
    confidence    REAL,               -- model probability if predicted
    model_version TEXT,               -- v6 for ML predictions
    created_at    TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_ec_uid  ON ec_evidence(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_ec_num  ON ec_evidence(ec_number);
CREATE INDEX IF NOT EXISTS idx_ec_tier ON ec_evidence(evidence_tier);

-- ── Km evidence ──────────────────────────────────────────────────────────────
-- All Km stored in mM. Tier 1=BRENDA experimental, 3=model_v6 predicted.
-- commentary: BRENDA condition notes (pH, temp, mutant flag, substrate).
CREATE TABLE IF NOT EXISTS km_evidence (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id    TEXT    NOT NULL,
    ec_number     TEXT,
    km_value_mM   REAL    NOT NULL,
    km_log10_mM   REAL,
    km_unit       TEXT    DEFAULT 'mM',
    substrate     TEXT,               -- CO2|HCO3-|bicarbonate
    source        TEXT    NOT NULL,   -- brenda|model_v6
    evidence_tier INTEGER NOT NULL,
    commentary    TEXT,
    model_version TEXT,
    created_at    TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_km_uid  ON km_evidence(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_km_tier ON km_evidence(evidence_tier);
CREATE INDEX IF NOT EXISTS idx_km_val  ON km_evidence(km_value_mM);

-- ── Best evidence view ───────────────────────────────────────────────────────
CREATE VIEW IF NOT EXISTS best_evidence AS
SELECT
    s.cdb_id, s.uniprot_id, s.ec_number, s.organism, s.length,
    s.label, s.reviewed, s.source,
    ec.ec_number    AS ec_best,
    ec.evidence_tier AS ec_tier,
    ec.source        AS ec_source,
    km.km_value_mM,
    km.km_log10_mM,
    km.evidence_tier AS km_tier,
    km.source        AS km_source
FROM sequences s
LEFT JOIN ec_evidence ec ON ec.uniprot_id = s.uniprot_id
    AND ec.evidence_tier = (SELECT MIN(e2.evidence_tier) FROM ec_evidence e2 WHERE e2.uniprot_id = s.uniprot_id)
LEFT JOIN km_evidence km ON km.uniprot_id = s.uniprot_id
    AND km.evidence_tier = (SELECT MIN(k2.evidence_tier) FROM km_evidence k2 WHERE k2.uniprot_id = s.uniprot_id);

-- ── Composition features (script 04a — DONE) ─────────────────────────────────
-- 489 features total. dp_* (400 dipeptides) + pse_* (30 PseAAC) as flat REAL columns.
CREATE TABLE IF NOT EXISTS features_composition (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    -- A1: AA composition (20)
    aac_A REAL, aac_C REAL, aac_D REAL, aac_E REAL, aac_F REAL,
    aac_G REAL, aac_H REAL, aac_I REAL, aac_K REAL, aac_L REAL,
    aac_M REAL, aac_N REAL, aac_P REAL, aac_Q REAL, aac_R REAL,
    aac_S REAL, aac_T REAL, aac_V REAL, aac_W REAL, aac_Y REAL,
    -- A4: Physicochemical (15)
    phys_length REAL, phys_length_log REAL, phys_mw REAL, phys_pi REAL,
    phys_charge_ph7 REAL, phys_gravy REAL, phys_aromaticity REAL,
    phys_instability REAL, phys_frac_glycine REAL, phys_frac_proline REAL,
    phys_frac_charged REAL, phys_frac_aromatic REAL, phys_frac_polar REAL,
    phys_frac_nonpolar REAL, phys_frac_small REAL,
    -- A5: Catalytic core (17)
    inv_cat_D REAL, inv_cat_E REAL, inv_cat_H REAL, inv_cat_K REAL,
    inv_cat_C REAL, inv_cat_S REAL, inv_cat_T REAL,
    inv_cat_mean_dist REAL, inv_cat_std_dist REAL,
    inv_cat_min_dist REAL, inv_cat_max_dist REAL, inv_cat_clustering REAL,
    inv_hydrophobic REAL, inv_charged REAL, inv_polar REAL,
    inv_aromatic REAL, inv_net_charge REAL,
    -- A6: EC-specific regex motifs (7)
    motif_rubisco_kk REAL, motif_rubisco_gk REAL,
    motif_ca_hh REAL, motif_ca_his_cluster REAL,
    motif_pepc_rr REAL, motif_biotin_mk REAL, motif_biotin_amk REAL,
    -- A2+A3: JSON blobs
    dp_AA REAL DEFAULT 0.0,
    dp_AC REAL DEFAULT 0.0,
    dp_AD REAL DEFAULT 0.0,
    dp_AE REAL DEFAULT 0.0,
    dp_AF REAL DEFAULT 0.0,
    dp_AG REAL DEFAULT 0.0,
    dp_AH REAL DEFAULT 0.0,
    dp_AI REAL DEFAULT 0.0,
    dp_AK REAL DEFAULT 0.0,
    dp_AL REAL DEFAULT 0.0,
    dp_AM REAL DEFAULT 0.0,
    dp_AN REAL DEFAULT 0.0,
    dp_AP REAL DEFAULT 0.0,
    dp_AQ REAL DEFAULT 0.0,
    dp_AR REAL DEFAULT 0.0,
    dp_AS REAL DEFAULT 0.0,
    dp_AT REAL DEFAULT 0.0,
    dp_AV REAL DEFAULT 0.0,
    dp_AW REAL DEFAULT 0.0,
    dp_AY REAL DEFAULT 0.0,
    dp_CA REAL DEFAULT 0.0,
    dp_CC REAL DEFAULT 0.0,
    dp_CD REAL DEFAULT 0.0,
    dp_CE REAL DEFAULT 0.0,
    dp_CF REAL DEFAULT 0.0,
    dp_CG REAL DEFAULT 0.0,
    dp_CH REAL DEFAULT 0.0,
    dp_CI REAL DEFAULT 0.0,
    dp_CK REAL DEFAULT 0.0,
    dp_CL REAL DEFAULT 0.0,
    dp_CM REAL DEFAULT 0.0,
    dp_CN REAL DEFAULT 0.0,
    dp_CP REAL DEFAULT 0.0,
    dp_CQ REAL DEFAULT 0.0,
    dp_CR REAL DEFAULT 0.0,
    dp_CS REAL DEFAULT 0.0,
    dp_CT REAL DEFAULT 0.0,
    dp_CV REAL DEFAULT 0.0,
    dp_CW REAL DEFAULT 0.0,
    dp_CY REAL DEFAULT 0.0,
    dp_DA REAL DEFAULT 0.0,
    dp_DC REAL DEFAULT 0.0,
    dp_DD REAL DEFAULT 0.0,
    dp_DE REAL DEFAULT 0.0,
    dp_DF REAL DEFAULT 0.0,
    dp_DG REAL DEFAULT 0.0,
    dp_DH REAL DEFAULT 0.0,
    dp_DI REAL DEFAULT 0.0,
    dp_DK REAL DEFAULT 0.0,
    dp_DL REAL DEFAULT 0.0,
    dp_DM REAL DEFAULT 0.0,
    dp_DN REAL DEFAULT 0.0,
    dp_DP REAL DEFAULT 0.0,
    dp_DQ REAL DEFAULT 0.0,
    dp_DR REAL DEFAULT 0.0,
    dp_DS REAL DEFAULT 0.0,
    dp_DT REAL DEFAULT 0.0,
    dp_DV REAL DEFAULT 0.0,
    dp_DW REAL DEFAULT 0.0,
    dp_DY REAL DEFAULT 0.0,
    dp_EA REAL DEFAULT 0.0,
    dp_EC REAL DEFAULT 0.0,
    dp_ED REAL DEFAULT 0.0,
    dp_EE REAL DEFAULT 0.0,
    dp_EF REAL DEFAULT 0.0,
    dp_EG REAL DEFAULT 0.0,
    dp_EH REAL DEFAULT 0.0,
    dp_EI REAL DEFAULT 0.0,
    dp_EK REAL DEFAULT 0.0,
    dp_EL REAL DEFAULT 0.0,
    dp_EM REAL DEFAULT 0.0,
    dp_EN REAL DEFAULT 0.0,
    dp_EP REAL DEFAULT 0.0,
    dp_EQ REAL DEFAULT 0.0,
    dp_ER REAL DEFAULT 0.0,
    dp_ES REAL DEFAULT 0.0,
    dp_ET REAL DEFAULT 0.0,
    dp_EV REAL DEFAULT 0.0,
    dp_EW REAL DEFAULT 0.0,
    dp_EY REAL DEFAULT 0.0,
    dp_FA REAL DEFAULT 0.0,
    dp_FC REAL DEFAULT 0.0,
    dp_FD REAL DEFAULT 0.0,
    dp_FE REAL DEFAULT 0.0,
    dp_FF REAL DEFAULT 0.0,
    dp_FG REAL DEFAULT 0.0,
    dp_FH REAL DEFAULT 0.0,
    dp_FI REAL DEFAULT 0.0,
    dp_FK REAL DEFAULT 0.0,
    dp_FL REAL DEFAULT 0.0,
    dp_FM REAL DEFAULT 0.0,
    dp_FN REAL DEFAULT 0.0,
    dp_FP REAL DEFAULT 0.0,
    dp_FQ REAL DEFAULT 0.0,
    dp_FR REAL DEFAULT 0.0,
    dp_FS REAL DEFAULT 0.0,
    dp_FT REAL DEFAULT 0.0,
    dp_FV REAL DEFAULT 0.0,
    dp_FW REAL DEFAULT 0.0,
    dp_FY REAL DEFAULT 0.0,
    dp_GA REAL DEFAULT 0.0,
    dp_GC REAL DEFAULT 0.0,
    dp_GD REAL DEFAULT 0.0,
    dp_GE REAL DEFAULT 0.0,
    dp_GF REAL DEFAULT 0.0,
    dp_GG REAL DEFAULT 0.0,
    dp_GH REAL DEFAULT 0.0,
    dp_GI REAL DEFAULT 0.0,
    dp_GK REAL DEFAULT 0.0,
    dp_GL REAL DEFAULT 0.0,
    dp_GM REAL DEFAULT 0.0,
    dp_GN REAL DEFAULT 0.0,
    dp_GP REAL DEFAULT 0.0,
    dp_GQ REAL DEFAULT 0.0,
    dp_GR REAL DEFAULT 0.0,
    dp_GS REAL DEFAULT 0.0,
    dp_GT REAL DEFAULT 0.0,
    dp_GV REAL DEFAULT 0.0,
    dp_GW REAL DEFAULT 0.0,
    dp_GY REAL DEFAULT 0.0,
    dp_HA REAL DEFAULT 0.0,
    dp_HC REAL DEFAULT 0.0,
    dp_HD REAL DEFAULT 0.0,
    dp_HE REAL DEFAULT 0.0,
    dp_HF REAL DEFAULT 0.0,
    dp_HG REAL DEFAULT 0.0,
    dp_HH REAL DEFAULT 0.0,
    dp_HI REAL DEFAULT 0.0,
    dp_HK REAL DEFAULT 0.0,
    dp_HL REAL DEFAULT 0.0,
    dp_HM REAL DEFAULT 0.0,
    dp_HN REAL DEFAULT 0.0,
    dp_HP REAL DEFAULT 0.0,
    dp_HQ REAL DEFAULT 0.0,
    dp_HR REAL DEFAULT 0.0,
    dp_HS REAL DEFAULT 0.0,
    dp_HT REAL DEFAULT 0.0,
    dp_HV REAL DEFAULT 0.0,
    dp_HW REAL DEFAULT 0.0,
    dp_HY REAL DEFAULT 0.0,
    dp_IA REAL DEFAULT 0.0,
    dp_IC REAL DEFAULT 0.0,
    dp_ID REAL DEFAULT 0.0,
    dp_IE REAL DEFAULT 0.0,
    dp_IF REAL DEFAULT 0.0,
    dp_IG REAL DEFAULT 0.0,
    dp_IH REAL DEFAULT 0.0,
    dp_II REAL DEFAULT 0.0,
    dp_IK REAL DEFAULT 0.0,
    dp_IL REAL DEFAULT 0.0,
    dp_IM REAL DEFAULT 0.0,
    dp_IN REAL DEFAULT 0.0,
    dp_IP REAL DEFAULT 0.0,
    dp_IQ REAL DEFAULT 0.0,
    dp_IR REAL DEFAULT 0.0,
    dp_IS REAL DEFAULT 0.0,
    dp_IT REAL DEFAULT 0.0,
    dp_IV REAL DEFAULT 0.0,
    dp_IW REAL DEFAULT 0.0,
    dp_IY REAL DEFAULT 0.0,
    dp_KA REAL DEFAULT 0.0,
    dp_KC REAL DEFAULT 0.0,
    dp_KD REAL DEFAULT 0.0,
    dp_KE REAL DEFAULT 0.0,
    dp_KF REAL DEFAULT 0.0,
    dp_KG REAL DEFAULT 0.0,
    dp_KH REAL DEFAULT 0.0,
    dp_KI REAL DEFAULT 0.0,
    dp_KK REAL DEFAULT 0.0,
    dp_KL REAL DEFAULT 0.0,
    dp_KM REAL DEFAULT 0.0,
    dp_KN REAL DEFAULT 0.0,
    dp_KP REAL DEFAULT 0.0,
    dp_KQ REAL DEFAULT 0.0,
    dp_KR REAL DEFAULT 0.0,
    dp_KS REAL DEFAULT 0.0,
    dp_KT REAL DEFAULT 0.0,
    dp_KV REAL DEFAULT 0.0,
    dp_KW REAL DEFAULT 0.0,
    dp_KY REAL DEFAULT 0.0,
    dp_LA REAL DEFAULT 0.0,
    dp_LC REAL DEFAULT 0.0,
    dp_LD REAL DEFAULT 0.0,
    dp_LE REAL DEFAULT 0.0,
    dp_LF REAL DEFAULT 0.0,
    dp_LG REAL DEFAULT 0.0,
    dp_LH REAL DEFAULT 0.0,
    dp_LI REAL DEFAULT 0.0,
    dp_LK REAL DEFAULT 0.0,
    dp_LL REAL DEFAULT 0.0,
    dp_LM REAL DEFAULT 0.0,
    dp_LN REAL DEFAULT 0.0,
    dp_LP REAL DEFAULT 0.0,
    dp_LQ REAL DEFAULT 0.0,
    dp_LR REAL DEFAULT 0.0,
    dp_LS REAL DEFAULT 0.0,
    dp_LT REAL DEFAULT 0.0,
    dp_LV REAL DEFAULT 0.0,
    dp_LW REAL DEFAULT 0.0,
    dp_LY REAL DEFAULT 0.0,
    dp_MA REAL DEFAULT 0.0,
    dp_MC REAL DEFAULT 0.0,
    dp_MD REAL DEFAULT 0.0,
    dp_ME REAL DEFAULT 0.0,
    dp_MF REAL DEFAULT 0.0,
    dp_MG REAL DEFAULT 0.0,
    dp_MH REAL DEFAULT 0.0,
    dp_MI REAL DEFAULT 0.0,
    dp_MK REAL DEFAULT 0.0,
    dp_ML REAL DEFAULT 0.0,
    dp_MM REAL DEFAULT 0.0,
    dp_MN REAL DEFAULT 0.0,
    dp_MP REAL DEFAULT 0.0,
    dp_MQ REAL DEFAULT 0.0,
    dp_MR REAL DEFAULT 0.0,
    dp_MS REAL DEFAULT 0.0,
    dp_MT REAL DEFAULT 0.0,
    dp_MV REAL DEFAULT 0.0,
    dp_MW REAL DEFAULT 0.0,
    dp_MY REAL DEFAULT 0.0,
    dp_NA REAL DEFAULT 0.0,
    dp_NC REAL DEFAULT 0.0,
    dp_ND REAL DEFAULT 0.0,
    dp_NE REAL DEFAULT 0.0,
    dp_NF REAL DEFAULT 0.0,
    dp_NG REAL DEFAULT 0.0,
    dp_NH REAL DEFAULT 0.0,
    dp_NI REAL DEFAULT 0.0,
    dp_NK REAL DEFAULT 0.0,
    dp_NL REAL DEFAULT 0.0,
    dp_NM REAL DEFAULT 0.0,
    dp_NN REAL DEFAULT 0.0,
    dp_NP REAL DEFAULT 0.0,
    dp_NQ REAL DEFAULT 0.0,
    dp_NR REAL DEFAULT 0.0,
    dp_NS REAL DEFAULT 0.0,
    dp_NT REAL DEFAULT 0.0,
    dp_NV REAL DEFAULT 0.0,
    dp_NW REAL DEFAULT 0.0,
    dp_NY REAL DEFAULT 0.0,
    dp_PA REAL DEFAULT 0.0,
    dp_PC REAL DEFAULT 0.0,
    dp_PD REAL DEFAULT 0.0,
    dp_PE REAL DEFAULT 0.0,
    dp_PF REAL DEFAULT 0.0,
    dp_PG REAL DEFAULT 0.0,
    dp_PH REAL DEFAULT 0.0,
    dp_PI REAL DEFAULT 0.0,
    dp_PK REAL DEFAULT 0.0,
    dp_PL REAL DEFAULT 0.0,
    dp_PM REAL DEFAULT 0.0,
    dp_PN REAL DEFAULT 0.0,
    dp_PP REAL DEFAULT 0.0,
    dp_PQ REAL DEFAULT 0.0,
    dp_PR REAL DEFAULT 0.0,
    dp_PS REAL DEFAULT 0.0,
    dp_PT REAL DEFAULT 0.0,
    dp_PV REAL DEFAULT 0.0,
    dp_PW REAL DEFAULT 0.0,
    dp_PY REAL DEFAULT 0.0,
    dp_QA REAL DEFAULT 0.0,
    dp_QC REAL DEFAULT 0.0,
    dp_QD REAL DEFAULT 0.0,
    dp_QE REAL DEFAULT 0.0,
    dp_QF REAL DEFAULT 0.0,
    dp_QG REAL DEFAULT 0.0,
    dp_QH REAL DEFAULT 0.0,
    dp_QI REAL DEFAULT 0.0,
    dp_QK REAL DEFAULT 0.0,
    dp_QL REAL DEFAULT 0.0,
    dp_QM REAL DEFAULT 0.0,
    dp_QN REAL DEFAULT 0.0,
    dp_QP REAL DEFAULT 0.0,
    dp_QQ REAL DEFAULT 0.0,
    dp_QR REAL DEFAULT 0.0,
    dp_QS REAL DEFAULT 0.0,
    dp_QT REAL DEFAULT 0.0,
    dp_QV REAL DEFAULT 0.0,
    dp_QW REAL DEFAULT 0.0,
    dp_QY REAL DEFAULT 0.0,
    dp_RA REAL DEFAULT 0.0,
    dp_RC REAL DEFAULT 0.0,
    dp_RD REAL DEFAULT 0.0,
    dp_RE REAL DEFAULT 0.0,
    dp_RF REAL DEFAULT 0.0,
    dp_RG REAL DEFAULT 0.0,
    dp_RH REAL DEFAULT 0.0,
    dp_RI REAL DEFAULT 0.0,
    dp_RK REAL DEFAULT 0.0,
    dp_RL REAL DEFAULT 0.0,
    dp_RM REAL DEFAULT 0.0,
    dp_RN REAL DEFAULT 0.0,
    dp_RP REAL DEFAULT 0.0,
    dp_RQ REAL DEFAULT 0.0,
    dp_RR REAL DEFAULT 0.0,
    dp_RS REAL DEFAULT 0.0,
    dp_RT REAL DEFAULT 0.0,
    dp_RV REAL DEFAULT 0.0,
    dp_RW REAL DEFAULT 0.0,
    dp_RY REAL DEFAULT 0.0,
    dp_SA REAL DEFAULT 0.0,
    dp_SC REAL DEFAULT 0.0,
    dp_SD REAL DEFAULT 0.0,
    dp_SE REAL DEFAULT 0.0,
    dp_SF REAL DEFAULT 0.0,
    dp_SG REAL DEFAULT 0.0,
    dp_SH REAL DEFAULT 0.0,
    dp_SI REAL DEFAULT 0.0,
    dp_SK REAL DEFAULT 0.0,
    dp_SL REAL DEFAULT 0.0,
    dp_SM REAL DEFAULT 0.0,
    dp_SN REAL DEFAULT 0.0,
    dp_SP REAL DEFAULT 0.0,
    dp_SQ REAL DEFAULT 0.0,
    dp_SR REAL DEFAULT 0.0,
    dp_SS REAL DEFAULT 0.0,
    dp_ST REAL DEFAULT 0.0,
    dp_SV REAL DEFAULT 0.0,
    dp_SW REAL DEFAULT 0.0,
    dp_SY REAL DEFAULT 0.0,
    dp_TA REAL DEFAULT 0.0,
    dp_TC REAL DEFAULT 0.0,
    dp_TD REAL DEFAULT 0.0,
    dp_TE REAL DEFAULT 0.0,
    dp_TF REAL DEFAULT 0.0,
    dp_TG REAL DEFAULT 0.0,
    dp_TH REAL DEFAULT 0.0,
    dp_TI REAL DEFAULT 0.0,
    dp_TK REAL DEFAULT 0.0,
    dp_TL REAL DEFAULT 0.0,
    dp_TM REAL DEFAULT 0.0,
    dp_TN REAL DEFAULT 0.0,
    dp_TP REAL DEFAULT 0.0,
    dp_TQ REAL DEFAULT 0.0,
    dp_TR REAL DEFAULT 0.0,
    dp_TS REAL DEFAULT 0.0,
    dp_TT REAL DEFAULT 0.0,
    dp_TV REAL DEFAULT 0.0,
    dp_TW REAL DEFAULT 0.0,
    dp_TY REAL DEFAULT 0.0,
    dp_VA REAL DEFAULT 0.0,
    dp_VC REAL DEFAULT 0.0,
    dp_VD REAL DEFAULT 0.0,
    dp_VE REAL DEFAULT 0.0,
    dp_VF REAL DEFAULT 0.0,
    dp_VG REAL DEFAULT 0.0,
    dp_VH REAL DEFAULT 0.0,
    dp_VI REAL DEFAULT 0.0,
    dp_VK REAL DEFAULT 0.0,
    dp_VL REAL DEFAULT 0.0,
    dp_VM REAL DEFAULT 0.0,
    dp_VN REAL DEFAULT 0.0,
    dp_VP REAL DEFAULT 0.0,
    dp_VQ REAL DEFAULT 0.0,
    dp_VR REAL DEFAULT 0.0,
    dp_VS REAL DEFAULT 0.0,
    dp_VT REAL DEFAULT 0.0,
    dp_VV REAL DEFAULT 0.0,
    dp_VW REAL DEFAULT 0.0,
    dp_VY REAL DEFAULT 0.0,
    dp_WA REAL DEFAULT 0.0,
    dp_WC REAL DEFAULT 0.0,
    dp_WD REAL DEFAULT 0.0,
    dp_WE REAL DEFAULT 0.0,
    dp_WF REAL DEFAULT 0.0,
    dp_WG REAL DEFAULT 0.0,
    dp_WH REAL DEFAULT 0.0,
    dp_WI REAL DEFAULT 0.0,
    dp_WK REAL DEFAULT 0.0,
    dp_WL REAL DEFAULT 0.0,
    dp_WM REAL DEFAULT 0.0,
    dp_WN REAL DEFAULT 0.0,
    dp_WP REAL DEFAULT 0.0,
    dp_WQ REAL DEFAULT 0.0,
    dp_WR REAL DEFAULT 0.0,
    dp_WS REAL DEFAULT 0.0,
    dp_WT REAL DEFAULT 0.0,
    dp_WV REAL DEFAULT 0.0,
    dp_WW REAL DEFAULT 0.0,
    dp_WY REAL DEFAULT 0.0,
    dp_YA REAL DEFAULT 0.0,
    dp_YC REAL DEFAULT 0.0,
    dp_YD REAL DEFAULT 0.0,
    dp_YE REAL DEFAULT 0.0,
    dp_YF REAL DEFAULT 0.0,
    dp_YG REAL DEFAULT 0.0,
    dp_YH REAL DEFAULT 0.0,
    dp_YI REAL DEFAULT 0.0,
    dp_YK REAL DEFAULT 0.0,
    dp_YL REAL DEFAULT 0.0,
    dp_YM REAL DEFAULT 0.0,
    dp_YN REAL DEFAULT 0.0,
    dp_YP REAL DEFAULT 0.0,
    dp_YQ REAL DEFAULT 0.0,
    dp_YR REAL DEFAULT 0.0,
    dp_YS REAL DEFAULT 0.0,
    dp_YT REAL DEFAULT 0.0,
    dp_YV REAL DEFAULT 0.0,
    dp_YW REAL DEFAULT 0.0,
    dp_YY REAL DEFAULT 0.0,
    pse_A REAL DEFAULT 0.0,
    pse_C REAL DEFAULT 0.0,
    pse_D REAL DEFAULT 0.0,
    pse_E REAL DEFAULT 0.0,
    pse_F REAL DEFAULT 0.0,
    pse_G REAL DEFAULT 0.0,
    pse_H REAL DEFAULT 0.0,
    pse_I REAL DEFAULT 0.0,
    pse_K REAL DEFAULT 0.0,
    pse_L REAL DEFAULT 0.0,
    pse_M REAL DEFAULT 0.0,
    pse_N REAL DEFAULT 0.0,
    pse_P REAL DEFAULT 0.0,
    pse_Q REAL DEFAULT 0.0,
    pse_R REAL DEFAULT 0.0,
    pse_S REAL DEFAULT 0.0,
    pse_T REAL DEFAULT 0.0,
    pse_V REAL DEFAULT 0.0,
    pse_W REAL DEFAULT 0.0,
    pse_Y REAL DEFAULT 0.0,
    pse_corr_1 REAL DEFAULT 0.0,
    pse_corr_2 REAL DEFAULT 0.0,
    pse_corr_3 REAL DEFAULT 0.0,
    pse_corr_4 REAL DEFAULT 0.0,
    pse_corr_5 REAL DEFAULT 0.0,
    pse_corr_6 REAL DEFAULT 0.0,
    pse_corr_7 REAL DEFAULT 0.0,
    pse_corr_8 REAL DEFAULT 0.0,
    pse_corr_9 REAL DEFAULT 0.0,
    pse_corr_10 REAL DEFAULT 0.0,
    feature_version TEXT DEFAULT 'v1'
);
CREATE INDEX IF NOT EXISTS idx_comp_uid ON features_composition(uniprot_id);

-- ── Pfam domain features (script 04b — running) ───────────────────────────────
-- HMMER vs Pfam-A.hmm. 18 carboxylase-relevant binary columns.
CREATE TABLE IF NOT EXISTS features_domains (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    pfam_hits_json  TEXT,               -- full JSON list of all Pfam hits
    pfam_PF00016 INTEGER DEFAULT 0,     -- RuBisCO large subunit
    pfam_PF02788 INTEGER DEFAULT 0,     -- RuBisCO small subunit
    pfam_PF00101 INTEGER DEFAULT 0,     -- RuBisCO-like
    pfam_PF00194 INTEGER DEFAULT 0,     -- Carbonic anhydrase
    pfam_PF03119 INTEGER DEFAULT 0,     -- Biotin carboxylase N-term
    pfam_PF00311 INTEGER DEFAULT 0,     -- PEPC
    pfam_PF00821 INTEGER DEFAULT 0,     -- Pyruvate carboxylase
    pfam_PF02785 INTEGER DEFAULT 0,     -- Biotin carboxylase C-term
    pfam_PF00364 INTEGER DEFAULT 0,     -- Biotin/lipoyl attachment
    pfam_PF01039 INTEGER DEFAULT 0,     -- Carbamoyl-phosphate synthase
    pfam_PF02786 INTEGER DEFAULT 0,     -- CPSase large N
    pfam_PF02787 INTEGER DEFAULT 0,     -- CPSase large C
    pfam_PF00289 INTEGER DEFAULT 0,     -- CPSase small
    pfam_PF01309 INTEGER DEFAULT 0,     -- PEPC C-terminal
    pfam_PF03599 INTEGER DEFAULT 0,     -- Acetyl-CoA carboxylase
    pfam_PF03590 INTEGER DEFAULT 0,     -- ACC central domain
    pfam_PF00384 INTEGER DEFAULT 0,     -- Molybdopterin oxidoreductase
    pfam_PF00682 INTEGER DEFAULT 0,     -- HMGL-like
    pfam_n_hits  INTEGER DEFAULT 0
);
CREATE INDEX IF NOT EXISTS idx_dom_uid  ON features_domains(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_dom_pf16 ON features_domains(pfam_PF00016);

-- ── InterProScan features (script 04c — running) ──────────────────────────────
-- Panther, Gene3D/CATH, TIGRFAM, SUPERFAMILY, CDD, HAMAP, ProSite.
-- JSON columns store full hit lists; scalar columns store primary hit + count for ML.
CREATE TABLE IF NOT EXISTS features_interpro (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    panther_family    TEXT,             -- PTHR accession
    panther_subfamily TEXT,
    n_panther         INTEGER DEFAULT 0,
    cath_superfamily  TEXT,             -- Gene3D primary hit
    gene3d_domains_json TEXT,
    n_gene3d          INTEGER DEFAULT 0,
    tigrfam_hits_json TEXT,
    n_tigrfam         INTEGER DEFAULT 0,
    superfamily_json  TEXT,
    scop_primary      TEXT,               -- primary SCOP superfamily (SSF accession)
    cdd_hits_json     TEXT,
    hamap_hits_json   TEXT,
    prosite_profiles_json TEXT,
    n_prosite_prof    INTEGER DEFAULT 0,
    prosite_patterns_json TEXT,
    n_prosite_pat     INTEGER DEFAULT 0,
    raw_ipr_json      TEXT              -- full raw output for extensibility
);
CREATE INDEX IF NOT EXISTS idx_ipr_uid     ON features_interpro(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_ipr_panther ON features_interpro(panther_family);
CREATE INDEX IF NOT EXISTS idx_ipr_cath    ON features_interpro(cath_superfamily);

-- ── Expert structural motifs (pending — structural biology collaborator) ───────
-- Replaces old features_meme. Motif library defined by expert, not automated.
-- motif_hits_json: {motif_id: 0/1} for all expert-defined motifs.
-- motif_library_version: bumped each time expert adds/revises motifs.
CREATE TABLE IF NOT EXISTS features_expert_motifs (
    id                    INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id           INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id            TEXT    UNIQUE NOT NULL,
    motif_hits_json       TEXT,
    n_motif_hits          INTEGER DEFAULT 0,
    motif_library_version TEXT    DEFAULT 'pending'
);
CREATE INDEX IF NOT EXISTS idx_motif_uid ON features_expert_motifs(uniprot_id);

-- ── ESM-2 embeddings (script 04e — running on A100) ──────────────────────────
-- esm2_t33_650M_UR50D. Mean-pooled over sequence. 1280-dim float32.
-- Stored as bytes blob: 1280 × 4 = 5,120 bytes per row.
CREATE TABLE IF NOT EXISTS features_esm2 (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    embedding_blob  BLOB,
    model_version   TEXT DEFAULT 'esm2_t33_650M_UR50D',
    computed_at     TEXT
);
CREATE INDEX IF NOT EXISTS idx_esm_uid ON features_esm2(uniprot_id);

-- ── Ankh embeddings (script 04d — running on CPU) ────────────────────────────
-- ankh-large. Mean-pooled over sequence. 1536-dim float32 (NOT 1024).
-- Stored as bytes blob: 1536 × 4 = 6,144 bytes per row.
CREATE TABLE IF NOT EXISTS features_ankh (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    embedding_blob  BLOB,
    model_version   TEXT DEFAULT 'ankh-large',
    computed_at     TEXT
);
CREATE INDEX IF NOT EXISTS idx_ankh_uid ON features_ankh(uniprot_id);

-- ── BLAST homology features (script 04f — pending) ───────────────────────────
-- BLASTp vs curated carboxylase reference DB. Used as ML feature + fast baseline.
CREATE TABLE IF NOT EXISTS features_blast (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id       INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id        TEXT    UNIQUE NOT NULL,
    blast_best_pident REAL    DEFAULT 0.0,
    blast_best_evalue REAL    DEFAULT 999.0,
    blast_best_ec     TEXT    DEFAULT '',
    blast_has_hit     INTEGER DEFAULT 0
);
CREATE INDEX IF NOT EXISTS idx_blast_uid ON features_blast(uniprot_id);

-- ── ML predictions ────────────────────────────────────────────────────────────
-- One row per sequence per model version. Multiple versions stored side by side.
-- model_version: v6 = current production model.
-- Future models (v7+) can be added without schema change.
CREATE TABLE IF NOT EXISTS predictions (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id    TEXT    NOT NULL,
    model_version TEXT    NOT NULL,   -- v6 | v7 | ...
    is_co2_pred   INTEGER,            -- 1=predicted carboxylase
    co2_prob      REAL,               -- binary classification probability
    ec_pred       TEXT,               -- predicted EC class
    ec_prob       REAL,               -- EC class confidence
    km_pred_mM    REAL,               -- predicted Km in mM
    km_pred_log10 REAL,
    created_at    TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_pred_uid   ON predictions(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_pred_model ON predictions(model_version);
CREATE INDEX IF NOT EXISTS idx_pred_ec    ON predictions(ec_pred);
CREATE INDEX IF NOT EXISTS idx_pred_km    ON predictions(km_pred_mM);

-- ── Confidence scores ─────────────────────────────────────────────────────────
-- Computed after all predictions stored. Used by web app for display.
-- confidence_label: high|medium|low|review
CREATE TABLE IF NOT EXISTS confidence_scores (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id      INTEGER UNIQUE REFERENCES sequences(id),
    uniprot_id       TEXT    UNIQUE NOT NULL,
    method_agreement INTEGER,   -- n independent methods agreeing on EC (0-4)
    ec_confidence    REAL,
    km_confidence    REAL,
    overall_score    REAL,
    confidence_label TEXT
);
CREATE INDEX IF NOT EXISTS idx_conf_uid   ON confidence_scores(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_conf_label ON confidence_scores(confidence_label);

-- ── Database metadata ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS db_metadata (key TEXT PRIMARY KEY, value TEXT);
INSERT OR REPLACE INTO db_metadata VALUES ('schema_version',  '2.1');
INSERT OR REPLACE INTO db_metadata VALUES ('km_unit',         'mM');
INSERT OR REPLACE INTO db_metadata VALUES ('created',         datetime('now'));
INSERT OR REPLACE INTO db_metadata VALUES ('label_system',    '0=negative|1=carboxylase|2=ancestral_CO2');
INSERT OR REPLACE INTO db_metadata VALUES ('model_current',   'v6');
INSERT OR REPLACE INTO db_metadata VALUES ('esm2_dim',        '1280');
INSERT OR REPLACE INTO db_metadata VALUES ('ankh_dim',        '1536');
INSERT OR REPLACE INTO db_metadata VALUES ('motif_status',    'pending_expert_review');
INSERT OR REPLACE INTO db_metadata VALUES ('n_sequences',     '2380446');
INSERT OR REPLACE INTO db_metadata VALUES ('n_km_gold',       '2971');

-- ── FIMO motif scan features (script 04f — to write) ─────────────────────────
-- FIMO scans all sequences against a MEME motif database built from known
-- carboxylase active sites and conserved regions. Automated, data-driven.
-- Separate from expert motifs which are structure-informed and hand-curated.
-- motif_hits_json: {motif_id: score} for all MEME/FIMO motifs found.
CREATE TABLE IF NOT EXISTS features_fimo (
    id                INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id       INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id        TEXT    UNIQUE NOT NULL,
    motif_hits_json   TEXT,               -- {motif_id: q_value} all hits
    motif_binary_json TEXT,               -- {motif_id: 0/1} binary presence
    n_motif_hits      INTEGER DEFAULT 0,
    motif_db_version  TEXT    DEFAULT 'pending'  -- MEME motif DB version used
);
CREATE INDEX IF NOT EXISTS idx_fimo_uid ON features_fimo(uniprot_id);
```

## per-table columns + sample row
### id_map
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|cdb_id|TEXT|1||0
2|uniprot_id|TEXT|1||0
id|cdb_id|uniprot_id
1|CDB000001|A0A009HG66
```
### sqlite_sequence
```
cid|name|type|notnull|dflt_value|pk
0|name||0||0
1|seq||0||0
name|seq
sequences|2380446
```
### sequences
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|cdb_id|TEXT|1||0
2|uniprot_id|TEXT|1||0
3|ec_number|TEXT|1||0
4|label|INTEGER|1||0
5|source|TEXT|1||0
6|sequence|TEXT|1||0
7|length|INTEGER|1||0
8|organism|TEXT|0||0
9|reviewed|INTEGER|0|0|0
10|km_best_mM|REAL|0||0
11|km_log10_mM|REAL|0||0
12|seq_valid|INTEGER|0|1|0
13|created_at|TEXT|0|datetime('now')|0
id|cdb_id|uniprot_id|ec_number|label|source|sequence|length|organism|reviewed|km_best_mM|km_log10_mM|seq_valid|created_at
1|CDB000001|A0A009HG66|1.15.1.1|0|brenda_neg|MTTITLPALPYGYDDLAPHISKETLEYHHDKHHNTYVVNLNNLIKGTDLEGKTLEEIIKATAGDASKAGIFNNAAQVWNHTFYWNSMKPNGGGKPTGAIAAKIDEAFGSYEKFAEEFTAAATTQFGSGWAWLVADEVNGKLSITKTSNADTPLAHGQIAVLTIDVWEHAYYIDFRNLRPKYIATFLENLVNWDYANAKLAGQPAGVEK|208|Acinetobacter baumannii (strain 1295743)|0|||1|2026-04-01 13:10:39
```
### ec_evidence
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|ec_number|TEXT|1||0
4|source|TEXT|1||0
5|evidence_tier|INTEGER|1||0
6|confidence|REAL|0||0
7|model_version|TEXT|0||0
8|created_at|TEXT|0|datetime('now')|0
id|sequence_id|uniprot_id|ec_number|source|evidence_tier|confidence|model_version|created_at
1|1|A0A009HG66|1.15.1.1|brenda_neg|2|||2026-04-01 13:10:40
```
### km_evidence
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|ec_number|TEXT|0||0
4|km_value_mM|REAL|1||0
5|km_log10_mM|REAL|0||0
6|km_unit|TEXT|0|'mM'|0
7|substrate|TEXT|0||0
8|source|TEXT|1||0
9|evidence_tier|INTEGER|1||0
10|commentary|TEXT|0||0
11|model_version|TEXT|0||0
12|created_at|TEXT|0|datetime('now')|0
id|sequence_id|uniprot_id|ec_number|km_value_mM|km_log10_mM|km_unit|substrate|source|evidence_tier|commentary|model_version|created_at
1|4562|A0A023L9N9|4.1.1.31|0.1|-1.0|mM|CO2|brenda|1|||2026-04-01 13:10:40
```
### features_composition
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|aac_A|REAL|0||0
4|aac_C|REAL|0||0
5|aac_D|REAL|0||0
6|aac_E|REAL|0||0
7|aac_F|REAL|0||0
8|aac_G|REAL|0||0
9|aac_H|REAL|0||0
10|aac_I|REAL|0||0
11|aac_K|REAL|0||0
12|aac_L|REAL|0||0
13|aac_M|REAL|0||0
14|aac_N|REAL|0||0
15|aac_P|REAL|0||0
16|aac_Q|REAL|0||0
17|aac_R|REAL|0||0
18|aac_S|REAL|0||0
19|aac_T|REAL|0||0
20|aac_V|REAL|0||0
21|aac_W|REAL|0||0
22|aac_Y|REAL|0||0
23|phys_length|REAL|0||0
24|phys_length_log|REAL|0||0
25|phys_mw|REAL|0||0
26|phys_pi|REAL|0||0
27|phys_charge_ph7|REAL|0||0
28|phys_gravy|REAL|0||0
29|phys_aromaticity|REAL|0||0
30|phys_instability|REAL|0||0
31|phys_frac_glycine|REAL|0||0
32|phys_frac_proline|REAL|0||0
33|phys_frac_charged|REAL|0||0
34|phys_frac_aromatic|REAL|0||0
35|phys_frac_polar|REAL|0||0
36|phys_frac_nonpolar|REAL|0||0
37|phys_frac_small|REAL|0||0
38|inv_cat_D|REAL|0||0
39|inv_cat_E|REAL|0||0
40|inv_cat_H|REAL|0||0
41|inv_cat_K|REAL|0||0
42|inv_cat_C|REAL|0||0
43|inv_cat_S|REAL|0||0
44|inv_cat_T|REAL|0||0
45|inv_cat_mean_dist|REAL|0||0
46|inv_cat_std_dist|REAL|0||0
47|inv_cat_min_dist|REAL|0||0
48|inv_cat_max_dist|REAL|0||0
49|inv_cat_clustering|REAL|0||0
50|inv_hydrophobic|REAL|0||0
51|inv_charged|REAL|0||0
52|inv_polar|REAL|0||0
53|inv_aromatic|REAL|0||0
54|inv_net_charge|REAL|0||0
55|motif_rubisco_kk|REAL|0||0
56|motif_rubisco_gk|REAL|0||0
57|motif_ca_hh|REAL|0||0
58|motif_ca_his_cluster|REAL|0||0
id|sequence_id|uniprot_id|aac_A|aac_C|aac_D|aac_E|aac_F|aac_G|aac_H|aac_I|aac_K|aac_L|aac_M|aac_N|aac_P|aac_Q|aac_R|aac_S|aac_T|aac_V|aac_W|aac_Y|phys_length|phys_length_log|phys_mw|phys_pi|phys_charge_ph7|phys_gravy|phys_aromaticity|phys_instability|phys_frac_glycine|phys_frac_proline|phys_frac_charged|phys_frac_aromatic|phys_frac_polar|phys_frac_nonpolar|phys_frac_small|inv_cat_D|inv_cat_E|inv_c
1|1|A0A009HG66|0.129807692307692|0.0|0.0528846153846153|0.0625|0.0384615384615384|0.0769230769230769|0.0384615384615384|0.0625|0.0721153846153846|0.0769230769230769|0.0096153846153846|0.0721153846153846|0.0384615384615384|0.0192307692307692|0.0096153846153846|0.0336538461538461|0.0865384615384615|0.0432692307692307|0.0288461538461538|0.048076923076923|208.0|5.34233425196481|22919.365|5.56558322906
```
### features_domains
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|pfam_hits_json|TEXT|0||0
4|pfam_PF00016|INTEGER|0|0|0
5|pfam_PF02788|INTEGER|0|0|0
6|pfam_PF00101|INTEGER|0|0|0
7|pfam_PF00194|INTEGER|0|0|0
8|pfam_PF03119|INTEGER|0|0|0
9|pfam_PF00311|INTEGER|0|0|0
10|pfam_PF00821|INTEGER|0|0|0
11|pfam_PF02785|INTEGER|0|0|0
12|pfam_PF00364|INTEGER|0|0|0
13|pfam_PF01039|INTEGER|0|0|0
14|pfam_PF02786|INTEGER|0|0|0
15|pfam_PF02787|INTEGER|0|0|0
16|pfam_PF00289|INTEGER|0|0|0
17|pfam_PF01309|INTEGER|0|0|0
18|pfam_PF03599|INTEGER|0|0|0
19|pfam_PF03590|INTEGER|0|0|0
20|pfam_PF00384|INTEGER|0|0|0
21|pfam_PF00682|INTEGER|0|0|0
22|pfam_n_hits|INTEGER|0|0|0
id|sequence_id|uniprot_id|pfam_hits_json|pfam_PF00016|pfam_PF02788|pfam_PF00101|pfam_PF00194|pfam_PF03119|pfam_PF00311|pfam_PF00821|pfam_PF02785|pfam_PF00364|pfam_PF01039|pfam_PF02786|pfam_PF02787|pfam_PF00289|pfam_PF01309|pfam_PF03599|pfam_PF03590|pfam_PF00384|pfam_PF00682|pfam_n_hits
1|1|A0A009HG66|["PF00081", "PF02777"]|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|2
```
### features_interpro
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|panther_family|TEXT|0||0
4|panther_subfamily|TEXT|0||0
5|n_panther|INTEGER|0|0|0
6|cath_superfamily|TEXT|0||0
7|gene3d_domains_json|TEXT|0||0
8|n_gene3d|INTEGER|0|0|0
9|tigrfam_hits_json|TEXT|0||0
10|n_tigrfam|INTEGER|0|0|0
11|superfamily_json|TEXT|0||0
12|scop_primary|TEXT|0||0
13|cdd_hits_json|TEXT|0||0
14|hamap_hits_json|TEXT|0||0
15|prosite_profiles_json|TEXT|0||0
16|n_prosite_prof|INTEGER|0|0|0
17|prosite_patterns_json|TEXT|0||0
18|n_prosite_pat|INTEGER|0|0|0
19|raw_ipr_json|TEXT|0||0
id|sequence_id|uniprot_id|panther_family|panther_subfamily|n_panther|cath_superfamily|gene3d_domains_json|n_gene3d|tigrfam_hits_json|n_tigrfam|superfamily_json|scop_primary|cdd_hits_json|hamap_hits_json|prosite_profiles_json|n_prosite_prof|prosite_patterns_json|n_prosite_pat|raw_ipr_json
1|1|A0A009HG66|PTHR42769||1|G3DSA:1.10.287.990|["G3DSA:1.10.287.990", "G3DSA:3.55.40.20"]|2|[]|0|["SSF54719", "SSF46609"]||[]|[]|[]|0|["PS00088"]|1|[{"database":"ProSitePatterns","accession":"PS00088","interpro_acc":"IPR019833"},{"database":"SUPERFAMILY","accession":"SSF54719","interpro_acc":"IPR036314"},{"database":"Pfam","accession":"PF00081","interpro_acc":"IPR019831"},{"database":"SUPERFAMILY"
```
### features_expert_motifs
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|motif_hits_json|TEXT|0||0
4|n_motif_hits|INTEGER|0|0|0
5|motif_library_version|TEXT|0|'pending'|0
```
### features_esm2
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|embedding_blob|BLOB|0||0
4|model_version|TEXT|0|'esm2_t33_650M_UR50D'|0
5|computed_at|TEXT|0||0
id|sequence_id|uniprot_id|embedding_blob|model_version|computed_at
1|1|A0A009HG66|�e���[�P�<��="x�cā��$���a��<��?=�K�����w'�=SK>��C)�=	#4=UL�;-�����=X\6�+">Tٔ��}�9t:_�3��,m=Y����ܼ�^��%?P���= w>�߅��MP��)F=��I��<=�E���GH���n�~����>h�a��w�~=��>=�|J>�=o,�A�^�v8
��h<�
ܽ׻�<����,ؽ#�J=�\y�{�=�ݪ��oB=�&<��W=��_���^>c�=����M=�`�='#p<��+�>>�<�@ּ�MF;��ѻ7�=R��=e��:Y��=�Zڽ�p��gC�G�=�]��� ����<@b=�oe<ϗ�=�V���|esm2_t33_650M_UR50D|20260401_151036
```
### features_ankh
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|embedding_blob|BLOB|0||0
4|model_version|TEXT|0|'ankh-large'|0
5|computed_at|TEXT|0||0
id|sequence_id|uniprot_id|embedding_blob|model_version|computed_at
1|1|A0A009HG66|)a<F�t���;��B;���:5���֌;�c��w��:zbp�tc�:�N�����*�;a�D;��(��F�:Ǽ;3�S;�P;+�:vF �yQe;|X;�'�ؒ�:Ʉ��<v;d�8;�
ڻ������m8=�
�7���ʻH�;:�:řS����p8�EJ�YU6�&�Ԕ�:i<w�4��A���9;[��:%;�Z;@����x;��:�ح9��=;�ND�]l^�*r:7����#�ڥ8����qd2���:S��9�;@�u;`���
���g3��ݰ;��|ankh-large|20260401_151036
```
### features_blast
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|blast_best_pident|REAL|0|0.0|0
4|blast_best_evalue|REAL|0|999.0|0
5|blast_best_ec|TEXT|0|''|0
6|blast_has_hit|INTEGER|0|0|0
```
### predictions
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|model_version|TEXT|1||0
4|is_co2_pred|INTEGER|0||0
5|co2_prob|REAL|0||0
6|ec_pred|TEXT|0||0
7|ec_prob|REAL|0||0
8|km_pred_mM|REAL|0||0
9|km_pred_log10|REAL|0||0
10|created_at|TEXT|0|datetime('now')|0
id|sequence_id|uniprot_id|model_version|is_co2_pred|co2_prob|ec_pred|ec_prob|km_pred_mM|km_pred_log10|created_at
1|1|A0A009HG66|v5|0|8.24207901928276e-08|4.2.1.1|0.909557521343231|||2026-04-08 14:22:02
```
### confidence_scores
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|0||0
2|uniprot_id|TEXT|1||0
3|method_agreement|INTEGER|0||0
4|ec_confidence|REAL|0||0
5|km_confidence|REAL|0||0
6|overall_score|REAL|0||0
7|confidence_label|TEXT|0||0
id|sequence_id|uniprot_id|method_agreement|ec_confidence|km_confidence|overall_score|confidence_label
1|1|A0A009HG66|1|0.909557521343231||0.454778790473938|review
```
### db_metadata
```
cid|name|type|notnull|dflt_value|pk
0|key|TEXT|0||1
1|value|TEXT|0||0
key|value
predictions_v5_count|2380446
```
### features_fimo
```
cid|name|type|notnull|dflt_value|pk
0|id|INTEGER|0||1
1|sequence_id|INTEGER|1||0
2|uniprot_id|TEXT|1||0
3|motif_hits_json|TEXT|0||0
4|motif_binary_json|TEXT|0||0
5|n_motif_hits|INTEGER|0|0|0
6|motif_db_version|TEXT|0|'pending'|0
```
### external_annotations_cache
```
cid|name|type|notnull|dflt_value|pk
0|uniprot_id|TEXT|0||1
1|data_json|TEXT|1||0
2|cached_at|INTEGER|1||0
3|fetch_status|TEXT|1||0
uniprot_id|data_json|cached_at|fetch_status
P00875|{"uniprot_id": "P00875", "uniprot": {"protein_name": "Ribulose bisphosphate carboxylase large chain", "gene_name": "rbcL", "organism": "Spinacia oleracea", "lineage": ["Eukaryota", "Viridiplantae", "Streptophyta", "Embryophyta", "Tracheophyta", "Spermatophyta", "Magnoliopsida", "eudicotyledons", "Gunneridae", "Pentapetalae", "Caryophyllales", "Chenopodiaceae", "Chenopodioideae", "Anserineae
```
