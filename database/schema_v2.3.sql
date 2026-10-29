-- =====================================================================
-- CarboDB release schema v2.3  (Paper 1 final database)
-- File: CarboDB/database/schema_v2.3.sql
-- Date: 2026-10-28
--
-- Builds ONE release file, carbodb_v5.2.sqlite, in a single pass.
-- Feature tables (composition, domains, interpro, esm2, ankh) are copied
-- unchanged from the v2.2 DB (CarboDB_v3/data/primary/carbodb.sqlite),
-- because the sequences did not change. Everything else is rebuilt.
--
-- Changes vs live DB (v2.2, metadata wrongly says 2.1):
--   * JOIN KEY = sequences.id (sequence_id) everywhere; uniprot_id kept
--     only as a convenience column (the old best_evidence view joined on it)
--   * sequences: + kingdom, taxon_id, cluster_id, split           (F3, TP/TN)
--   * km_evidence: + Km cleaning columns                           (F1)
--   * NEW domain_hits + feature_catalog: e-value, bit score, name  (F2)
--   * NEW motif_catalog + motif_hits, replacing the empty
--     features_expert_motifs and features_fimo tables              (F4)
--   * features_blast: filled this time                             (F3)
--   * predictions: global vs within-EC Km, model variant, EC/Km only
--     for predicted carboxylases                                   (F5, F6)
--   * NEW models, benchmark_results, benchmark_calls (TP/TN lists)
--   * NEW sources, runs: provenance                                (F8)
--   * external_annotations_cache removed (app-owned, separate DB)
-- =====================================================================
PRAGMA foreign_keys = ON;

-- ---------------------------------------------------------------------
-- 0. Provenance
-- ---------------------------------------------------------------------
CREATE TABLE sources (
    source_id      INTEGER PRIMARY KEY,
    name           TEXT NOT NULL UNIQUE,   -- brenda | uniprot_swissprot | uniprot_trembl | uniprot_neg | pfam | interpro | esm2 | ankh | gmgc | motifs_aparna
    release        TEXT NOT NULL,          -- e.g. 2026-03-27, Pfam 37.1, InterProScan 5.72-103.0
    url            TEXT,
    retrieved_on   TEXT,
    file_path      TEXT,                   -- path under CarboDB_data/raw
    sha256         TEXT,
    licence        TEXT                    -- CC BY 4.0 (BRENDA, UniProt) ...
);

CREATE TABLE runs (
    run_id         INTEGER PRIMARY KEY,
    stage          TEXT NOT NULL,          -- 03_merge | 04a_composition | 04b_hmmer | ... | 08_train | 10_predict
    tool           TEXT,
    tool_version   TEXT,
    params_json    TEXT,
    config_sha256  TEXT,
    code_commit    TEXT,                   -- git commit of CarboDB repo
    host           TEXT,                   -- wbbi206
    started_at     TEXT,
    finished_at    TEXT,
    n_rows         INTEGER,
    note           TEXT
);

-- ---------------------------------------------------------------------
-- 0b. Artifacts: big files live on disk, the DB holds a pointer + checksum
--     (embeddings as .npy, domtblout, models, matrices, figures, FASTA)
-- ---------------------------------------------------------------------
CREATE TABLE artifacts (
    artifact_id   INTEGER PRIMARY KEY,
    kind          TEXT NOT NULL CHECK (kind IN ('embedding','annotation_raw','fasta','matrix','model','prediction','benchmark','figure','other')),
    name          TEXT NOT NULL,          -- esm2_t33_650M | pfam_domtblout | binary_full_v5.2 ...
    path          TEXT NOT NULL,          -- relative to CarboDB_data/ (e.g. primary/embeddings/esm2/esm2.npy)
    row_index     TEXT,                   -- relative path of the sequence_id order file (embeddings, matrices)
    shape         TEXT,                   -- e.g. 2380446x1280
    dtype         TEXT,                   -- float32
    sha256        TEXT NOT NULL,
    size_bytes    INTEGER,
    run_id        INTEGER REFERENCES runs(run_id),
    created_at    TEXT DEFAULT (datetime('now')),
    UNIQUE (kind, name, path)
);

CREATE TABLE db_metadata (key TEXT PRIMARY KEY, value TEXT);
-- required keys: schema_version=2.3, release=v5.2, model_version=v5.2,
-- n_sequences, n_label0/1/2, n_km_use, built_at, code_commit, zenodo_doi

-- ---------------------------------------------------------------------
-- 1. Identity and sequences
-- ---------------------------------------------------------------------
CREATE TABLE id_map (                      -- append-only, never reassigned
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    cdb_id     TEXT UNIQUE NOT NULL,       -- CDB000001
    uniprot_id TEXT UNIQUE NOT NULL
);

CREATE TABLE sequences (
    id           INTEGER PRIMARY KEY,      -- = sequence_id everywhere
    cdb_id       TEXT UNIQUE NOT NULL REFERENCES id_map(cdb_id),
    uniprot_id   TEXT UNIQUE NOT NULL,
    ec_number    TEXT NOT NULL,            -- source EC (BRENDA/UniProt), not predicted
    label        INTEGER NOT NULL CHECK (label IN (0,1,2)),  -- 0 negative | 1 carboxylase | 2 ancestral/CO2-related
    source       TEXT NOT NULL,            -- brenda | swissprot | trembl | uniprot_neg | brenda_neg
    reviewed     INTEGER DEFAULT 0,
    sequence     TEXT NOT NULL,
    length       INTEGER NOT NULL,
    organism     TEXT,
    taxon_id     INTEGER,                  -- NCBI taxon            (new)
    kingdom      TEXT,                     -- bacteria | archaea | eukaryota | virus | unknown  (new, feature)
    seq_valid    INTEGER DEFAULT 1,
    cluster_id   INTEGER,                  -- CD-HIT 40 % cluster   (new)
    split        TEXT CHECK (split IN ('train','val','test','none')),  -- (new) ML split; 'none' for label 2
    km_best_mM   REAL,                     -- best 'use' Km after cleaning; NULL if none
    km_log10_mM  REAL
);
CREATE INDEX idx_seq_uid      ON sequences(uniprot_id);
CREATE INDEX idx_seq_label_ec ON sequences(label, ec_number);
CREATE INDEX idx_seq_split    ON sequences(split);
CREATE INDEX idx_seq_organism ON sequences(organism);

-- ---------------------------------------------------------------------
-- 2. Evidence (what the sources say)
-- ---------------------------------------------------------------------
CREATE TABLE ec_evidence (
    id            INTEGER PRIMARY KEY,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    ec_number     TEXT NOT NULL,
    source        TEXT NOT NULL,           -- brenda | swissprot | trembl
    evidence_tier INTEGER NOT NULL CHECK (evidence_tier BETWEEN 1 AND 4),  -- 1 experimental 2 curated 3 predicted 4 inferred
    source_id     INTEGER REFERENCES sources(source_id)
);
CREATE INDEX idx_ec_seq ON ec_evidence(sequence_id);
CREATE INDEX idx_ec_num ON ec_evidence(ec_number);

CREATE TABLE km_evidence (                 -- experimental Km only; predictions live in predictions
    id                   INTEGER PRIMARY KEY,
    sequence_id          INTEGER NOT NULL REFERENCES sequences(id),
    ec_number            TEXT NOT NULL,
    km_value_mM          REAL NOT NULL,    -- value used for training (after conversion)
    km_log10_mM          REAL,
    km_raw_value         REAL,             -- as reported in BRENDA
    km_raw_unit          TEXT,
    substrate_raw        TEXT,             -- BRENDA substrate string
    co2_form             TEXT CHECK (co2_form IN ('CO2','HCO3-','unspecified')),
    condition_class      TEXT NOT NULL CHECK (condition_class IN ('use','convert','exclude')),
    exclude_reason       TEXT,             -- mutant | ambiguous_substrate | non_physiological | broadcast_duplicate | ...
    conversion_note      TEXT,             -- e.g. 'HCO3- -> CO2 at pH 7.5, pKa 6.35'
    is_mutant            INTEGER DEFAULT 0,
    mutation             TEXT,             -- D117H ...
    ph                   REAL,
    temperature_c        REAL,
    organism_brenda      TEXT,
    measurement_group_id TEXT NOT NULL,    -- one literature measurement; broadcast rows share it
    match_level          TEXT CHECK (match_level IN ('uniprot','organism_ec')),  -- how the BRENDA value reached this sequence
    brenda_ref           TEXT,             -- BRENDA literature id / PubMed
    commentary           TEXT,
    source_id            INTEGER REFERENCES sources(source_id)
);
CREATE INDEX idx_km_seq   ON km_evidence(sequence_id);
CREATE INDEX idx_km_class ON km_evidence(condition_class);
CREATE INDEX idx_km_group ON km_evidence(measurement_group_id);

-- ---------------------------------------------------------------------
-- 3. Annotation: long tables + catalog (tool output, normalised)
-- ---------------------------------------------------------------------
CREATE TABLE feature_catalog (             -- one row per external accession
    accession     TEXT PRIMARY KEY,        -- PF00016 | PTHR42769 | G3DSA:... | PS00088 | MF_01338 ...
    database      TEXT NOT NULL,           -- Pfam | PANTHER | Gene3D | SUPERFAMILY | PROSITE | HAMAP | CDD | TIGRFAM
    name          TEXT,                    -- RuBisCO_large
    description   TEXT,
    interpro_acc  TEXT,
    db_release    TEXT
);

CREATE TABLE domain_hits (                 -- every Pfam / InterPro hit with scores
    id            INTEGER PRIMARY KEY,
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    accession     TEXT NOT NULL REFERENCES feature_catalog(accession),
    tool          TEXT NOT NULL,           -- hmmscan | interproscan
    seq_from      INTEGER,
    seq_to        INTEGER,
    hmm_from      INTEGER,
    hmm_to        INTEGER,
    i_evalue      REAL,
    bitscore      REAL,
    run_id        INTEGER REFERENCES runs(run_id)
);
CREATE INDEX idx_dh_seq ON domain_hits(sequence_id);
CREATE INDEX idx_dh_acc ON domain_hits(accession);

-- ---------------------------------------------------------------------
-- 4. ML feature tables (wide; copied from v2.2 where unchanged)
-- ---------------------------------------------------------------------
CREATE TABLE features_composition (        -- 489 features, script 04a (unchanged)
    sequence_id INTEGER PRIMARY KEY REFERENCES sequences(id),
    aac_A REAL, aac_C REAL, aac_D REAL, aac_E REAL, aac_F REAL, aac_G REAL, aac_H REAL, aac_I REAL, aac_K REAL, aac_L REAL, aac_M REAL, aac_N REAL, aac_P REAL, aac_Q REAL, aac_R REAL, aac_S REAL, aac_T REAL, aac_V REAL, aac_W REAL, aac_Y REAL,
    phys_length REAL, phys_length_log REAL, phys_mw REAL, phys_pi REAL,
    phys_charge_ph7 REAL, phys_gravy REAL, phys_aromaticity REAL,
    phys_instability REAL, phys_frac_glycine REAL, phys_frac_proline REAL,
    phys_frac_charged REAL, phys_frac_aromatic REAL, phys_frac_polar REAL,
    phys_frac_nonpolar REAL, phys_frac_small REAL,
    inv_cat_D REAL, inv_cat_E REAL, inv_cat_H REAL, inv_cat_K REAL,
    inv_cat_C REAL, inv_cat_S REAL, inv_cat_T REAL,
    inv_cat_mean_dist REAL, inv_cat_std_dist REAL,
    inv_cat_min_dist REAL, inv_cat_max_dist REAL, inv_cat_clustering REAL,
    inv_hydrophobic REAL, inv_charged REAL, inv_polar REAL,
    inv_aromatic REAL, inv_net_charge REAL,
    motif_rubisco_kk REAL, motif_rubisco_gk REAL,
    motif_ca_hh REAL, motif_ca_his_cluster REAL,
    motif_pepc_rr REAL, motif_biotin_mk REAL, motif_biotin_amk REAL,
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
    run_id INTEGER REFERENCES runs(run_id)
);

CREATE TABLE features_domains (            -- ML view of domain_hits: 18 carboxylase Pfam flags
    sequence_id INTEGER PRIMARY KEY REFERENCES sequences(id),
    pfam_PF00016 INTEGER DEFAULT 0,   -- RuBisCO large subunit
    pfam_PF02788 INTEGER DEFAULT 0,   -- RuBisCO small subunit
    pfam_PF00101 INTEGER DEFAULT 0,   -- RuBisCO-like
    pfam_PF00194 INTEGER DEFAULT 0,   -- Carbonic anhydrase
    pfam_PF03119 INTEGER DEFAULT 0,   -- Biotin carboxylase N-term
    pfam_PF00311 INTEGER DEFAULT 0,   -- PEPC
    pfam_PF00821 INTEGER DEFAULT 0,   -- Pyruvate carboxylase
    pfam_PF02785 INTEGER DEFAULT 0,   -- Biotin carboxylase C-term
    pfam_PF00364 INTEGER DEFAULT 0,   -- Biotin/lipoyl attachment
    pfam_PF01039 INTEGER DEFAULT 0,   -- Carbamoyl-phosphate synthase
    pfam_PF02786 INTEGER DEFAULT 0,   -- CPSase large N
    pfam_PF02787 INTEGER DEFAULT 0,   -- CPSase large C
    pfam_PF00289 INTEGER DEFAULT 0,   -- CPSase small
    pfam_PF01309 INTEGER DEFAULT 0,   -- PEPC C-terminal
    pfam_PF03599 INTEGER DEFAULT 0,   -- Acetyl-CoA carboxylase
    pfam_PF03590 INTEGER DEFAULT 0,   -- ACC central domain
    pfam_PF00384 INTEGER DEFAULT 0,   -- Molybdopterin oxidoreductase
    pfam_PF00682 INTEGER DEFAULT 0,   -- HMGL-like
    pfam_n_hits INTEGER DEFAULT 0,
    run_id INTEGER REFERENCES runs(run_id)
);

CREATE TABLE features_interpro (           -- unchanged from v2.2 (JSON lists kept for the app)
    sequence_id           INTEGER PRIMARY KEY REFERENCES sequences(id),
    panther_family        TEXT, panther_subfamily TEXT, n_panther INTEGER DEFAULT 0,
    cath_superfamily      TEXT, gene3d_domains_json TEXT, n_gene3d INTEGER DEFAULT 0,
    tigrfam_hits_json     TEXT, n_tigrfam INTEGER DEFAULT 0,
    superfamily_json      TEXT, scop_primary TEXT,
    cdd_hits_json         TEXT, hamap_hits_json TEXT,
    prosite_profiles_json TEXT, n_prosite_prof INTEGER DEFAULT 0,
    prosite_patterns_json TEXT, n_prosite_pat INTEGER DEFAULT 0,
    raw_ipr_json          TEXT,
    run_id INTEGER REFERENCES runs(run_id)
);

CREATE TABLE features_esm2 (               -- 1280 x float32 = 5120 bytes
    sequence_id    INTEGER PRIMARY KEY REFERENCES sequences(id),
    embedding_blob BLOB NOT NULL,
    model_version  TEXT DEFAULT 'esm2_t33_650M_UR50D',
    run_id INTEGER REFERENCES runs(run_id)
);

CREATE TABLE features_ankh (               -- 1536 x float32 = 6144 bytes; stored, not used in v5.2 models
    sequence_id    INTEGER PRIMARY KEY REFERENCES sequences(id),
    embedding_blob BLOB NOT NULL,
    model_version  TEXT DEFAULT 'ankh-large',
    run_id INTEGER REFERENCES runs(run_id)
);

CREATE TABLE features_blast (              -- BLASTp vs training positives (leave-split-out)
    sequence_id       INTEGER PRIMARY KEY REFERENCES sequences(id),
    blast_best_hit    INTEGER REFERENCES sequences(id),
    blast_best_pident REAL DEFAULT 0.0,
    blast_best_evalue REAL DEFAULT 999.0,
    blast_best_bits   REAL,
    blast_best_ec     TEXT DEFAULT '',
    blast_has_hit     INTEGER DEFAULT 0,
    run_id INTEGER REFERENCES runs(run_id)
);

-- Motifs (Aparna handover: motifs_regex.tsv -> motif_catalog; scan -> motif_hits)
CREATE TABLE motif_catalog (
    motif_id        TEXT PRIMARY KEY,      -- e.g. RBCL_Q149, CARBOX_GEN_01
    kind            TEXT NOT NULL CHECK (kind IN ('regex','pwm','hmm','residue')),
    pattern         TEXT,                  -- regex / consensus; residue: 'Q@149'
    reference       TEXT,                  -- residue motifs: reference HMM/sequence (Form I HMM, P00875)
    ec_scope        TEXT,                  -- EC class it applies to, or 'all'
    role            TEXT,                  -- general_carboxylase | low_km_indicator | ...
    discovery       TEXT,                  -- MEME run id / literature
    evidence        TEXT,                  -- controls passed (C-recovery, C-null, ...)
    library_version TEXT NOT NULL,         -- motifs_aparna_vN
    author          TEXT DEFAULT 'Aparna P K'
);

CREATE TABLE motif_hits (                  -- long; only present hits stored
    sequence_id INTEGER NOT NULL REFERENCES sequences(id),
    motif_id    TEXT    NOT NULL REFERENCES motif_catalog(motif_id),
    present     INTEGER NOT NULL DEFAULT 1,
    score       REAL,
    position    INTEGER,
    run_id      INTEGER REFERENCES runs(run_id),
    PRIMARY KEY (sequence_id, motif_id)
);

-- ---------------------------------------------------------------------
-- 5. Models, predictions, confidence
-- ---------------------------------------------------------------------
CREATE TABLE models (
    model_id        TEXT PRIMARY KEY,      -- binary_full_v5.2 | ec_full_v5.2 | km_global_v5.2 | km_ec_4.1.1.39_v5.2 | binary_noesm_v5.2 ...
    model_version   TEXT NOT NULL,         -- v5.2
    target          TEXT NOT NULL CHECK (target IN ('binary','ec','km_global','km_within_ec')),
    variant         TEXT NOT NULL,         -- full | no_esm | pfam_comp | esm_only | motifs_only | rf | mlp
    ec_scope        TEXT,                  -- km_within_ec only
    learner         TEXT,                  -- xgboost | random_forest | mlp
    feature_set     TEXT,                  -- manifest groups used
    n_features      INTEGER,
    feature_order   TEXT,                  -- path to feature_names_<model>.json shipped with the model
    model_file      TEXT,
    sha256          TEXT,
    trained_run_id  INTEGER REFERENCES runs(run_id),
    in_app          INTEGER DEFAULT 0      -- 1 = served by the webapp as a mode
);

CREATE TABLE predictions (                 -- one row per sequence per variant
    sequence_id      INTEGER NOT NULL REFERENCES sequences(id),
    model_version    TEXT NOT NULL,        -- v5.2
    variant          TEXT NOT NULL DEFAULT 'full',
    co2_prob         REAL NOT NULL,
    is_co2_pred      INTEGER NOT NULL,     -- co2_prob >= threshold in db_metadata
    -- EC and Km are filled ONLY when is_co2_pred = 1 (fixes EC for non-carboxylases)
    ec_pred          TEXT,
    ec_prob          REAL,
    ec_top3_json     TEXT,
    km_global_mM     REAL,
    km_within_ec_mM  REAL,                 -- NULL if EC not trainable for within-EC
    km_reported_mM   REAL,                 -- value shown to users
    km_model_used    TEXT CHECK (km_model_used IN ('global','within_ec')),
    km_in_domain     INTEGER,              -- 0 = predicted EC outside Km training classes
    explanation_json TEXT,                 -- grouped local SHAP (app)
    run_id           INTEGER REFERENCES runs(run_id),
    PRIMARY KEY (sequence_id, model_version, variant)
);
CREATE INDEX idx_pred_ec ON predictions(ec_pred);
CREATE INDEX idx_pred_km ON predictions(km_reported_mM);

CREATE TABLE confidence_scores (
    sequence_id      INTEGER PRIMARY KEY REFERENCES sequences(id),
    method_agreement INTEGER,              -- n of {model, BLAST, Pfam, PANTHER} agreeing on EC
    ec_confidence    REAL,
    km_confidence    REAL,
    overall_score    REAL,
    confidence_label TEXT CHECK (confidence_label IN ('high','medium','low','review'))
);

-- ---------------------------------------------------------------------
-- 6. Benchmark (paper tables + TP/TN lists come only from here)
-- ---------------------------------------------------------------------
CREATE TABLE benchmark_results (
    id            INTEGER PRIMARY KEY,
    run_id        INTEGER REFERENCES runs(run_id),
    method        TEXT NOT NULL,           -- carbodb_full | esm_only | pfam_rule | blast_nn | prosite | motifs_only | unikp | catpred | rf | mlp
    target        TEXT NOT NULL,           -- binary | ec | km_global | km_within_ec
    stratum_type  TEXT NOT NULL DEFAULT 'all',  -- all | ec | kingdom | identity_bin | fragment | no_pfam
    stratum_value TEXT NOT NULL DEFAULT 'all',
    metric        TEXT NOT NULL,           -- auroc | auprc | f1 | top1 | f1_macro | r2 | pearson | within_2x
    value         REAL,
    ci_low        REAL,
    ci_high       REAL,
    n             INTEGER,
    n_groups      INTEGER                  -- Km: distinct measurement groups
);
CREATE INDEX idx_bench ON benchmark_results(target, method, stratum_type);

CREATE TABLE benchmark_calls (             -- the TP/TN lists (supplement), test split only
    method        TEXT NOT NULL,
    target        TEXT NOT NULL,           -- binary | ec
    sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    true_value    TEXT NOT NULL,
    pred_value    TEXT,
    score         REAL,
    outcome       TEXT CHECK (outcome IN ('TP','FP','TN','FN','correct','wrong')),
    identity_bin  TEXT,
    PRIMARY KEY (method, target, sequence_id)
);

-- ---------------------------------------------------------------------
-- 7. Views for app and supplement
-- ---------------------------------------------------------------------
CREATE VIEW v_entry AS
SELECT s.cdb_id, s.uniprot_id, s.organism, s.kingdom, s.length, s.label, s.source, s.reviewed,
       s.ec_number AS ec_source, s.km_best_mM AS km_experimental_mM,
       p.co2_prob, p.is_co2_pred, p.ec_pred, p.ec_prob,
       p.km_reported_mM, p.km_model_used, p.km_in_domain,
       c.confidence_label
FROM sequences s
LEFT JOIN predictions p ON p.sequence_id = s.id AND p.variant = 'full'
                       AND p.model_version = (SELECT value FROM db_metadata WHERE key='model_version')
LEFT JOIN confidence_scores c ON c.sequence_id = s.id;

CREATE VIEW v_domain_hits AS               -- what the app shows: name, description, e-value, bit score
SELECT s.cdb_id, s.uniprot_id, h.accession, f.database, f.name, f.description,
       h.seq_from, h.seq_to, h.i_evalue, h.bitscore
FROM domain_hits h
JOIN sequences s ON s.id = h.sequence_id
LEFT JOIN feature_catalog f ON f.accession = h.accession;

CREATE VIEW v_km_training AS               -- exactly what the Km models train on
SELECT k.sequence_id, k.ec_number, k.km_value_mM, k.km_log10_mM, k.measurement_group_id
FROM km_evidence k
JOIN sequences s ON s.id = k.sequence_id
WHERE k.condition_class IN ('use','convert') AND s.label = 1;

CREATE VIEW v_tp_tn AS                     -- supplement export
SELECT b.method, b.target, s.cdb_id, s.uniprot_id, s.organism, b.true_value, b.pred_value,
       b.score, b.outcome, b.identity_bin
FROM benchmark_calls b JOIN sequences s ON s.id = b.sequence_id;
