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
