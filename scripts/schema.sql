-- CarboxyDB SQLite Schema v2
-- Evidence-tiered, matches working CarboxyPred_v2 database from January 2026.
-- All Km values stored in mM with explicit unit column.
-- ═══════════════════════════════════════════════════════════════════════════

PRAGMA journal_mode = WAL;
PRAGMA foreign_keys = ON;

-- ── Core sequence table ──────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS sequences (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    uniprot_id       TEXT    UNIQUE NOT NULL,
    entry_name       TEXT,
    sequence         TEXT    NOT NULL,
    length           INTEGER NOT NULL,
    organism         TEXT,
    taxonomy_id      INTEGER,
    lineage          TEXT,
    gene_name        TEXT,
    protein_name     TEXT,
    annotation_score REAL,
    reviewed         INTEGER DEFAULT 0,   -- 1 = SwissProt, 0 = TrEMBL/other
    label            INTEGER NOT NULL,    -- 1 = positive, 0 = negative
    seq_valid        INTEGER DEFAULT 1,   -- 0 = failed AA validation
    source           TEXT,                -- brenda|swissport|trembl|negative
    created_at       TEXT DEFAULT (datetime('now')),
    updated_at       TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_seq_uid    ON sequences(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_seq_label  ON sequences(label);
CREATE INDEX IF NOT EXISTS idx_seq_org    ON sequences(organism);
CREATE INDEX IF NOT EXISTS idx_seq_len    ON sequences(length);
CREATE INDEX IF NOT EXISTS idx_seq_source ON sequences(source);

-- ── EC evidence — multiple sources per UID ───────────────────────────────────
-- Evidence tiers:
--   1 = experimental  (BRENDA entry with measured Km)
--   2 = curated       (SwissProt manually reviewed)
--   3 = predicted     (TrEMBL computational / ML model output)
--   4 = inferred      (BLAST best-hit / Pfam domain)
CREATE TABLE IF NOT EXISTS ec_evidence (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    NOT NULL,
    ec_number       TEXT    NOT NULL,
    source          TEXT    NOT NULL,  -- brenda|swissport|trembl|model_v3|model_v5|pfam|blast
    evidence_type   TEXT    NOT NULL,  -- experimental|curated|predicted|inferred
    evidence_tier   INTEGER NOT NULL,
    confidence      REAL,              -- model probability (for predicted rows)
    model_version   TEXT,              -- v3|v5
    created_at      TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_ec_uid    ON ec_evidence(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_ec_num    ON ec_evidence(ec_number);
CREATE INDEX IF NOT EXISTS idx_ec_tier   ON ec_evidence(evidence_tier);
CREATE INDEX IF NOT EXISTS idx_ec_src    ON ec_evidence(source);

-- ── Km evidence — multiple sources per UID ───────────────────────────────────
-- IMPORTANT: km_value_mM is in mM. The old DB stored µM — we now use mM.
-- Web app multiplies by 1000 to display µM.
CREATE TABLE IF NOT EXISTS km_evidence (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    NOT NULL,
    ec_number       TEXT,
    km_value_mM     REAL,              -- Km in millimolar (mM)
    km_log10_mM     REAL,              -- log10(km_value_mM) for regression
    km_unit         TEXT DEFAULT 'mM', -- always mM; explicit for clarity
    substrate       TEXT,              -- CO2|HCO3-|bicarbonate|etc
    organism        TEXT,
    source          TEXT    NOT NULL,  -- brenda|model_v3|model_v5
    evidence_type   TEXT    NOT NULL,  -- experimental|predicted
    evidence_tier   INTEGER NOT NULL,
    reference       TEXT,
    commentary      TEXT,
    model_version   TEXT,
    created_at      TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_km_uid    ON km_evidence(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_km_ec     ON km_evidence(ec_number);
CREATE INDEX IF NOT EXISTS idx_km_val    ON km_evidence(km_value_mM);
CREATE INDEX IF NOT EXISTS idx_km_tier   ON km_evidence(evidence_tier);

-- ── Best evidence view ───────────────────────────────────────────────────────
-- Returns the highest-priority EC + Km per sequence (lowest tier number = best).
CREATE VIEW IF NOT EXISTS best_evidence AS
SELECT
    s.uniprot_id,
    s.organism,
    s.length,
    s.label,
    s.reviewed,
    s.source,
    ec.ec_number,
    ec.evidence_tier    AS ec_tier,
    ec.source           AS ec_source,
    km.km_value_mM,
    km.km_log10_mM,
    km.km_unit,
    km.evidence_tier    AS km_tier,
    km.source           AS km_source
FROM sequences s
LEFT JOIN ec_evidence ec ON ec.uniprot_id = s.uniprot_id
    AND ec.evidence_tier = (
        SELECT MIN(e2.evidence_tier) FROM ec_evidence e2
        WHERE e2.uniprot_id = s.uniprot_id
    )
LEFT JOIN km_evidence km ON km.uniprot_id = s.uniprot_id
    AND km.evidence_tier = (
        SELECT MIN(k2.evidence_tier) FROM km_evidence k2
        WHERE k2.uniprot_id = s.uniprot_id
    );

-- ── Sequence composition features ────────────────────────────────────────────
-- AA composition, dipeptides (JSON blob), PseAAC (JSON blob), physicochemical
CREATE TABLE IF NOT EXISTS features_composition (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    -- AA composition (20 values)
    aac_A REAL, aac_C REAL, aac_D REAL, aac_E REAL, aac_F REAL,
    aac_G REAL, aac_H REAL, aac_I REAL, aac_K REAL, aac_L REAL,
    aac_M REAL, aac_N REAL, aac_P REAL, aac_Q REAL, aac_R REAL,
    aac_S REAL, aac_T REAL, aac_V REAL, aac_W REAL, aac_Y REAL,
    -- Physicochemical
    phys_length       REAL,
    phys_length_log   REAL,
    phys_mw           REAL,
    phys_pi           REAL,
    phys_charge_ph7   REAL,
    phys_gravy        REAL,
    phys_aromaticity  REAL,
    phys_instability  REAL,
    phys_frac_glycine REAL,
    phys_frac_proline REAL,
    phys_frac_charged REAL,
    phys_frac_aromatic REAL,
    phys_frac_polar   REAL,
    phys_frac_nonpolar REAL,
    phys_frac_small   REAL,
    -- Catalytic core features
    inv_cat_D REAL, inv_cat_E REAL, inv_cat_H REAL, inv_cat_K REAL,
    inv_cat_C REAL, inv_cat_S REAL, inv_cat_T REAL,
    inv_cat_mean_dist REAL, inv_cat_std_dist REAL,
    inv_cat_min_dist  REAL, inv_cat_max_dist REAL,
    inv_cat_clustering REAL,
    inv_hydrophobic REAL, inv_charged REAL, inv_polar REAL,
    inv_aromatic    REAL, inv_net_charge REAL,
    -- EC motifs
    motif_rubisco_kk REAL, motif_rubisco_gk REAL,
    motif_ca_hh      REAL, motif_ca_his_cluster REAL,
    motif_pepc_rr    REAL, motif_biotin_mk REAL, motif_biotin_amk REAL,
    -- Bulk JSON blobs (400 dipeptides + 30 PseAAC stored as JSON)
    dipeptide_json  TEXT,
    pseudo_aac_json TEXT,
    feature_version TEXT DEFAULT 'v3'
);
CREATE INDEX IF NOT EXISTS idx_comp_uid ON features_composition(uniprot_id);

-- ── Pfam domain features ──────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS features_domains (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    pfam_hits_json TEXT,           -- JSON array of all hit Pfam accessions
    pfam_PF00016 INTEGER DEFAULT 0, pfam_PF02788 INTEGER DEFAULT 0,
    pfam_PF00101 INTEGER DEFAULT 0, pfam_PF00194 INTEGER DEFAULT 0,
    pfam_PF03119 INTEGER DEFAULT 0, pfam_PF00311 INTEGER DEFAULT 0,
    pfam_PF00821 INTEGER DEFAULT 0, pfam_PF02785 INTEGER DEFAULT 0,
    pfam_PF00364 INTEGER DEFAULT 0, pfam_PF01039 INTEGER DEFAULT 0,
    pfam_PF02786 INTEGER DEFAULT 0, pfam_PF02787 INTEGER DEFAULT 0,
    pfam_PF00289 INTEGER DEFAULT 0, pfam_PF01309 INTEGER DEFAULT 0,
    pfam_PF03599 INTEGER DEFAULT 0, pfam_PF03590 INTEGER DEFAULT 0,
    pfam_PF00384 INTEGER DEFAULT 0, pfam_PF00682 INTEGER DEFAULT 0,
    pfam_n_hits  INTEGER DEFAULT 0,
    -- PROSITE binary (14 patterns × count + present)
    prosite_PS00157_count   INTEGER DEFAULT 0,
    prosite_PS00157_present INTEGER DEFAULT 0,
    prosite_PS00158_count   INTEGER DEFAULT 0,
    prosite_PS00158_present INTEGER DEFAULT 0,
    prosite_PS00162_count   INTEGER DEFAULT 0,
    prosite_PS00162_present INTEGER DEFAULT 0,
    prosite_PS00188_count   INTEGER DEFAULT 0,
    prosite_PS00188_present INTEGER DEFAULT 0,
    prosite_PS00781_count   INTEGER DEFAULT 0,
    prosite_PS00781_present INTEGER DEFAULT 0,
    prosite_PS00393_count   INTEGER DEFAULT 0,
    prosite_PS00393_present INTEGER DEFAULT 0,
    prosite_PS00017_count   INTEGER DEFAULT 0,
    prosite_PS00017_present INTEGER DEFAULT 0,
    prosite_other_json      TEXT   -- remaining 7 patterns as JSON
);
CREATE INDEX IF NOT EXISTS idx_dom_uid ON features_domains(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_dom_pf16 ON features_domains(pfam_PF00016);
CREATE INDEX IF NOT EXISTS idx_dom_pf194 ON features_domains(pfam_PF00194);

-- ── MEME motif features (pending subproject) ──────────────────────────────────
CREATE TABLE IF NOT EXISTS features_meme (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id  TEXT    UNIQUE NOT NULL,
    meme_hits_json TEXT,  -- {motif_name: 0/1} for all 65 motifs
    meme_n_hits INTEGER DEFAULT 0,
    motif_version TEXT DEFAULT 'pending'  -- filled when MEME subproject completes
);
CREATE INDEX IF NOT EXISTS idx_meme_uid ON features_meme(uniprot_id);

-- ── BLAST homology features ────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS features_blast (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    blast_best_pident REAL DEFAULT 0.0,
    blast_best_evalue REAL DEFAULT 999.0,
    blast_best_ec     TEXT DEFAULT '',
    blast_has_hit     INTEGER DEFAULT 0
);
CREATE INDEX IF NOT EXISTS idx_blast_uid ON features_blast(uniprot_id);

-- ── ESM-2 embeddings (large — GPU job) ────────────────────────────────────────
CREATE TABLE IF NOT EXISTS features_esm2 (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER UNIQUE NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    UNIQUE NOT NULL,
    embedding_blob  BLOB,   -- float32 numpy array, 1280-dim, stored as bytes
    model_version   TEXT DEFAULT 'esm2_t33_650M_UR50D',
    computed_at     TEXT
);
CREATE INDEX IF NOT EXISTS idx_esm_uid ON features_esm2(uniprot_id);

-- ── ML predictions ─────────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS predictions (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id     INTEGER NOT NULL REFERENCES sequences(id),
    uniprot_id      TEXT    NOT NULL,
    model_version   TEXT    NOT NULL,  -- v3|v5
    is_co2_pred     INTEGER,           -- 1 = predicted carboxylase
    co2_prob        REAL,              -- binary classification probability
    ec_pred         TEXT,              -- predicted EC class
    ec_prob         REAL,              -- EC class confidence
    km_pred_mM      REAL,              -- predicted Km in mM
    km_pred_log10   REAL,
    created_at      TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_pred_uid   ON predictions(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_pred_model ON predictions(model_version);
CREATE INDEX IF NOT EXISTS idx_pred_ec    ON predictions(ec_pred);
CREATE INDEX IF NOT EXISTS idx_pred_km    ON predictions(km_pred_mM);

-- ── Confidence scores ──────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS confidence_scores (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence_id      INTEGER UNIQUE REFERENCES sequences(id),
    uniprot_id       TEXT    UNIQUE NOT NULL,
    method_agreement INTEGER,     -- n methods agreeing on EC class (0-4)
    ec_confidence    REAL,        -- mean probability across agreeing methods
    km_confidence    REAL,        -- v3/v5 agreement on Km
    overall_score    REAL,        -- composite score
    confidence_label TEXT         -- high|medium|low|review
);
CREATE INDEX IF NOT EXISTS idx_conf_uid   ON confidence_scores(uniprot_id);
CREATE INDEX IF NOT EXISTS idx_conf_label ON confidence_scores(confidence_label);

-- ── Database metadata ─────────────────────────────────────────────────────────
CREATE TABLE IF NOT EXISTS db_metadata (
    key   TEXT PRIMARY KEY,
    value TEXT
);
INSERT OR REPLACE INTO db_metadata VALUES ('schema_version', '2.0');
INSERT OR REPLACE INTO db_metadata VALUES ('km_unit',        'mM');
INSERT OR REPLACE INTO db_metadata VALUES ('created',        datetime('now'));
INSERT OR REPLACE INTO db_metadata VALUES ('meme_status',    'pending');
