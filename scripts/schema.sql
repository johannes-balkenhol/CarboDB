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
-- 489 features total. dipeptide_json (400) + pseudo_aac_json (30) as JSON blobs
-- to avoid 430 extra columns (~60MB vs ~900MB flat for 2.4M rows).
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
    dipeptide_json  TEXT,   -- 400 dipeptide frequencies
    pseudo_aac_json TEXT,   -- 30 PseAAC correlation factors
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
