#!/usr/bin/env python3
"""
05_migrate_v22_to_v23.py -- CarboDB v5.2, Stage A ("schema DB" for the webapp)

Builds a NEW database with schema v2.3 from the live v2.2 database, in one pass:
  * identity, sequences, EC/Km evidence, all feature tables, v5 predictions are copied
  * Pfam domain hits with e-value / bit score / name come from the rescan parser output
  * EC and Km predictions are kept only where the sequence is predicted a carboxylase
  * Km cleaning is NOT applied yet (db_metadata km_cleaning = pending); every Km row is
    marked condition_class='use' as a placeholder and gets its own measurement_group_id
The old DB is opened read-only. Build on the local SSD, then copy to /storage.

Usage:
  python scripts/05_migrate_v22_to_v23.py \
     --old    /storage/users/job37yv/Projects/CarboDB_v3/data/primary/carbodb.sqlite \
     --schema database/schema_v2.3.sql \
     --out    /var/tmp/job37yv/carbodb_v2.3_dev.sqlite \
     --pfam-dir /storage/users/job37yv/Projects/CarboDB_data/primary/annotations/pfam   (optional)
"""
import argparse, csv, gzip, json, os, socket, sqlite3, subprocess, sys, time

def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)

def cols(con, schema, table):
    return [r[1] for r in con.execute(f"PRAGMA {schema}.table_info({table})")]

def copy_features(con, table, run_id):
    new = [c for c in cols(con, "main", table) if c != "run_id"]
    old = set(cols(con, "old", table))
    missing = [c for c in new if c not in old]
    if missing:
        log(f"  {table}: columns not in old DB, left NULL/default: {missing[:8]}{'...' if len(missing)>8 else ''}")
    use = [c for c in new if c in old]
    sel = ", ".join(f"o.{c}" for c in use)
    con.execute(f"INSERT OR IGNORE INTO main.{table} ({', '.join(use)}, run_id) "
                f"SELECT {sel}, ? FROM old.{table} o JOIN main.sequences s ON s.id = o.sequence_id", (run_id,))
    n = con.execute(f"SELECT COUNT(*) FROM main.{table}").fetchone()[0]
    log(f"  {table}: {n:,} rows")
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--schema", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--pfam-dir")
    ap.add_argument("--force", action="store_true", help="overwrite --out if it exists")
    a = ap.parse_args()
    if os.path.exists(a.out):
        if not a.force: sys.exit(f"{a.out} exists (use --force)")
        os.remove(a.out)
    t0 = time.time()
    con = sqlite3.connect(a.out, isolation_level=None)
    con.executescript(open(a.schema).read())
    con.close()
    # re-open with URI support so the old DB is attached strictly read-only
    con = sqlite3.connect(f"file:{a.out}", uri=True, isolation_level=None)
    con.executescript("PRAGMA foreign_keys=OFF; PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; "
                      "PRAGMA cache_size=-8000000; PRAGMA temp_store=MEMORY;")
    con.execute(f"ATTACH DATABASE 'file:{a.old}?mode=ro' AS old")
    con.execute("BEGIN")

    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        commit = ""
    oldmeta = dict(con.execute("SELECT key, value FROM old.db_metadata").fetchall())
    con.execute("INSERT INTO runs(stage, tool, tool_version, params_json, code_commit, host, started_at, note) "
                "VALUES ('05_migrate','05_migrate_v22_to_v23.py','v2.3',?,?,?,datetime('now'),?)",
                (json.dumps({"old": a.old, "pfam_dir": a.pfam_dir}), commit, socket.gethostname(),
                 "Stage A: v5 content in schema v2.3"))
    run_id = con.execute("SELECT last_insert_rowid()").fetchone()[0]

    # sources (releases to be completed by hand where TODO)
    for name, rel in [("brenda", "2026-03-27"), ("uniprot_swissprot", "2026-03 (TODO exact release)"),
                      ("uniprot_trembl", "2026-03 (TODO exact release)"), ("pfam", "TODO (see raw/dbs/pfam/SOURCE.md)"),
                      ("interpro", "InterProScan 5.72-103.0"), ("esm2", "esm2_t33_650M_UR50D"), ("ankh", "ankh-large")]:
        con.execute("INSERT INTO sources(name, release, licence) VALUES (?,?,?)",
                    (name, rel, "CC BY 4.0" if name.startswith(("brenda", "uniprot")) else None))

    log("identity + sequences")
    con.execute("INSERT INTO main.id_map(id, cdb_id, uniprot_id) SELECT id, cdb_id, uniprot_id FROM old.id_map")
    con.execute("""INSERT INTO main.sequences(id, cdb_id, uniprot_id, ec_number, label, source, reviewed, sequence,
                     length, organism, seq_valid, km_best_mM, km_log10_mM)
                   SELECT id, cdb_id, uniprot_id, ec_number, label, source, reviewed, sequence,
                     length, organism, seq_valid, km_best_mM, km_log10_mM FROM old.sequences""")
    nseq = con.execute("SELECT COUNT(*) FROM main.sequences").fetchone()[0]; log(f"  sequences: {nseq:,}")

    log("evidence")
    con.execute("""INSERT INTO main.ec_evidence(sequence_id, ec_number, source, evidence_tier)
                   SELECT e.sequence_id, e.ec_number, e.source, e.evidence_tier FROM old.ec_evidence e
                   WHERE e.evidence_tier BETWEEN 1 AND 4""")
    con.execute("""INSERT INTO main.km_evidence(sequence_id, ec_number, km_value_mM, km_log10_mM, km_raw_value,
                     km_raw_unit, substrate_raw, co2_form, condition_class, exclude_reason, measurement_group_id, commentary)
                   SELECT k.sequence_id, COALESCE(k.ec_number, s.ec_number), k.km_value_mM, k.km_log10_mM,
                     k.km_value_mM, k.km_unit, k.substrate,
                     CASE WHEN lower(k.substrate) LIKE '%hco3%' OR lower(k.substrate) LIKE '%bicarbonate%'
                               OR lower(COALESCE(k.commentary,'')) LIKE '%hco3%' THEN 'HCO3-'
                          WHEN upper(k.substrate) = 'CO2' THEN 'CO2' ELSE 'unspecified' END,
                     CASE WHEN k.source = 'brenda' THEN 'use' ELSE 'exclude' END,
                     CASE WHEN k.source = 'brenda' THEN NULL
                          ELSE 'origin_unverified (v5 source=' || k.source || '); only BRENDA Km is used' END,
                     'v5_km_' || k.id, k.commentary
                   FROM old.km_evidence k JOIN old.sequences s ON s.id = k.sequence_id""")
    log(f"  ec_evidence: {con.execute('SELECT COUNT(*) FROM ec_evidence').fetchone()[0]:,}  "
        f"km_evidence: {con.execute('SELECT COUNT(*) FROM km_evidence').fetchone()[0]:,}")

    log("feature tables (copied unchanged)")
    nf = {t: copy_features(con, t, run_id) for t in
          ["features_composition", "features_domains", "features_interpro", "features_esm2", "features_ankh"]}

    nhits = 0
    if a.pfam_dir and os.path.exists(os.path.join(a.pfam_dir, "domain_hits_pfam.tsv.gz")):
        log("Pfam catalog + domain hits")
        with open(os.path.join(a.pfam_dir, "feature_catalog_pfam.tsv")) as f:
            r = csv.DictReader(f, delimiter="\t")
            con.executemany("INSERT OR IGNORE INTO feature_catalog(accession, database, name, description, db_release) "
                            "VALUES (?,?,?,?,?)",
                            ((x["accession"], "Pfam", x["name"], x["description"], x["acc_version"]) for x in r))
        idmap = dict(con.execute("SELECT cdb_id, id FROM sequences"))
        def rows():
            with gzip.open(os.path.join(a.pfam_dir, "domain_hits_pfam.tsv.gz"), "rt") as f:
                for x in csv.DictReader(f, delimiter="\t"):
                    sid = idmap.get(x["cdb_id"])
                    if sid is None: continue
                    yield (sid, x["accession"], "hmmscan", x["seq_from"], x["seq_to"], x["hmm_from"], x["hmm_to"],
                           float(x["i_evalue"]), float(x["bitscore"]), run_id)
        con.executemany("INSERT INTO domain_hits(sequence_id, accession, tool, seq_from, seq_to, hmm_from, hmm_to, "
                        "i_evalue, bitscore, run_id) VALUES (?,?,?,?,?,?,?,?,?,?)", rows())
        nhits = con.execute("SELECT COUNT(*) FROM domain_hits").fetchone()[0]
        log(f"  feature_catalog: {con.execute('SELECT COUNT(*) FROM feature_catalog').fetchone()[0]:,}  domain_hits: {nhits:,}")
    else:
        log("no Pfam rescan output given -> domain_hits empty (app falls back to features_domains.pfam_hits_json)")

    log("predictions (v5; EC/Km only for predicted carboxylases)")
    con.execute("""INSERT OR IGNORE INTO main.predictions(sequence_id, model_version, variant, co2_prob, is_co2_pred,
                     ec_pred, ec_prob, km_global_mM, km_reported_mM, km_model_used, km_in_domain, run_id)
                   SELECT p.sequence_id, p.model_version, 'full', COALESCE(p.co2_prob, 0), COALESCE(p.is_co2_pred, 0),
                     CASE WHEN p.is_co2_pred = 1 THEN p.ec_pred END,
                     CASE WHEN p.is_co2_pred = 1 THEN p.ec_prob END,
                     CASE WHEN p.is_co2_pred = 1 THEN p.km_pred_mM END,
                     CASE WHEN p.is_co2_pred = 1 THEN p.km_pred_mM END,
                     CASE WHEN p.is_co2_pred = 1 AND p.km_pred_mM IS NOT NULL THEN 'global' END,
                     CASE WHEN p.is_co2_pred = 1 THEN
                          (p.ec_pred IN (SELECT DISTINCT ec_number FROM old.km_evidence WHERE source='brenda')) END,
                     ?
                   FROM old.predictions p JOIN main.sequences s ON s.id = p.sequence_id""", (run_id,))
    con.execute("""INSERT OR IGNORE INTO main.confidence_scores(sequence_id, method_agreement, ec_confidence,
                     km_confidence, overall_score, confidence_label)
                   SELECT c.sequence_id, c.method_agreement, c.ec_confidence, c.km_confidence, c.overall_score,
                     CASE WHEN c.confidence_label IN ('high','medium','low','review') THEN c.confidence_label END
                   FROM old.confidence_scores c JOIN main.sequences s ON s.id = c.sequence_id""")

    mdir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(a.old))), "data", "models")
    for mid, tgt, f in [("binary_full_v5", "binary", "binary_v5.json"), ("ec_full_v5", "ec", "ec_v5.json"),
                        ("km_global_v5", "km_global", "km_v5_final.json")]:
        con.execute("INSERT INTO models(model_id, model_version, target, variant, learner, model_file, in_app) "
                    "VALUES (?,?,?,?,?,?,1)", (mid, "v5", tgt, "full", "xgboost", os.path.join(mdir, f)))

    stats = {r[0]: r[1] for r in con.execute("SELECT 'n_label'||label, COUNT(*) FROM sequences GROUP BY label")}
    meta = {"schema_version": "2.3", "release": "v5.2-dev-stageA", "model_version": "v5",
            "content": "v5 predictions and features in schema v2.3; NOT the final v5.2 release",
            "km_cleaning": "pending", "n_sequences": nseq, "built_at": time.strftime("%Y-%m-%d %H:%M"),
            "code_commit": commit, "source_db": a.old, "source_schema_version_claimed": oldmeta.get("schema_version"),
            "n_domain_hits": nhits, **stats}
    con.executemany("INSERT OR REPLACE INTO db_metadata VALUES (?,?)", [(k, str(v)) for k, v in meta.items()])
    con.execute("COMMIT")

    log("integrity checks")
    chk = {
        "sequences without composition": "SELECT COUNT(*) FROM sequences s LEFT JOIN features_composition f ON f.sequence_id=s.id WHERE f.sequence_id IS NULL",
        "sequences without esm2": "SELECT COUNT(*) FROM sequences s LEFT JOIN features_esm2 f ON f.sequence_id=s.id WHERE f.sequence_id IS NULL",
        "sequences without features_domains": "SELECT COUNT(*) FROM sequences s LEFT JOIN features_domains f ON f.sequence_id=s.id WHERE f.sequence_id IS NULL",
        "sequences without predictions": "SELECT COUNT(*) FROM sequences s LEFT JOIN predictions p ON p.sequence_id=s.id WHERE p.sequence_id IS NULL",
        "predicted carboxylases": "SELECT COUNT(*) FROM predictions WHERE is_co2_pred=1",
        "EC/Km removed from non-carboxylases": "SELECT COUNT(*) FROM old.predictions WHERE COALESCE(is_co2_pred,0)=0 AND (ec_pred IS NOT NULL OR km_pred_mM IS NOT NULL)",
        "predicted carboxylases with Pfam scores": "SELECT COUNT(DISTINCT h.sequence_id) FROM domain_hits h JOIN predictions p ON p.sequence_id=h.sequence_id WHERE p.is_co2_pred=1",
        "foreign-key violations": "SELECT COUNT(*) FROM pragma_foreign_key_check",
    }
    res = {k: con.execute(q).fetchone()[0] for k, q in chk.items()}
    con.execute("DETACH DATABASE old"); con.execute("ANALYZE")
    con.close()
    rep = [f"# Stage A build report ({time.strftime('%Y-%m-%d %H:%M')}, {(time.time()-t0)/60:.1f} min)", "",
           f"- output: {a.out} ({os.path.getsize(a.out)/1e9:.1f} GB)", f"- sequences: {nseq:,}"]
    rep += [f"- {k}: {v:,}" for k, v in nf.items()] + [f"- domain_hits: {nhits:,}", "", "## Checks"]
    rep += [f"- {k}: {v:,}" for k, v in res.items()]
    txt = "\n".join(rep) + "\n"
    open(os.path.splitext(a.out)[0] + "_build_report.md", "w").write(txt)
    print(txt)

if __name__ == "__main__":
    main()
