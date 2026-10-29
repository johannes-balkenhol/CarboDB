#!/usr/bin/env python3
"""04b_parse_domtblout.py -- CarboDB v5.2, fix F2 (Pfam e-values, bit scores, names).
Reads existing hmmscan --domtblout files (no rerun). Read-only on inputs and DB; writes only into --out:
domain_hits_pfam.tsv.gz, feature_catalog_pfam.tsv, missing_features_domains.tsv, parse_report.md"""
import argparse, glob, gzip, json, os, random, sqlite3, sys, time
from collections import Counter
COLS = ["target_name","target_acc","tlen","query_name","query_acc","qlen","full_evalue","full_score",
        "full_bias","dom_n","dom_of","c_evalue","i_evalue","dom_score","dom_bias","hmm_from","hmm_to",
        "ali_from","ali_to","env_from","env_to","acc"]
def parse_line(line):
    f = line.rstrip("\n").split(None, 22)
    if len(f) < 22: return None
    r = dict(zip(COLS, f[:22])); r["description"] = f[22] if len(f) > 22 else ""; return r
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domtbl-dir", required=True); ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--check-sample", type=int, default=20000)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.domtbl_dir, "*.domtblout")))
    if not files: sys.exit(f"no *.domtblout in {a.domtbl_dir}")
    t0 = time.time(); catalog = {}; hits = {}; best_ie = {}; n_hits = bad = 0; per_file = []
    with gzip.open(os.path.join(a.out, "domain_hits_pfam.tsv.gz"), "wt") as fo:
        fo.write("cdb_id\tuniprot_id\taccession\tacc_version\tseq_from\tseq_to\tenv_from\tenv_to\thmm_from\thmm_to\t"
                 "i_evalue\tc_evalue\tbitscore\tfull_evalue\tfull_score\tdom_n\tdom_of\tsource_file\n")
        for fp in files:
            n_file = 0
            for line in open(fp):
                if line.startswith("#") or not line.strip(): continue
                r = parse_line(line)
                if r is None: bad += 1; continue
                q = r["query_name"].split("|"); cdb, uid = q[0], (q[1] if len(q) > 1 else "")
                acc_v = r["target_acc"]; acc = acc_v.split(".")[0]
                if acc not in catalog:
                    catalog[acc] = (acc_v.split(".")[1] if "." in acc_v else "", r["target_name"], r["description"], r["tlen"])
                fo.write("\t".join([cdb, uid, acc, acc_v, r["ali_from"], r["ali_to"], r["env_from"], r["env_to"],
                    r["hmm_from"], r["hmm_to"], r["i_evalue"], r["c_evalue"], r["dom_score"], r["full_evalue"],
                    r["full_score"], r["dom_n"], r["dom_of"], os.path.basename(fp)]) + "\n")
                hits.setdefault(cdb, set()).add(acc); k = (cdb, acc); ie = float(r["i_evalue"])
                if ie < best_ie.get(k, 1e9): best_ie[k] = ie
                n_hits += 1; n_file += 1
            per_file.append(n_file)
    with open(os.path.join(a.out, "feature_catalog_pfam.tsv"), "w") as fc:
        fc.write("accession\tdatabase\tacc_version\tname\tdescription\tmodel_length\n")
        for acc in sorted(catalog):
            v, name, desc, L = catalog[acc]; fc.write(f"{acc}\tPfam\t{v}\t{name}\t{desc}\t{L}\n")
    con = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    all_cdb = {r[0]: r[1] for r in con.execute("SELECT cdb_id, id FROM sequences")}
    in_fd = {r[0] for r in con.execute("SELECT s.cdb_id FROM features_domains d JOIN sequences s ON s.id=d.sequence_id")}
    missing = set(all_cdb) - in_fd; missing_hits = missing & set(hits); unknown = set(hits) - set(all_cdb)
    random.seed(1); sample = random.sample(sorted(in_fd), min(a.check_sample, len(in_fd))); db_sets = {}
    for i in range(0, len(sample), 900):
        ch = sample[i:i+900]
        for cdb, js in con.execute("SELECT s.cdb_id, d.pfam_hits_json FROM features_domains d JOIN sequences s "
                                   "ON s.id=d.sequence_id WHERE s.cdb_id IN (%s)" % ",".join("?"*len(ch)), ch):
            db_sets[cdb] = set(json.loads(js)) if js else set()
    agree = {}
    for cut in (None, 1e-3, 1e-5, 1e-10):
        ok = sum(({x for x in hits.get(c, ()) if cut is None or best_ie[(c, x)] < cut} == s) for c, s in db_sets.items())
        agree["none" if cut is None else f"{cut:g}"] = ok / max(1, len(db_sets))
    con.close()
    with open(os.path.join(a.out, "missing_features_domains.tsv"), "w") as fm:
        fm.write("cdb_id\tsequence_id\thas_parsed_hits\n")
        for c in sorted(missing): fm.write(f"{c}\t{all_cdb[c]}\t{int(c in hits)}\n")
    top = Counter(x for s in hits.values() for x in s).most_common(10)
    rep = (f"# 04b parse report ({time.strftime('%Y-%m-%d %H:%M')}, {time.time()-t0:.0f} s)\n\n"
           f"- domtblout files: {len(files)} (hits/file min {min(per_file)}, max {max(per_file)})\n"
           f"- domain hits parsed: {n_hits:,}   malformed lines: {bad}\n"
           f"- sequences with >=1 Pfam hit: {len(hits):,} of {len(all_cdb):,}\n"
           f"- distinct Pfam families: {len(catalog):,}\n- CDB ids in domtblout but not in DB: {len(unknown)}\n\n"
           f"## Missing in features_domains: {len(missing):,}\n"
           f"- with parsed hits (insert bug, fixable): {len(missing_hits):,}\n"
           f"- without hits (no Pfam hit or never scanned): {len(missing)-len(missing_hits):,}\n\n"
           f"## DB pfam_hits_json vs domtblout (n={len(db_sets):,})\n| i-evalue cutoff | exact agreement |\n|---|---|\n"
           + "".join(f"| {k} | {v:.3%} |\n" for k, v in agree.items())
           + "\n## Most frequent families\n" + "".join(f"- {x} {catalog[x][1]}: {n:,}\n" for x, n in top))
    open(os.path.join(a.out, "parse_report.md"), "w").write(rep); print(rep)
if __name__ == "__main__": main()
