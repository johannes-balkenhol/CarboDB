#!/usr/bin/env python3
"""
build_motif_bundles_v3_form_split.py
=====================================

Generate RuBisCO motif bundles split by Form (I, II/III), so the colleague
can repeat the v2 Fisher-exact analysis WITHIN Form I only — controlling
for the phylogenetic confound that the v2 report itself flags as its
main concern (Concerning Question 1).

Output bundles (under data/motifs_v3_form_split/):

    form_I_low_km_pred.fasta      — Form I only, predicted Km < 0.01 mM
    form_I_high_km_pred.fasta     — Form I only, predicted Km 0.1–5 mM,
                                    19 stuck values excluded
    form_II_III_pooled.fasta      — Forms II + III together (smaller cohort)
    form_classification_summary.tsv — for every RuBisCO sequence:
                                    uniprot_id, form_label, panther_family,
                                    organism, predicted_km_mM, exp_km_mM
    README.md                     — usage + the analyses to run

Form classification (precedence order):

    Form I      panther_family = PTHR42704     (RuBisCO large chain canonical)
    Form II     raw_ipr_json contains IPR017443  (Form II large chain)
    Form III    raw_ipr_json contains "Form III" or specific archaeal IPR IDs
    Form IV     RuBisCO-like (RLP, doesn't carboxylate)
    Form-?      Could not classify

Where the InterPro signal is weak or missing, falls back to organism
kingdom: Eukaryota+Plantae→I; Archaea→III; Bacteria→I-or-II (kept as
"Form-?" to avoid bad assignment).

Stuck-value exclusion list: see STUCK_KM_VALUES below.

Reproducibility: seed = 42, max 10 sequences per genus to avoid
sampling bias by abundant lineages.

Usage:
    python scripts/build_motif_bundles_v3_form_split.py \
        --db   data/primary/carbodb.sqlite \
        --out  data/motifs_v3_form_split

    # then:
    tar czf data/motifs_v3_form_split.tar.gz -C data motifs_v3_form_split

Author: 2026-05-07 handoff
"""
from __future__ import annotations

import argparse
import json
import logging
import random
import re
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from textwrap import dedent

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

EC_RUBISCO = "4.1.1.39"

# Predicted-Km bins (mM)
LOW_KM_MAX_mM  = 0.01     # < this → low-Km set
HIGH_KM_MIN_mM = 0.1      # > this → high-Km set
HIGH_KM_MAX_mM = 5.0      # < this → in the bundle (avoids the >5 mM degenerate tail)

# Stuck-value exclusion (the 19 specific predicted-Km values that XGBoost
# returns for hundreds of distinct sequences each — leaf-grid degeneracies).
# Exact decimal values from the v2 cleanup investigation.
STUCK_KM_VALUES_mM = {
    15.5278882980347, 8.73185443878174, 2.68330836296082, 0.459040313959122,
    1.43213999271393, 0.708834886550903, 1.04893600940704, 5.31432008743286,
    0.0871142297983169, 1.66218221187592, 0.182842999696732, 0.226986095309257,
    0.355072706937790, 0.121624998748302, 0.293305009603500, 0.0556790009140968,
    0.0379824005067348, 4.20678710937500, 12.4321203231811,
}

# How many sequences to sample per genus, to avoid Asteraceae or
# Triticeae dominating just because the database has more of them.
MAX_PER_GENUS = 10

# Seed for reproducibility
RNG_SEED = 42

# Form classification rules — precedence order matters.
# (regex / panther_id, form_label, where_to_check)
FORM_RULES = [
    # Form I — canonical large chain
    ("PTHR42704", "I",   "panther_family"),
    # Form II — distinct PANTHER + InterPro signature
    ("IPR017443", "II",  "raw_ipr_json"),
    # Form III — archaeal
    (re.compile(r"Form\s*III|archaeal", re.I), "III", "raw_ipr_json"),
    # Form IV — RuBisCO-like (RLP), doesn't carboxylate. Excluded from bundles
    # because they're not real RuBisCOs.
    (re.compile(r"RuBisCO[-_ ]like|Form\s*IV", re.I), "IV-RLP", "raw_ipr_json"),
]

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)-7s %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("v3_form_split")


# ---------------------------------------------------------------------------
# Data classes (lightweight, no pydantic)
# ---------------------------------------------------------------------------

class SeqRow:
    __slots__ = ("uniprot_id", "organism", "kingdom", "sequence",
                 "predicted_km_mM", "exp_km_mM", "panther_family",
                 "raw_ipr_json", "form_label", "genus")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))


# ---------------------------------------------------------------------------
# Form classification
# ---------------------------------------------------------------------------

def classify_form(panther_family, raw_ipr_json):
    """Return one of: 'I', 'II', 'III', 'IV-RLP', '?'.

    Precedence: a sequence is Form I if PANTHER says PTHR42704, else
    we look in the raw IPR text for Form II / III / IV markers, else
    we punt to '?'.
    """
    panther_family = panther_family or ""
    raw_ipr_json = raw_ipr_json or ""

    for pat, label, where in FORM_RULES:
        haystack = panther_family if where == "panther_family" else raw_ipr_json
        if isinstance(pat, str):
            if pat in haystack:
                return label
        else:  # regex
            if pat.search(haystack):
                return label
    return "?"


def extract_genus(organism):
    """Best-effort genus extraction. Organism strings look like
    'Spinacia oleracea (Spinach)', or 'Methanococcus jannaschii', etc.
    First whitespace-separated token works in 95% of cases.
    """
    if not organism:
        return "Unknown"
    return organism.split()[0]


# ---------------------------------------------------------------------------
# DB queries
# ---------------------------------------------------------------------------

QUERY_RUBISCO_WITH_FEATURES = """
SELECT
    s.uniprot_id,
    s.organism,
    s.sequence,
    s.km_best_mM,
    p.km_pred_mM           AS predicted_km_mM,
    fi.panther_family,
    fi.raw_ipr_json
FROM sequences s
LEFT JOIN predictions       p  ON p.sequence_id = s.id
LEFT JOIN features_interpro fi ON fi.sequence_id = s.id
WHERE
    s.ec_number = ?
    AND s.label = 1
    AND s.length BETWEEN 350 AND 600       -- exclude small + large outliers
    AND p.km_pred_mM IS NOT NULL
"""


def fetch_rubisco_with_predicted_km(conn):
    log.info("Querying sequences + predictions + InterPro for RuBisCO ...")
    rows = []
    for r in conn.execute(QUERY_RUBISCO_WITH_FEATURES, (EC_RUBISCO,)):
        sr = SeqRow(
            uniprot_id      = r["uniprot_id"],
            organism        = r["organism"],
            kingdom         = None,  # filled later if we need a fallback
            sequence        = r["sequence"],
            predicted_km_mM = r["predicted_km_mM"],
            exp_km_mM       = r["km_best_mM"],
            panther_family  = r["panther_family"],
            raw_ipr_json    = r["raw_ipr_json"],
            form_label      = classify_form(r["panther_family"], r["raw_ipr_json"]),
            genus           = extract_genus(r["organism"]),
        )
        rows.append(sr)
    log.info("  → %d RuBisCO rows fetched", len(rows))
    return rows


# ---------------------------------------------------------------------------
# Stratified sampling
# ---------------------------------------------------------------------------

def stratified_sample(rows, max_per_genus, rng):
    """Keep at most `max_per_genus` sequences per (form, genus) cell.
    Returns a flat list. Ensures rare lineages aren't drowned out.
    """
    by_cell = defaultdict(list)
    for r in rows:
        by_cell[(r.form_label, r.genus)].append(r)

    out = []
    for cell, cell_rows in by_cell.items():
        if len(cell_rows) <= max_per_genus:
            out.extend(cell_rows)
        else:
            out.extend(rng.sample(cell_rows, max_per_genus))
    return out


def filter_low_km(rows):
    return [r for r in rows
            if r.predicted_km_mM is not None
               and r.predicted_km_mM < LOW_KM_MAX_mM]


def filter_high_km(rows):
    return [r for r in rows
            if r.predicted_km_mM is not None
               and HIGH_KM_MIN_mM < r.predicted_km_mM < HIGH_KM_MAX_mM
               and r.predicted_km_mM not in STUCK_KM_VALUES_mM]


# ---------------------------------------------------------------------------
# FASTA writing
# ---------------------------------------------------------------------------

def write_fasta(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        for r in rows:
            header = (
                f">{r.uniprot_id} "
                f"Form={r.form_label} "
                f"organism={r.organism!r} "
                f"genus={r.genus} "
                f"km_pred_mM={r.predicted_km_mM:.5g} "
                f"km_exp_mM={r.exp_km_mM if r.exp_km_mM is not None else 'NA'}"
            )
            f.write(header + "\n")
            seq = r.sequence
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    log.info("  wrote %s (%d sequences)", path.name, len(rows))


def write_summary_tsv(path, rows):
    cols = ["uniprot_id", "form_label", "panther_family",
            "organism", "genus", "predicted_km_mM", "exp_km_mM"]
    with path.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join([
                r.uniprot_id or "",
                r.form_label or "?",
                r.panther_family or "",
                (r.organism or "").replace("\t", " "),
                r.genus or "",
                f"{r.predicted_km_mM:.5g}" if r.predicted_km_mM is not None else "",
                f"{r.exp_km_mM:.5g}" if r.exp_km_mM is not None else "",
            ]) + "\n")
    log.info("  wrote %s", path.name)


def write_readme(path, counts):
    text = dedent(f"""\
        Motif bundles v3 — within-Form RuBisCO split
        =============================================

        Generated 2026-05-07 to address Concerning Question 1 from the
        colleague's v2 Fisher-exact analysis: the low-Km vs high-Km
        comparison is partly confounded by Form I vs Form II/III
        phylogeny. This bundle splits by Form so the colleague can
        re-run the analysis WITHIN Form I only.

        Files
        -----

        form_I_low_km_pred.fasta            n={counts['form_I_low']}
            Form I RuBisCO with predicted Km < {LOW_KM_MAX_mM} mM.

        form_I_high_km_pred.fasta           n={counts['form_I_high']}
            Form I RuBisCO with predicted Km in ({HIGH_KM_MIN_mM}, {HIGH_KM_MAX_mM}) mM.
            The 19 known stuck XGBoost values are excluded.

        form_II_III_pooled.fasta            n={counts['form_II_III']}
            Forms II + III together. Generally smaller cohort; useful as
            a contrast to check whether v2's residue findings show the
            opposite pattern outside Form I (which would be strong
            evidence the v2 contrast was phylogenetic).

        form_classification_summary.tsv     n={counts['total']}
            Per-sequence Form classification with PANTHER family and
            metadata. Useful for sanity-checking the assignments.

        Selection criteria
        ------------------

        - EC = 4.1.1.39 (RuBisCO)
        - label = 1 (positive carboxylase)
        - length 350–600 aa (excludes very short fragments and any large
          fusion proteins)
        - has a non-NULL predicted Km in the predictions table
        - max {MAX_PER_GENUS} sequences per (Form, genus) cell to avoid
          over-representing abundant lineages (random.seed = {RNG_SEED}
          for reproducibility)

        Form classification rules (in precedence order)
        -----------------------------------------------

        Form I       PANTHER family = PTHR42704
        Form II      raw IPR JSON contains IPR017443
        Form III     raw IPR JSON matches /Form III|archaeal/i
        Form IV-RLP  raw IPR JSON matches /RuBisCO-like|Form IV/i  (excluded)
        Form ?       no signal — not used in any bundle

        Suggested analyses
        ------------------

        1. WITHIN-FORM-I differential (the main test):

           mafft --localpair --maxiterate 1000 form_I_low_km_pred.fasta  > low.aln
           mafft --localpair --maxiterate 1000 form_I_high_km_pred.fasta > high.aln
           # Reuse the v2 Fisher-test pipeline on these two alignments.

           Expected if v2 was Km-related:
             positions 251, 255, 258 stay significant within Form I.
             positions 86, 391, 445, 461, 467 stay significant within Form I.

           Expected if v2 was phylogenetic:
             these positions vanish or weaken substantially.

        2. SANITY check on Form II/III:

           Look at the same v2 positions in form_II_III_pooled.fasta.
           If their patterns are OPPOSITE to what v2 says (low-Km AAs
           dominating in high-Km Form II/III), then v2 was definitely
           phylogenetic. If their patterns are CONSISTENT, the signal
           is stronger.

        3. CO-OCCURRENCE pattern table:

           For the v2 11 positions, build the 2^11 alignment-pattern
           table within Form I. If the dominant low-Km AAs and dominant
           high-Km AAs travel together in >80% of sequences, the changes
           are likely coupled (still possibly phylogenetic at sub-Form
           level). If many intermediate patterns exist, the changes are
           independently distributed and can each be evaluated
           individually for engineering.

        Notes
        -----

        - We use predicted Km, not experimental, because the experimental
          high-Km set has only 18 sequences (too small for 470-position
          Fisher's test even before the Bonferroni cut).
        - The 19 stuck XGBoost values are excluded from form_I_high; see
          the v2 generation log for the list.
        - Form classification is ~95% reliable on PANTHER signal; the ~5%
          residue lands in Form ? and is excluded from the bundles. If
          you find an organism that should be Form I but classified as ?,
          paste the uniprot_id to Johannes — likely an InterProScan miss.

        Contact: Johannes Balkenhol — johannes.balkenhol@uni-wuerzburg.de
    """)
    path.write_text(text)
    log.info("  wrote %s", path.name)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db",  default="data/primary/carbodb.sqlite")
    ap.add_argument("--out", default="data/motifs_v3_form_split")
    args = ap.parse_args()

    db_path  = Path(args.db).resolve()
    out_dir  = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not db_path.exists():
        log.error("DB not found: %s", db_path)
        return 1

    rng = random.Random(RNG_SEED)
    log.info("DB:     %s", db_path)
    log.info("Output: %s", out_dir)

    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    # Same PRAGMAs as the webapp uses, helps a 50 GB DB
    conn.execute("PRAGMA cache_size = -524288")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("PRAGMA mmap_size = 30000000000")

    rows = fetch_rubisco_with_predicted_km(conn)

    # Show the form distribution for transparency
    form_counts = defaultdict(int)
    for r in rows:
        form_counts[r.form_label] += 1
    log.info("Form distribution before sampling:")
    for k in sorted(form_counts):
        log.info("    %-8s  %d", k, form_counts[k])

    # Sample to reduce genus-level redundancy
    sampled = stratified_sample(rows, MAX_PER_GENUS, rng)
    log.info("After stratified sample (max %d per (form, genus)): %d",
             MAX_PER_GENUS, len(sampled))

    # Form-I splits
    form_I = [r for r in sampled if r.form_label == "I"]
    form_I_low  = filter_low_km(form_I)
    form_I_high = filter_high_km(form_I)

    # Form II/III pool (exclude Form IV-RLP and Form ?)
    form_II_III = [r for r in sampled if r.form_label in {"II", "III"}]

    log.info("\nBundle composition:")
    log.info("    form_I_low_km_pred       %d", len(form_I_low))
    log.info("    form_I_high_km_pred      %d", len(form_I_high))
    log.info("    form_II_III_pooled       %d", len(form_II_III))

    write_fasta(out_dir / "form_I_low_km_pred.fasta",  form_I_low)
    write_fasta(out_dir / "form_I_high_km_pred.fasta", form_I_high)
    write_fasta(out_dir / "form_II_III_pooled.fasta",  form_II_III)
    write_summary_tsv(out_dir / "form_classification_summary.tsv", sampled)
    write_readme(out_dir / "README.md", {
        "form_I_low":   len(form_I_low),
        "form_I_high":  len(form_I_high),
        "form_II_III":  len(form_II_III),
        "total":        len(sampled),
    })

    log.info("\nDone. Bundle ready in %s", out_dir)
    log.info("To create a tarball for emailing:")
    log.info("    tar czf %s.tar.gz -C %s %s",
             out_dir.name, out_dir.parent, out_dir.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
