#!/usr/bin/env python3
"""
build_motif_bundles_v3_1_per_form_balanced.py
==============================================

Generate THREE within-Form Km contrasts (Form I, II, III) using each
Form's OWN median Km as the split point. This avoids the absolute-threshold
problem (Form II/III have no sequences below 0.01 mM in our predictions
because their entire distribution sits at moderate Km), and gives the
colleague three independent within-Form Fisher's-test runs:

    1. Within Form I:   does the v2 result survive?
    2. Within Form II:  do the same residues differ on this lineage?
    3. Within Form III: same question, archaeal lineage.

Cross-comparison of the three results identifies WHICH of the v2 residues
are universal Km-tuners vs which are Form-specific.

Output bundles (under data/motifs_v3_1_per_form_balanced/):

    form_I_low.fasta              — Form I, bottom 50% of Form I Km
    form_I_high.fasta             — Form I, top 50% of Form I Km
    form_II_low.fasta             — Form II, bottom 50% of Form II Km
    form_II_high.fasta            — Form II, top 50% of Form II Km
    form_III_low.fasta            — Form III, bottom 50% of Form III Km
    form_III_high.fasta           — Form III, top 50% of Form III Km
    form_I_experimental_only.fasta — Form I sequences WITH BRENDA/SwissProt Km
    form_classification_summary.tsv — full per-seq audit
    README.md                      — methodology + analysis protocol

Form classification (precedence order):

    Form I      Hamap MF_01338 or MF_00132 in raw_ipr_json
    Form II     Hamap MF_01339 OR CDD cd08211 in raw_ipr_json
    Form III    Hamap MF_01133 OR CDD cd08213 OR TIGR03326 in raw_ipr_json
    Form-?      No marker matched (often Form I missing Hamap annotation)

PANTHER PTHR42704 is NOT used — it matches Forms I, II and III equally.
CDD cd08212 is NOT used — it is a generic RuBisCO_large entry, not Form-III-specific.
IPR017443 is NOT used — it matches Forms II/III/IV-RLP indiscriminately.

Verified on textbook entries: P00875 spinach=I, P04718 R.rubrum=II,
Q58632 M.jannaschii=III.

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
MAX_PER_GENUS = 5             # diversity cap inside any one cohort
MAX_PER_FORM_CELL = 500       # max size for any low/high cell (caps Form I)

# Seed for reproducibility
RNG_SEED = 42

# Form classification rules — precedence order matters.
#
# IMPORTANT: PANTHER (PTHR42704) does NOT distinguish Form — both Form II
# and Form III sequences carry it. CDD cd08212 is also generic (matches
# all Forms). The reliable markers are Hamap families and a few specific
# CDD entries.
#
# Confirmed via diagnostic on textbook entries (2026-05-07):
#   P00875 spinach    → I (Hamap MF_01338)
#   P04718 R.rubrum   → II (Hamap MF_01339 + CDD cd08211)
#   Q58632 M.jannas.  → III (Hamap MF_01133 + CDD cd08213 + TIGR03326)
#
# Distribution before sampling on full DB:
#   Form I: ~55,000 (Hamap MF_01338)
#   Form II: ~400  (Hamap MF_01339 + CDD cd08211)
#   Form III: ~440 (Hamap MF_01133 + CDD cd08213 + TIGR03326)
#   ?: ~93,500 (likely Form I missing Hamap annotation)
#
# (regex_or_string, form_label, where_to_check)
FORM_RULES = [
    # === Hamap — most reliable single marker ===
    ("MF_01338", "I",   "raw_ipr_json"),   # Form I large chain (cyano/plant)
    ("MF_00132", "I",   "raw_ipr_json"),   # Form I cyanobacterial alt
    ("MF_01339", "II",  "raw_ipr_json"),   # Form II large chain
    ("MF_01133", "III", "raw_ipr_json"),   # Form III archaeal

    # === CDD — adds entries Hamap missed ===
    ("cd08211",  "II",  "raw_ipr_json"),   # Form II specific (NOT cd08212!)
    ("cd08213",  "III", "raw_ipr_json"),   # Form III specific subtype

    # === TIGRFAM — additional archaeal Form III ===
    ("TIGR03326", "III", "raw_ipr_json"),

    # NOTE: deliberately NO PTHR42704 (matches all Forms)
    # NOTE: deliberately NO cd08212 (matches all Forms)
    # NOTE: deliberately NO IPR017443 (matches Forms II/III/IV-RLP all)
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




def genus_stratified_take(rows, target_n, rng):
    """Pick up to target_n rows, with at most MAX_PER_GENUS per genus.
    Round-robin across genera so we get diversity early in the sample.
    """
    if not rows or target_n <= 0:
        return []

    by_g = defaultdict(list)
    for r in rows:
        by_g[r.genus].append(r)
    for g in by_g:
        rng.shuffle(by_g[g])

    if len(rows) <= target_n:
        out = []
        for g, gs in by_g.items():
            out.extend(gs[:MAX_PER_GENUS])
        rng.shuffle(out)
        return out[:target_n]

    out = []
    genera = list(by_g.keys())
    rng.shuffle(genera)
    used = defaultdict(int)
    while len(out) < target_n:
        progressed = False
        for g in genera:
            if used[g] >= MAX_PER_GENUS or used[g] >= len(by_g[g]):
                continue
            out.append(by_g[g][used[g]])
            used[g] += 1
            progressed = True
            if len(out) >= target_n:
                break
        if not progressed:
            break
    return out[:target_n]

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

        Form I       Hamap MF_01338 or MF_00132 in raw_ipr_json
        Form II      Hamap MF_01339 OR CDD cd08211 in raw_ipr_json
        Form III     Hamap MF_01133 OR CDD cd08213 OR TIGR03326 in raw_ipr_json
        Form ?       no marker matched (often Form I with missing Hamap)

        PANTHER (PTHR42704) is NOT used — it matches all Forms.
        CDD cd08212 is NOT used — it is a generic RuBisCO_large entry.
        IPR017443 is NOT used — it matches Forms II/III/IV-RLP.

        Verified on textbook entries: P00875=I, P04718=II, Q58632=III.

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




def write_readme_v31(path, bundles, medians):
    lines = []
    lines.append("Motif bundles v3.1 — within-Form balanced Km contrasts")
    lines.append("=========================================================")
    lines.append("")
    lines.append("Generated 2026-05-08. Three independent within-Form Km contrasts")
    lines.append("(Form I, II, III), each split at that Form's OWN median predicted Km.")
    lines.append("Cross-comparison identifies universal vs Form-specific Km signal.")
    lines.append("")
    lines.append("Per-Form median Km used as the split:")
    for f, m in medians.items():
        lines.append(f"    Form {f:<4} median = {m:.4f} mM")
    lines.append("")
    lines.append("Bundle composition:")
    for name, lst in bundles.items():
        lines.append(f"    {name + '.fasta':<40} n={len(lst)}")
    lines.append("")
    lines.append("Selection criteria")
    lines.append("------------------")
    lines.append("- EC = 4.1.1.39, label=1, length 350-600")
    lines.append("- non-NULL predicted Km")
    lines.append("- 19 stuck XGBoost values excluded")
    lines.append(f"- max {MAX_PER_GENUS} sequences per genus per cohort")
    lines.append(f"- target {MAX_PER_FORM_CELL} per Form-I cell; full half for II/III")
    lines.append(f"- random seed {RNG_SEED}")
    lines.append("")
    lines.append("Form classification (in precedence order)")
    lines.append("-----------------------------------------")
    lines.append("    Form I     Hamap MF_01338 or MF_00132 in raw_ipr_json")
    lines.append("    Form II    Hamap MF_01339 OR CDD cd08211")
    lines.append("    Form III   Hamap MF_01133 OR CDD cd08213 OR TIGR03326")
    lines.append("")
    lines.append("Verified textbook entries: P00875=I, P04718=II, Q58632=III.")
    lines.append("")
    lines.append("Suggested analyses")
    lines.append("------------------")
    lines.append("")
    lines.append("Three independent Fisher's-test runs:")
    lines.append("")
    lines.append("    mafft --localpair --maxiterate 1000 form_I_low.fasta   > I_low.aln")
    lines.append("    mafft --localpair --maxiterate 1000 form_I_high.fasta  > I_high.aln")
    lines.append("    # Run Fisher's exact at each alignment column.")
    lines.append("    # Repeat for Form II and Form III pairs.")
    lines.append("")
    lines.append("Compare results across the three Forms:")
    lines.append("    - Universal Km tuner: significant in all three")
    lines.append("    - Form-I-specific:    significant only in I")
    lines.append("    - Phylogenetic noise: was significant in v2 but not in any single Form here")
    lines.append("")
    lines.append("CAVEATS")
    lines.append("-------")
    lines.append("- Form II/III cohorts are small (hundreds, not thousands).")
    lines.append("  Fisher's-test power is reduced; only large effects will reach")
    lines.append("  Bonferroni-corrected significance.")
    lines.append("- Form III low-Km cohort has very few distinct genera (~7).")
    lines.append("  Watch for genus-confounded results in Form III.")
    lines.append("- The Km split is by PREDICTED Km. The 'within-Form low vs high'")
    lines.append("  contrast assumes the model's Km predictions reflect real")
    lines.append("  biological Km variation. For ~283 Form I sequences with")
    lines.append("  experimental Km, see form_I_experimental_only.fasta.")
    lines.append("")
    path.write_text("\n".join(lines))
    log.info("  wrote %s", path.name)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db",  default="data/primary/carbodb.sqlite")
    ap.add_argument("--out", default="data/motifs_v3_1_per_form_balanced")
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

    # ──────────────────────────────────────────────────────────────
    # Split by Form, then within each Form do MEDIAN-Km split
    # ──────────────────────────────────────────────────────────────
    import statistics

    by_form = {"I": [], "II": [], "III": []}
    for r in rows:
        if r.form_label in by_form:
            if r.predicted_km_mM in STUCK_KM_VALUES_mM:
                continue
            by_form[r.form_label].append(r)

    log.info("Sequences per Form (after stuck-value exclusion):")
    for form_lbl, seqs in by_form.items():
        log.info("    Form %-4s %d", form_lbl, len(seqs))

    bundles = {}
    medians = {}
    for form_lbl, seqs in by_form.items():
        if len(seqs) < 4:
            log.warning("Form %s has only %d sequences — skipping",
                        form_lbl, len(seqs))
            bundles[f"form_{form_lbl}_low"]  = []
            bundles[f"form_{form_lbl}_high"] = []
            continue

        kms = [s.predicted_km_mM for s in seqs]
        med = statistics.median(kms)
        medians[form_lbl] = med
        log.info("    Form %s median Km = %.4f mM (n=%d)",
                 form_lbl, med, len(seqs))

        low_pool  = [s for s in seqs if s.predicted_km_mM <  med]
        high_pool = [s for s in seqs if s.predicted_km_mM >= med]

        target = min(len(low_pool), len(high_pool), MAX_PER_FORM_CELL)

        bundles[f"form_{form_lbl}_low"]  = genus_stratified_take(low_pool,  target, rng)
        bundles[f"form_{form_lbl}_high"] = genus_stratified_take(high_pool, target, rng)

    form_I_with_exp = [r for r in by_form["I"] if r.exp_km_mM is not None]
    bundles["form_I_experimental_only"] = form_I_with_exp

    sampled_all = []
    for v in bundles.values():
        sampled_all.extend(v)

    log.info("\nFinal bundle composition:")
    for name, lst in bundles.items():
        log.info("    %-32s n=%d", name, len(lst))

    for name, lst in bundles.items():
        write_fasta(out_dir / f"{name}.fasta", lst)

    write_summary_tsv(out_dir / "form_classification_summary.tsv", sampled_all)
    write_readme_v31(out_dir / "README.md", bundles, medians)

    log.info("\nDone. Bundle ready in %s", out_dir)
    log.info("To create a tarball for emailing:")
    log.info("    tar czf %s.tar.gz -C %s %s",
             out_dir.name, out_dir.parent, out_dir.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
