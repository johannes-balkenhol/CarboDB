#!/usr/bin/env python3
"""
03_merge_all_sources.py
=======================
CarboxyDB — Step 3: Merge all sources into a single master dataset.

What this does:
  1. Load BRENDA positives + negatives
  2. Load UniProt SwissProt + TrEMBL positives + UniProt negatives (if present)
  3. Validate sequences (AA alphabet, length, min unique residues)
  4. Assign labels:
       label=1  true carboxylase    (BRENDA/UniProt CO2 EC, not in CO2_RELATED_EC)
       label=2  ancestral CO2-rel.  (EC in CFG.CO2_RELATED_EC — decarboxylases etc.)
       label=0  true negative       (BRENDA negatives + UniProt negatives)
  5. Fix Km join: BRENDA Km TSV has empty uniprot_id — join on ec_number+organism
     against brenda_positives TSV which has uniprot_id populated
  6. Deduplicate globally on uniprot_id (keep highest evidence tier)
  7. Assign CDB_IDs — permanent internal primary keys (CDB000001 format)
     Existing id_map.tsv is respected: only new sequences get new IDs
  8. Write outputs:
       data/primary/master.tsv        — full merged dataset
       data/primary/master.fasta      — FASTA with >CDB_ID|uniprot_id|ec|label
       data/primary/id_map.tsv        — CDB_ID <-> uniprot_id mapping
       data/primary/km_gold.tsv       — sequences with experimental Km only
       data/interim/merge_report.txt  — summary statistics

Usage:
    python scripts/03_merge_all_sources.py
    python scripts/03_merge_all_sources.py --dry-run   # stats only, no write

Annual update behaviour:
    Re-running this script on new BRENDA/UniProt downloads will:
    - Keep all existing CDB_IDs unchanged (read from id_map.tsv)
    - Assign new CDB_IDs only to sequences not yet in id_map.tsv
    - Overwrite master.tsv with the full updated dataset
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Config import ─────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import CFG, PATHS, TS, make_cdb_id, latest_file, setup_logging

log = setup_logging("03_merge")

# ── Constants ─────────────────────────────────────────────────────────────────
VALID_AA   = CFG.SEQ_VALID_AA
MIN_LEN    = CFG.SEQ_MIN_LEN
MAX_LEN    = CFG.SEQ_MAX_LEN
MIN_UNIQUE = CFG.SEQ_MIN_UNIQUE


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Sequence validation
# ═══════════════════════════════════════════════════════════════════════════════

def clean_seq(seq) -> str:
    return str(seq).upper().strip() if seq else ""

def is_valid_seq(seq: str) -> bool:
    if len(seq) < MIN_LEN or len(seq) > MAX_LEN:
        return False
    chars = set(seq)
    invalid = chars - VALID_AA
    if invalid:
        return False
    return len(chars) >= MIN_UNIQUE

def validate_df(df: pd.DataFrame, name: str) -> pd.DataFrame:
    before = len(df)
    df = df.copy()
    df["sequence"] = df["sequence"].apply(clean_seq)
    mask = df["sequence"].apply(is_valid_seq)
    df = df[mask].reset_index(drop=True)
    removed = before - len(df)
    log.info("  %s: %d → %d after validation (%d removed)", name, before, len(df), removed)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Label assignment
# ═══════════════════════════════════════════════════════════════════════════════

def assign_label(ec_number: str, base_label: int) -> int:
    """
    Assign label=2 if EC is in the ancestral CO2-related set,
    otherwise keep the base_label (1 for positives, 0 for negatives).
    """
    if base_label == 1 and str(ec_number).strip() in CFG.CO2_RELATED_EC:
        return CFG.LABEL_ANCESTRAL
    return base_label


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Load BRENDA data
# ═══════════════════════════════════════════════════════════════════════════════

def load_brenda() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load BRENDA positives, negatives, and Km data.
    Fixes the missing uniprot_id in Km TSV by joining on ec_number+organism.
    """
    log.info("── Loading BRENDA ──")

    pos_file = latest_file(PATHS.RAW_BRENDA, "brenda_positives_*.tsv")
    neg_file = latest_file(PATHS.RAW_BRENDA, "brenda_negatives_*.tsv")
    km_file  = latest_file(PATHS.RAW_BRENDA, "brenda_co2_km_*.tsv")

    log.info("  positives:  %s", pos_file.name)
    log.info("  negatives:  %s", neg_file.name)
    log.info("  Km raw:     %s", km_file.name)

    df_pos = pd.read_csv(pos_file, sep="\t", dtype=str).fillna("")
    df_neg = pd.read_csv(neg_file, sep="\t", dtype=str).fillna("")
    df_km  = pd.read_csv(km_file,  sep="\t", dtype=str).fillna("")

    log.info("  Raw: %d positives, %d negatives, %d Km entries",
             len(df_pos), len(df_neg), len(df_km))

    # ── Fix Km uniprot_id: join on ec_number + organism ──────────────────────
    # BRENDA's getKmValue() doesn't return UniProt IDs directly.
    # brenda_positives has uniprot_id; brenda_co2_km has ec_number+organism.
    # Match them to attach uniprot_id to Km entries.

    if "uniprot_id" not in df_km.columns or df_km["uniprot_id"].str.strip().eq("").all():
        log.info("  Fixing Km uniprot_id via ec_number+organism join...")

        def genus_species(s):
            """Extract first two words (genus + species) for fuzzy organism match."""
            parts = str(s).strip().split()
            return " ".join(parts[:2]).lower()

        df_pos["_gs"] = df_pos["organism"].apply(genus_species)
        df_km["_gs"]  = df_km["organism"].apply(genus_species)

        # Build lookup: (ec_number, genus_species) → uniprot_id
        # genus+species matching recovers ~1092/1241 vs ~857/1241 for exact match
        lookup = (
            df_pos[["ec_number", "_gs", "uniprot_id"]]
            .drop_duplicates(subset=["ec_number", "_gs"])
            .set_index(["ec_number", "_gs"])["uniprot_id"]
        )

        df_km["uniprot_id"] = df_km.apply(
            lambda r: lookup.get((r["ec_number"], r["_gs"]), ""),
            axis=1
        )

        n_matched = (df_km["uniprot_id"] != "").sum()
        n_unmatched = len(df_km) - n_matched
        log.info("  Km join result: %d / %d entries matched (genus+species)",
                 n_matched, len(df_km))
        log.info("  Unmatched %d — organisms absent from UniProt (expected)", n_unmatched)

        df_pos.drop(columns=["_gs"], inplace=True)
        df_km.drop(columns=["_gs"], inplace=True)

    # ── Build curated Km: one best Km per uniprot_id ──────────────────────────
    df_km_valid = df_km[df_km["uniprot_id"] != ""].copy()
    df_km_valid["km_value_mM"] = pd.to_numeric(
        df_km_valid["km_value_mM"], errors="coerce"
    )
    df_km_valid = df_km_valid[
        df_km_valid["km_value_mM"].between(CFG.KM_MIN_VALID, CFG.KM_MAX_VALID)
    ]

    def score_km_row(commentary: str) -> int:
        """
        Score a Km measurement for quality — higher = more physiological.
        Penalise mutants, inhibitors, non-standard conditions.
        Reward wild-type, physiological pH/temp.

        Returns integer score; higher is better.
        """
        c = str(commentary).lower()
        score = 0

        # Strong penalties — these are not native enzyme Km values
        if "mutant" in c:           score -= 100
        if "mutation" in c:         score -= 100
        if "variant" in c:          score -= 50
        if "inhibit" in c:          score -= 40
        if "in the presence of" in c: score -= 30  # effector present
        if "absence of zn" in c:    score -= 50
        if "absence of mg" in c:    score -= 30
        if "recombinant" in c:      score -= 10   # mild penalty only

        # Rewards — physiological wild-type conditions
        if "wild-type" in c or "wildtype" in c:  score += 50
        if "native" in c:           score += 30
        if "physiological" in c:    score += 20

        # pH preference: 7.0-8.0 is physiological
        import re
        ph_match = re.search(r'ph\s*(\d+\.?\d*)', c)
        if ph_match:
            ph = float(ph_match.group(1))
            if 7.0 <= ph <= 8.0:   score += 20
            elif 6.5 <= ph <= 8.5: score += 10
            else:                   score -= 10

        # Temperature preference: 25°C or 37°C
        temp_match = re.search(r'(\d+)\s*[°?]?\s*c\b', c)
        if temp_match:
            temp = int(temp_match.group(1))
            if temp in (25, 37):   score += 15
            elif 20 <= temp <= 40: score += 5

        return score

    df_km_valid["_quality_score"] = df_km_valid["commentary"].apply(score_km_row)

    # For each uniprot_id + ec_number: pick the row with highest quality score
    # If tie: take the median of tied rows
    df_km_valid = df_km_valid.sort_values(
        ["uniprot_id", "ec_number", "_quality_score"],
        ascending=[True, True, False]
    )

    # Get best score per uniprot_id+ec
    best_scores = df_km_valid.groupby(
        ["uniprot_id", "ec_number"]
    )["_quality_score"].max().reset_index()
    best_scores.columns = ["uniprot_id", "ec_number", "_best_score"]

    df_km_best = df_km_valid.merge(best_scores, on=["uniprot_id", "ec_number"])
    df_km_best = df_km_best[
        df_km_best["_quality_score"] == df_km_best["_best_score"]
    ]

    # If multiple rows still tied at best score, take median of those values
    km_curated = (
        df_km_best.groupby(["uniprot_id", "ec_number"])["km_value_mM"]
        .median()
        .reset_index()
        .rename(columns={"km_value_mM": "km_best_mM"})
    )
    km_curated["km_log10_mM"] = np.log10(km_curated["km_best_mM"])

    log.info("  Curated Km: %d unique uniprot_id+ec pairs (quality-filtered)",
             len(km_curated))

    # Assign labels to positives
    df_pos["label"]  = df_pos["ec_number"].apply(
        lambda ec: assign_label(ec, CFG.LABEL_POSITIVE)
    )
    df_pos["source"] = "brenda"
    df_pos["reviewed"] = 0

    # Assign labels to negatives
    df_neg["label"]  = CFG.LABEL_NEGATIVE
    df_neg["source"] = "brenda_neg"
    df_neg["reviewed"] = 0

    return df_pos, df_neg, km_curated


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Load UniProt data
# ═══════════════════════════════════════════════════════════════════════════════

def load_uniprot() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load SwissProt positives, TrEMBL positives, and UniProt negatives.
    Returns empty DataFrames if files are not yet present (graceful skip).
    """
    log.info("── Loading UniProt ──")

    sp_rows, tb_rows, neg_rows = [], [], []

    # SwissProt positives
    sp_files = sorted(PATHS.RAW_SWISSPROT.glob("swissport_positives_*.tsv")) + \
               sorted(PATHS.RAW_SWISSPROT.glob("swissprot_positives_*.tsv"))
    if sp_files:
        df_sp = pd.read_csv(sp_files[-1], sep="\t", dtype=str).fillna("")
        df_sp["source"]   = "swissprot"
        df_sp["reviewed"] = 1
        df_sp["label"] = df_sp["ec_number"].apply(
            lambda ec: assign_label(ec, CFG.LABEL_POSITIVE)
        )
        sp_rows.append(df_sp)
        log.info("  SwissProt: %d rows from %s", len(df_sp), sp_files[-1].name)
    else:
        log.warning("  SwissProt: no files found in %s — skipping", PATHS.RAW_SWISSPROT)

    # TrEMBL positives
    tb_files = sorted(PATHS.RAW_TREMBL.glob("trembl_positives_*.tsv"))
    if tb_files:
        df_tb = pd.read_csv(tb_files[-1], sep="\t", dtype=str).fillna("")
        df_tb["source"]   = "trembl"
        df_tb["reviewed"] = 0
        df_tb["label"] = df_tb["ec_number"].apply(
            lambda ec: assign_label(ec, CFG.LABEL_POSITIVE)
        )
        tb_rows.append(df_tb)
        log.info("  TrEMBL: %d rows from %s", len(df_tb), tb_files[-1].name)
    else:
        log.warning("  TrEMBL: no files found in %s — skipping", PATHS.RAW_TREMBL)

    # UniProt negatives
    neg_files = sorted(PATHS.RAW_NEGATIVES.glob("negatives_*.tsv"))
    if neg_files:
        df_neg_uni = pd.read_csv(neg_files[-1], sep="\t", dtype=str).fillna("")
        df_neg_uni["source"]   = "uniprot_neg"
        df_neg_uni["reviewed"] = 0
        df_neg_uni["label"]    = CFG.LABEL_NEGATIVE
        neg_rows.append(df_neg_uni)
        log.info("  UniProt negatives: %d rows from %s",
                 len(df_neg_uni), neg_files[-1].name)
    else:
        log.warning("  UniProt negatives: no files found — skipping")

    df_sp  = pd.concat(sp_rows,  ignore_index=True) if sp_rows  else pd.DataFrame()
    df_tb  = pd.concat(tb_rows,  ignore_index=True) if tb_rows  else pd.DataFrame()
    df_neg = pd.concat(neg_rows, ignore_index=True) if neg_rows else pd.DataFrame()

    return df_sp, df_tb, df_neg


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Harmonise columns across sources
# ═══════════════════════════════════════════════════════════════════════════════

# Master column set — every source gets these columns (missing → empty string)
MASTER_COLS = [
    "uniprot_id", "ec_number", "enzyme_name", "organism", "taxonomy_id",
    "sequence", "length", "label", "source", "reviewed",
    "protein_name", "gene_name", "annotation_score", "lineage_ids",
]

def harmonise(df: pd.DataFrame, source_tag: str) -> pd.DataFrame:
    """Ensure all master columns exist; normalise types."""
    if df.empty:
        return df
    df = df.copy()

    # Column name aliases across sources
    aliases = {
        "accession":      "uniprot_id",
        "organism_name":  "organism",
        "protein_names":  "protein_name",
        "gene_names":     "gene_name",
        "lineage":        "lineage_ids",
    }
    for old, new in aliases.items():
        if old in df.columns and new not in df.columns:
            df = df.rename(columns={old: new})

    # Add missing columns
    for col in MASTER_COLS:
        if col not in df.columns:
            df[col] = ""

    # Normalise types
    df["length"] = pd.to_numeric(df["length"], errors="coerce").fillna(0).astype(int)
    df["label"]  = pd.to_numeric(df["label"],  errors="coerce").fillna(0).astype(int)
    df["reviewed"] = pd.to_numeric(df["reviewed"], errors="coerce").fillna(0).astype(int)

    # Fill length from sequence if missing
    mask = df["length"] == 0
    df.loc[mask, "length"] = df.loc[mask, "sequence"].str.len()

    if not df.empty:
        log.info("  Harmonised %s: %d rows, labels: %s",
                 source_tag,
                 len(df),
                 dict(df["label"].value_counts().sort_index()))
    return df[MASTER_COLS + [c for c in df.columns if c not in MASTER_COLS]]


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Deduplicate globally
# ═══════════════════════════════════════════════════════════════════════════════

# Evidence tier priority: lower = better
SOURCE_TIER = {
    "brenda":      1,
    "swissprot":   2,
    "trembl":      3,
    "brenda_neg":  4,
    "uniprot_neg": 4,
}

def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep one row per uniprot_id, prioritising by evidence tier.
    Sequences with no uniprot_id are kept as-is (BRENDA-only entries).
    """
    before = len(df)

    # Rows without uniprot_id — keep all (can't dedup without ID)
    df_no_uid  = df[df["uniprot_id"].str.strip() == ""].copy()
    df_has_uid = df[df["uniprot_id"].str.strip() != ""].copy()

    # Assign tier for sorting
    df_has_uid["_tier"] = df_has_uid["source"].map(SOURCE_TIER).fillna(9)

    # Sort: best tier first, reviewed first within tier
    df_has_uid = df_has_uid.sort_values(
        ["uniprot_id", "_tier", "reviewed"],
        ascending=[True, True, False]
    )

    # Keep first (= best tier) per uniprot_id
    df_dedup = df_has_uid.drop_duplicates(subset="uniprot_id", keep="first")
    df_dedup = df_dedup.drop(columns=["_tier"])

    result = pd.concat([df_dedup, df_no_uid], ignore_index=True)
    log.info("  Dedup: %d → %d rows (%d removed)", before, len(result), before - len(result))
    return result.reset_index(drop=True)


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Assign CDB_IDs
# ═══════════════════════════════════════════════════════════════════════════════

def assign_cdb_ids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Assign permanent CDB_IDs.
    - Read existing id_map.tsv if it exists (annual update: preserve old IDs)
    - Assign new sequential IDs only to sequences not yet in the map
    - Write updated id_map.tsv
    """
    log.info("── Assigning CDB_IDs ──")

    # Load existing map if present
    existing_map = {}   # uniprot_id → cdb_id
    max_n = 0
    if PATHS.ID_MAP.exists():
        id_map_df = pd.read_csv(PATHS.ID_MAP, sep="\t", dtype=str)
        existing_map = dict(zip(id_map_df["uniprot_id"], id_map_df["cdb_id"]))
        # Find highest existing number
        for cid in existing_map.values():
            try:
                n = int(cid.replace(CFG.CDB_ID_PREFIX, ""))
                max_n = max(max_n, n)
            except ValueError:
                pass
        log.info("  Loaded %d existing CDB_IDs (max = %s)",
                 len(existing_map), make_cdb_id(max_n))

    cdb_ids = []
    counter = max_n
    new_count = 0

    for _, row in df.iterrows():
        uid = str(row.get("uniprot_id", "")).strip()
        if uid and uid in existing_map:
            cdb_ids.append(existing_map[uid])
        else:
            counter += 1
            cid = make_cdb_id(counter)
            cdb_ids.append(cid)
            if uid:
                existing_map[uid] = cid
            new_count += 1

    df = df.copy()
    df.insert(0, "cdb_id", cdb_ids)
    log.info("  Assigned %d new CDB_IDs (%d total)", new_count, counter)

    # Save updated id_map
    PATHS.PRIMARY.mkdir(parents=True, exist_ok=True)
    id_map_out = pd.DataFrame(
        [(v, k) for k, v in existing_map.items()],
        columns=["cdb_id", "uniprot_id"]
    ).sort_values("cdb_id")
    id_map_out.to_csv(PATHS.ID_MAP, sep="\t", index=False)
    log.info("  id_map.tsv written: %d entries", len(id_map_out))

    return df


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Attach Km values to master
# ═══════════════════════════════════════════════════════════════════════════════

def attach_km(df: pd.DataFrame, km_curated: pd.DataFrame) -> pd.DataFrame:
    """Merge curated Km onto master by uniprot_id + ec_number."""
    if km_curated.empty:
        df["km_best_mM"]   = np.nan
        df["km_log10_mM"]  = np.nan
        return df

    df = df.merge(
        km_curated[["uniprot_id", "ec_number", "km_best_mM", "km_log10_mM"]],
        on=["uniprot_id", "ec_number"],
        how="left"
    )
    n_km = df["km_best_mM"].notna().sum()
    log.info("  Sequences with experimental Km: %d", n_km)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Write outputs
# ═══════════════════════════════════════════════════════════════════════════════

def write_fasta(df: pd.DataFrame, path: Path):
    """Write master FASTA with >CDB_ID|uniprot_id|ec_number|label in header."""
    n = 0
    with open(path, "w") as fh:
        for _, row in df.iterrows():
            seq = str(row.get("sequence", "")).strip()
            if not seq:
                continue
            cdb  = row.get("cdb_id", "")
            uid  = row.get("uniprot_id", "")
            ec   = row.get("ec_number", "")
            lbl  = int(row.get("label", 0))
            km   = row.get("km_best_mM", np.nan)
            km_tag = f"|km={km:.6f}" if pd.notna(km) else ""
            fh.write(f">{cdb}|{uid}|{ec}|label={lbl}{km_tag}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
            n += 1
    log.info("  FASTA written: %d sequences → %s", n, path.name)


def write_outputs(df: pd.DataFrame, dry_run: bool = False):
    PATHS.PRIMARY.mkdir(parents=True, exist_ok=True)
    PATHS.INTERIM.mkdir(parents=True, exist_ok=True)

    if dry_run:
        log.info("  [DRY RUN] — no files written")
        return

    # master.tsv
    df.to_csv(PATHS.MASTER_TSV, sep="\t", index=False)
    log.info("  master.tsv written: %d rows → %s", len(df), PATHS.MASTER_TSV)

    # master.fasta
    write_fasta(df, PATHS.MASTER_FASTA)

    # km_gold.tsv — only sequences with experimental Km
    km_cols = ["cdb_id", "uniprot_id", "ec_number", "organism",
               "km_best_mM", "km_log10_mM", "label", "source"]
    km_cols = [c for c in km_cols if c in df.columns]
    df_km = df[df["km_best_mM"].notna()][km_cols].copy()
    df_km.to_csv(PATHS.PRIMARY / "km_gold.tsv", sep="\t", index=False)
    log.info("  km_gold.tsv written: %d rows", len(df_km))


# ═══════════════════════════════════════════════════════════════════════════════
# 10. Summary report
# ═══════════════════════════════════════════════════════════════════════════════

def print_summary(df: pd.DataFrame):
    log.info("\n%s", "=" * 70)
    log.info("MERGE SUMMARY")
    log.info("%s", "=" * 70)
    log.info("  Total sequences:      %d", len(df))
    log.info("  label=1 (carboxylase): %d", (df["label"] == 1).sum())
    log.info("  label=2 (ancestral):   %d", (df["label"] == 2).sum())
    log.info("  label=0 (negative):    %d", (df["label"] == 0).sum())
    log.info("  With Km:               %d", df["km_best_mM"].notna().sum())

    if "km_best_mM" in df.columns and df["km_best_mM"].notna().any():
        km = df["km_best_mM"].dropna()
        log.info("  Km range:              %.4f – %.2f mM", km.min(), km.max())

    log.info("\n  By source:")
    for src, cnt in df["source"].value_counts().items():
        log.info("    %-16s %d", src, cnt)

    log.info("\n  Top EC classes (label=1):")
    pos = df[df["label"] == 1]
    for ec, cnt in pos["ec_number"].value_counts().head(15).items():
        km_n = pos[pos["ec_number"] == ec]["km_best_mM"].notna().sum()
        log.info("    %-12s %6d seqs  (%d with Km)", ec, cnt, km_n)

    log.info("\n  label=2 EC classes (ancestral CO2-related):")
    anc = df[df["label"] == 2]
    for ec, cnt in anc["ec_number"].value_counts().items():
        log.info("    %-12s %6d seqs", ec, cnt)

    log.info("\n  Next step: python scripts/04_annotate_features.py")


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="Print stats only, do not write files")
    args = ap.parse_args()

    log.info("=" * 70)
    log.info("CarboxyDB — Step 03: Merge All Sources")
    log.info("=" * 70)

    # ── Load ──────────────────────────────────────────────────────────────────
    brenda_pos, brenda_neg, km_curated = load_brenda()
    sp, tb, uni_neg = load_uniprot()

    # ── Harmonise columns ─────────────────────────────────────────────────────
    log.info("── Harmonising columns ──")
    frames = []
    for df, tag in [
        (brenda_pos, "brenda_pos"),
        (brenda_neg, "brenda_neg"),
        (sp,         "swissprot"),
        (tb,         "trembl"),
        (uni_neg,    "uniprot_neg"),
    ]:
        if not df.empty:
            frames.append(harmonise(df, tag))

    df_all = pd.concat(frames, ignore_index=True)
    log.info("  Combined: %d rows before dedup", len(df_all))

    # ── Validate sequences ────────────────────────────────────────────────────
    log.info("── Validating sequences ──")
    df_all = validate_df(df_all, "all sources")

    # ── Remove label=1/2 sequences from negative pool ─────────────────────────
    log.info("── Removing cross-contamination ──")
    pos_uids = set(df_all[df_all["label"] >= 1]["uniprot_id"].str.strip())
    before = len(df_all)
    mask_neg_contaminated = (
        (df_all["label"] == 0) &
        (df_all["uniprot_id"].str.strip().isin(pos_uids))
    )
    df_all = df_all[~mask_neg_contaminated].reset_index(drop=True)
    log.info("  Removed %d negatives that also appear as positives",
             before - len(df_all))

    # ── Deduplicate globally ──────────────────────────────────────────────────
    log.info("── Deduplicating globally ──")
    df_all = deduplicate(df_all)

    # ── Assign CDB_IDs ────────────────────────────────────────────────────────
    df_all = assign_cdb_ids(df_all)

    # ── Attach Km ─────────────────────────────────────────────────────────────
    log.info("── Attaching Km values ──")
    df_all = attach_km(df_all, km_curated)

    # ── Write outputs ─────────────────────────────────────────────────────────
    log.info("── Writing outputs ──")
    write_outputs(df_all, dry_run=args.dry_run)

    # ── Summary ───────────────────────────────────────────────────────────────
    print_summary(df_all)


if __name__ == "__main__":
    main()
