#!/usr/bin/env python3
"""
01_brenda_download.py
=====================
CarboxyDB — Step 1: Download everything from BRENDA.

Strategy (confirmed from working January runs):
  1. Call getEcNumbersFromKmValue() → get ALL ~7000 EC numbers
     that have ANY Km data in BRENDA
  2. For each EC, call getKmValue() and check substrate string
     for true CO2/HCO3- terms → builds our positive EC list
  3. For each positive EC, call getSequence() → collect sequences
  4. Call getSequence() on 75 confirmed-negative EC classes
     → negative dataset

Real results from working run (December 2025):
  - 39 EC classes found with true CO2/HCO3- Km values
  - 14,501 Km entries (gold standard)
  - 3,001 unique sequences with Km (training Km subset)
  - 695,384 positive sequences total (all 39 EC classes)
  - 3,324,901 negative sequences

Usage:
    export BRENDA_EMAIL="you@uni.de"
    export BRENDA_PASSWORD="yourpassword"
    python 01_brenda_download.py

    # Or interactive:
    python 01_brenda_download.py --interactive

    # Or inline:
    python 01_brenda_download.py --email you@uni.de --password xxx
"""

import argparse
import getpass
import hashlib
import json
import os
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from zeep import Client, Settings
from zeep.helpers import serialize_object

# ── Output directory ─────────────────────────────────────────────────────────
OUT = Path("data/raw/brenda")
OUT.mkdir(parents=True, exist_ok=True)

TS = datetime.now().strftime("%Y%m%d_%H%M%S")

# ── TRUE CO2/HCO3- substrate terms (lowercase, exact or prefix match) ────────
# Confirmed working — avoids Co2+ (cobalt ion) false positives
TRUE_CO2_TERMS = [
    "co2",
    "carbon dioxide",
    "bicarbonate",
    "hco3-",
    "hco3",
    "hydrogen carbonate",
    "hydrogencarbonate",
]

# ── Confirmed-negative EC classes for negative dataset ───────────────────────
# These were confirmed NOT to have CO2/HCO3- Km in any working BRENDA scan.
# 75 classes covering all 6 EC main divisions.
NEGATIVE_EC = {
    # EC 1 — Oxidoreductases (no CO2)
    "1.1.1.1":   "Alcohol dehydrogenase",
    "1.1.1.27":  "Lactate dehydrogenase",
    "1.1.1.37":  "Malate dehydrogenase",
    "1.1.1.42":  "Isocitrate dehydrogenase",
    "1.1.1.44":  "Phosphogluconate dehydrogenase",
    "1.1.1.49":  "Glucose-6-phosphate 1-dehydrogenase",
    "1.2.1.12":  "Glyceraldehyde-3-phosphate dehydrogenase",
    "1.4.1.2":   "Glutamate dehydrogenase (NAD+)",
    "1.4.1.3":   "Glutamate dehydrogenase (NAD(P)+)",
    "1.6.5.3":   "NADH:ubiquinone oxidoreductase",
    "1.8.1.4":   "Dihydrolipoyl dehydrogenase",
    "1.9.3.1":   "Cytochrome-c oxidase",
    "1.11.1.6":  "Catalase",
    "1.11.1.7":  "Peroxidase",
    "1.14.11.2": "Procollagen-proline 4-dioxygenase",
    "1.14.13.39":"Nitric-oxide synthase",
    "1.15.1.1":  "Superoxide dismutase",
    # EC 2 — Transferases (no CO2)
    "2.1.1.37":  "DNA (cytosine-5-)-methyltransferase",
    "2.1.2.1":   "Glycine hydroxymethyltransferase",
    "2.2.1.1":   "Transketolase",
    "2.3.1.9":   "Acetyl-CoA acetyltransferase",
    "2.3.2.2":   "Gamma-glutamyltransferase",
    "2.4.1.1":   "Glycogen phosphorylase",
    "2.4.1.11":  "Glycogen synthase",
    "2.4.2.1":   "Purine nucleoside phosphorylase",
    "2.5.1.18":  "Glutathione S-transferase",
    "2.6.1.1":   "Aspartate transaminase",
    "2.6.1.2":   "Alanine transaminase",
    "2.7.1.1":   "Hexokinase",
    "2.7.1.11":  "6-Phosphofructokinase",
    "2.7.1.40":  "Pyruvate kinase",
    "2.7.1.69":  "Protein kinase",
    "2.7.4.3":   "Adenylate kinase",
    "2.7.7.7":   "DNA-directed DNA polymerase",
    "2.7.7.48":  "RNA-directed RNA polymerase",
    "2.8.1.1":   "Thiosulfate sulfurtransferase",
    # EC 3 — Hydrolases (no CO2)
    "3.1.1.3":   "Triacylglycerol lipase",
    "3.1.3.1":   "Alkaline phosphatase",
    "3.1.3.2":   "Acid phosphatase",
    "3.1.3.16":  "Protein-serine phosphatase",
    "3.2.1.1":   "Alpha-amylase",
    "3.2.1.17":  "Lysozyme",
    "3.2.1.21":  "Beta-glucosidase",
    "3.2.1.23":  "Beta-galactosidase",
    "3.4.11.1":  "Leucyl aminopeptidase",
    "3.4.21.4":  "Trypsin",
    "3.4.21.1":  "Chymotrypsin",
    "3.4.23.1":  "Pepsin A",
    "3.4.24.3":  "Microbial collagenase",
    "3.5.1.1":   "Asparaginase",
    "3.5.1.2":   "Glutaminase",
    "3.5.4.4":   "Adenosine deaminase",
    "3.6.1.3":   "Apyrase",
    "3.6.4.12":  "DNA helicase",
    # EC 4 — Lyases (ONLY those confirmed not CO2-related)
    "4.1.2.13":  "Fructose-bisphosphate aldolase",
    "4.2.1.2":   "Fumarate hydratase",
    "4.2.1.3":   "Aconitate hydratase",
    "4.2.1.11":  "Phosphopyruvate hydratase (enolase)",
    "4.3.1.3":   "Histidine ammonia-lyase",
    # EC 5 — Isomerases (no CO2)
    "5.1.3.1":   "Ribulose-phosphate 3-epimerase",
    "5.3.1.1":   "Triose-phosphate isomerase",
    "5.3.1.9":   "Glucose-6-phosphate isomerase",
    "5.4.2.1":   "Phosphoglycerate mutase",
    "5.4.2.2":   "Phosphoglucomutase",
    "5.4.99.2":  "Methylmalonyl-CoA mutase",
    # EC 6 — Ligases (only those confirmed not CO2-related)
    "6.1.1.1":   "Tyrosine--tRNA ligase",
    "6.1.1.2":   "Tryptophan--tRNA ligase",
    "6.1.1.5":   "Isoleucine--tRNA ligase",
    "6.2.1.1":   "Acetate--CoA ligase",
    "6.3.2.1":   "Glutamate--cysteine ligase",
    "6.3.2.2":   "Glutathione synthase",
    "6.5.1.1":   "DNA ligase (ATP)",
    "6.5.1.2":   "DNA ligase (NAD+)",
}


# ─────────────────────────────────────────────────────────────────────────────
# BRENDA client
# ─────────────────────────────────────────────────────────────────────────────

class BRENDAClient:

    WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"

    def __init__(self, email: str, password: str):
        self.email = email
        self.pw = hashlib.sha256(password.encode("utf-8")).hexdigest()
        print("Connecting to BRENDA SOAP API...", flush=True)
        settings = Settings(strict=False, xml_huge_tree=True)
        self.client = Client(self.WSDL, settings=settings)
        print("Connected.\n")

    # ── core call — exactly the param format that worked ──────────────────────
    def _call(self, method: str, params: list):
        """Call a BRENDA SOAP method with a list of param strings."""
        try:
            result = getattr(self.client.service, method)(*params)
            return [serialize_object(r) for r in result] if result else []
        except Exception as e:
            print(f"  SOAP error [{method}]: {e}")
            return []

    def get_all_ec_with_km(self) -> list:
        """Get every EC number in BRENDA that has any Km data."""
        print("Getting all EC numbers with Km data from BRENDA...")
        params = [self.email, self.pw]
        result = self.client.service.getEcNumbersFromKmValue(*params)
        ec_list = list(result) if result else []
        print(f"  → {len(ec_list)} EC numbers have Km data in BRENDA")
        return ec_list

    def get_km_values(self, ec: str) -> list:
        """All Km entries for one EC number (all substrates)."""
        params = [
            self.email, self.pw,
            f"ecNumber*{ec}", "organism*", "kmValue*", "kmValueMaximum*",
            "substrate*", "commentary*", "ligandStructureId*", "literature*",
        ]
        return self._call("getKmValue", params)

    def get_sequences(self, ec: str) -> list:
        """All sequences for one EC number."""
        params = [
            self.email, self.pw,
            f"ecNumber*{ec}", "organism*", "sequence*",
            "noOfAminoAcids*", "firstAccessionCode*", "source*", "id*",
        ]
        return self._call("getSequence", params)

    def get_enzyme_name(self, ec: str) -> str:
        """Recommended name for an EC number."""
        try:
            params = [self.email, self.pw, f"ecNumber*{ec}"]
            result = self.client.service.getEnzymeName(*params)
            return serialize_object(result[0]).get("enzymeName", "") if result else ""
        except Exception:
            return ""


# ─────────────────────────────────────────────────────────────────────────────
# Substrate filter
# ─────────────────────────────────────────────────────────────────────────────

def is_true_co2(substrate: str) -> bool:
    """
    Return True only if substrate is CO2 or HCO3-, not Co2+ (cobalt).
    Uses exact or prefix match on lowercased string.
    """
    if not substrate:
        return False
    s = str(substrate).lower().strip()
    # Quick cobalt / false-positive exclusion
    if any(fp in s for fp in ["co2+", "co(2+)", "cobalt", "cobalamin"]):
        return False
    return any(
        s == term or s.startswith(term + " ") or s.startswith(term + "/")
        for term in TRUE_CO2_TERMS
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

def scan_for_co2_enzymes(client: BRENDAClient):
    """
    Scan ALL BRENDA EC numbers for true CO2/HCO3- Km entries.
    Returns: (km_rows, seq_rows, ec_found)
    """
    print("=" * 70)
    print("PHASE 1 — SCANNING ALL EC NUMBERS FOR CO2/HCO3- Km")
    print("=" * 70)

    all_ec = client.get_all_ec_with_km()

    km_rows  = []   # one row per Km measurement
    seq_rows = []   # one row per unique sequence
    ec_found = []   # EC classes that had ≥1 CO2/HCO3- entry
    seen_seqs = set()

    for ec in tqdm(all_ec, desc="Scanning EC classes"):
        time.sleep(0.10)
        km_list = client.get_km_values(ec)

        co2_kms = [e for e in km_list if is_true_co2(e.get("substrate", ""))]
        if not co2_kms:
            continue

        # Record this EC as positive
        ec_found.append(ec)
        name = client.get_enzyme_name(ec)
        print(f"  ✓ EC {ec} ({name}): {len(co2_kms)} CO2 Km entries")

        # Collect Km rows
        for entry in co2_kms:
            km_rows.append({
                "ec_number":       ec,
                "enzyme_name":     name,
                "organism":        str(entry.get("organism", "")),
                "uniprot_id":      str(entry.get("uniprotId", "") or "").strip(),
                "substrate":       str(entry.get("substrate", "")),
                "km_value_mM":     entry.get("kmValue"),
                "km_max_mM":       entry.get("kmValueMaximum"),
                "commentary":      str(entry.get("commentary", "")),
                "literature":      str(entry.get("literature", "")),
            })

        # Collect sequences (download for this EC)
        time.sleep(0.20)
        seq_list = client.get_sequences(ec)
        for s in seq_list:
            uid = str(s.get("firstAccessionCode", "") or "").strip()
            seq = str(s.get("sequence", "") or "").strip()
            if uid and seq and uid not in seen_seqs:
                seen_seqs.add(uid)
                seq_rows.append({
                    "uniprot_id": uid,
                    "ec_number":  ec,
                    "enzyme_name": name,
                    "organism":   str(s.get("organism", "")),
                    "length":     int(s.get("noOfAminoAcids") or len(seq)),
                    "sequence":   seq,
                    "source":     str(s.get("source", "")),
                    "label":      1,
                })
        time.sleep(0.20)

    print(f"\n  EC classes with CO2/HCO3- Km: {len(ec_found)}")
    print(f"  Total Km entries:             {len(km_rows)}")
    print(f"  Unique sequences from EC scan: {len(seq_rows)}")

    # ── Third pass: fetch sequences for Km UniProt IDs not yet in seq_rows ──
    # BRENDA's getSequence() sometimes misses sequences that appear in Km entries.
    # This ensures every sequence with a Km value is downloaded, making the
    # pipeline fully reproducible — no manual patches needed on future runs.
    #
    # Strategy: for each Km entry that has a UniProt ID not yet in seen_seqs,
    # fetch the sequence directly from UniProt REST API.

    import requests as _requests

    # First do the genus+species join to get UniProt IDs for Km entries
    def _gs(s):
        parts = str(s).strip().split()
        return " ".join(parts[:2]).lower()

    seq_lookup = {}
    for row in seq_rows:
        key = (row["ec_number"], _gs(row["organism"]))
        if key not in seq_lookup:
            seq_lookup[key] = row["uniprot_id"]

    # Find Km entries whose UniProt ID is not yet in seq_rows
    km_uids_needed = set()
    for row in km_rows:
        uid = str(row.get("uniprot_id", "") or "").strip()
        if not uid:
            # Try genus+species lookup
            uid = seq_lookup.get((row["ec_number"], _gs(row["organism"])), "")
        if uid and uid not in seen_seqs:
            km_uids_needed.add((uid, row["ec_number"], row["organism"]))

    if km_uids_needed:
        print(f"\n  Fetching {len(km_uids_needed)} sequences for Km entries "
              f"not in BRENDA sequence download...")

        UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        fetched = 0
        failed  = 0

        for uid, ec, organism in tqdm(km_uids_needed,
                                      desc="Fetching Km sequences from UniProt"):
            if uid in seen_seqs:
                continue
            try:
                r = _requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{uid}.fasta",
                    timeout=30
                )
                if r.status_code != 200:
                    failed += 1
                    continue
                lines = r.text.strip().splitlines()
                if not lines or not lines[0].startswith(">"):
                    failed += 1
                    continue
                seq = "".join(lines[1:]).strip()
                if not seq:
                    failed += 1
                    continue

                seen_seqs.add(uid)
                seq_rows.append({
                    "uniprot_id":  uid,
                    "ec_number":   ec,
                    "enzyme_name": "",
                    "organism":    organism,
                    "length":      len(seq),
                    "sequence":    seq,
                    "source":      "uniprot_rest",
                    "label":       1,
                })
                fetched += 1
                time.sleep(0.2)

            except Exception as e:
                failed += 1
                time.sleep(0.5)

        print(f"  Fetched {fetched} additional sequences "
              f"({failed} not found in UniProt)")

    print(f"  Total unique sequences (final): {len(seq_rows)}")
    return km_rows, seq_rows, ec_found


def download_negatives(client: BRENDAClient) -> list:
    """
    Download sequences for confirmed-negative EC classes.
    Returns list of sequence dicts with label=0.
    """
    print("\n" + "=" * 70)
    print(f"PHASE 2 — DOWNLOADING NEGATIVE DATASET ({len(NEGATIVE_EC)} EC classes)")
    print("=" * 70)

    neg_rows = []
    seen = set()
    successful = 0

    for ec, name in tqdm(NEGATIVE_EC.items(), desc="Negative EC classes"):
        time.sleep(0.30)
        seq_list = client.get_sequences(ec)
        n = 0
        for s in seq_list:
            uid = str(s.get("firstAccessionCode", "") or "").strip()
            seq = str(s.get("sequence", "") or "").strip()
            if uid and seq and uid not in seen:
                seen.add(uid)
                neg_rows.append({
                    "uniprot_id": uid,
                    "ec_number":  ec,
                    "enzyme_name": name,
                    "organism":   str(s.get("organism", "")),
                    "length":     int(s.get("noOfAminoAcids") or len(seq)),
                    "sequence":   seq,
                    "source":     str(s.get("source", "")),
                    "label":      0,
                })
                n += 1
        if n:
            successful += 1

    print(f"\n  Successful EC downloads: {successful}/{len(NEGATIVE_EC)}")
    print(f"  Total negative sequences: {len(neg_rows)}")
    return neg_rows


def save_results(km_rows, seq_rows, neg_rows, ec_found):
    """Save all outputs."""
    print("\n" + "=" * 70)
    print("SAVING RESULTS")
    print("=" * 70)

    # ── 1. Km table (gold standard) ───────────────────────────────────────────
    df_km = pd.DataFrame(km_rows)
    print(f"  Km entries:      {len(df_km):>8,}  → brenda_co2_km_{TS}.tsv")

    # ── Build df_pos early — needed for Km→UniProt ID join ───────────────────
    df_pos = pd.DataFrame(seq_rows)

    # ── 2. Curate Km: one row per UniProt ID ──────────────────────────────────
    # PRIMARY KEY IS uniprot_id — Km belongs to a sequence, not a species.
    # Organism is secondary metadata only.
    #
    # Strategy:
    #   a) Use uniprot_id directly from getKmValue() if present (rare)
    #   b) Cross-reference with sequences via genus+species match (main path)
    #   c) Keep unmatched rows with empty uniprot_id for reference only
    #
    # This ensures brenda_km_curated_*.tsv always has uniprot_id populated
    # and annual reruns produce consistent results without manual patching.

    import re as _re

    def _genus_species(s):
        parts = str(s).strip().split()
        return " ".join(parts[:2]).lower()

    # Build lookup from sequences: (ec_number, genus_species) -> uniprot_id
    # If multiple UniProt IDs share same ec+org, keep all (comma-separated)
    # We prefer the first accession code from getSequence()
    seq_lookup = {}
    for _, row in df_pos.iterrows():
        key = (row["ec_number"], _genus_species(row["organism"]))
        if key not in seq_lookup:
            seq_lookup[key] = row["uniprot_id"]

    # Fix Km rows: fill uniprot_id from sequence lookup where missing
    df_km = df_km.copy()
    df_km["_gs"] = df_km["organism"].apply(_genus_species)

    def _get_uid(row):
        # Use direct uniprot_id if already present
        if str(row.get("uniprot_id", "")).strip():
            return str(row["uniprot_id"]).strip()
        # Fall back to sequence lookup
        return seq_lookup.get((row["ec_number"], row["_gs"]), "")

    df_km["uniprot_id"] = df_km.apply(_get_uid, axis=1)
    df_km = df_km.drop(columns=["_gs"])

    # Save the full Km table with uniprot_ids now populated
    df_km.to_csv(OUT / f"brenda_co2_km_{TS}.tsv", sep="\t", index=False)

    n_uid_filled = (df_km["uniprot_id"].str.len() > 0).sum()
    print(f"  Km with UniProt ID: {n_uid_filled:>6,} / {len(df_km)} "
          f"({100*n_uid_filled/max(len(df_km),1):.1f}%)")

    # Curate: one best Km per uniprot_id
    # Best = quality-scored: penalise mutants/inhibitors, reward wild-type/physiological pH
    def _score_km(commentary: str) -> int:
        c = str(commentary).lower()
        score = 0
        if "mutant" in c or "mutation" in c:   score -= 100
        if "variant" in c:                      score -= 50
        if "inhibit" in c:                      score -= 40
        if "in the presence of" in c:           score -= 30
        if "absence of zn" in c:                score -= 50
        if "wild-type" in c or "wildtype" in c: score += 50
        if "native" in c:                       score += 30
        ph = _re.search(r'ph\s*(\d+\.?\d*)', c)
        if ph:
            v = float(ph.group(1))
            score += 20 if 7.0 <= v <= 8.0 else (10 if 6.5 <= v <= 8.5 else -10)
        tmp = _re.search(r'(\d+)\s*[°?]?\s*c\b', c)
        if tmp:
            t = int(tmp.group(1))
            score += 15 if t in (25, 37) else (5 if 20 <= t <= 40 else 0)
        return score

    df_km_uid = df_km[df_km["uniprot_id"].str.len() > 0].copy()
    df_km_uid["km_value_mM"] = pd.to_numeric(df_km_uid["km_value_mM"], errors="coerce")
    df_km_uid = df_km_uid[df_km_uid["km_value_mM"].between(1e-5, 1000)]
    df_km_uid["_score"] = df_km_uid["commentary"].apply(_score_km)

    # Best row per uniprot_id = highest quality score, tie-break by median
    df_km_uid = df_km_uid.sort_values(
        ["uniprot_id", "_score"], ascending=[True, False]
    )
    best_scores = df_km_uid.groupby("uniprot_id")["_score"].max().reset_index()
    df_km_uid = df_km_uid.merge(best_scores.rename(columns={"_score":"_best"}),
                                on="uniprot_id")
    df_km_uid = df_km_uid[df_km_uid["_score"] == df_km_uid["_best"]]

    df_km_curated = (
        df_km_uid.groupby(["uniprot_id", "ec_number"])["km_value_mM"]
        .median()
        .reset_index()
        .rename(columns={"km_value_mM": "km_best_mM"})
    )
    # Attach organism from positives for reference
    org_map = df_pos.set_index("uniprot_id")["organism"].to_dict()
    df_km_curated["organism"] = df_km_curated["uniprot_id"].map(org_map).fillna("")
    df_km_curated["km_log10_mM"] = df_km_curated["km_best_mM"].apply(
        lambda x: __import__("math").log10(x) if x > 0 else float("nan")
    )

    df_km_curated.to_csv(OUT / f"brenda_km_curated_{TS}.tsv", sep="\t", index=False)
    print(f"  Unique seqs+Km:  {len(df_km_curated):>8,}  → brenda_km_curated_{TS}.tsv")

    # ── 3. Positive sequences ─────────────────────────────────────────────────
    # df_pos already built above for Km join
    df_pos.to_csv(OUT / f"brenda_positives_{TS}.tsv", sep="\t", index=False)

    pos_fasta = OUT / f"positives_{TS}.fasta"
    _write_fasta(df_pos, pos_fasta, label=1)
    print(f"  Positive seqs:   {len(df_pos):>8,}  → positives_{TS}.fasta")

    # ── 4. Negative sequences ─────────────────────────────────────────────────
    df_neg = pd.DataFrame(neg_rows)
    df_neg.to_csv(OUT / f"brenda_negatives_{TS}.tsv", sep="\t", index=False)

    neg_fasta = OUT / f"negatives_{TS}.fasta"
    _write_fasta(df_neg, neg_fasta, label=0)
    print(f"  Negative seqs:   {len(df_neg):>8,}  → negatives_{TS}.fasta")

    # ── 5. Combined FASTA ────────────────────────────────────────────────────
    comb_fasta = OUT / f"all_sequences_{TS}.fasta"
    with open(comb_fasta, "wb") as out:
        for f in [pos_fasta, neg_fasta]:
            out.write(f.read_bytes())
    print(f"  Combined:        {len(df_pos)+len(df_neg):>8,}  → all_sequences_{TS}.fasta")

    # ── 6. EC list ────────────────────────────────────────────────────────────
    ec_path = OUT / f"co2_ec_classes_{TS}.txt"
    with open(ec_path, "w") as f:
        for ec in sorted(ec_found):
            f.write(ec + "\n")
    print(f"  CO2 EC classes:  {len(ec_found):>8,}  → co2_ec_classes_{TS}.txt")

    # ── 7. Summary JSON ───────────────────────────────────────────────────────
    summary = {
        "timestamp":           TS,
        "km_entries":          len(df_km),
        "km_curated":          len(df_km_curated),
        "km_range_mM":         [float(df_km_curated["km_best_mM"].min()),
                                 float(df_km_curated["km_best_mM"].max())],
        "ec_classes_found":    sorted(ec_found),
        "n_ec_classes":        len(ec_found),
        "positive_sequences":  len(df_pos),
        "negative_sequences":  len(df_neg),
        "total_sequences":     len(df_pos) + len(df_neg),
        "ratio_neg_pos":       round(len(df_neg) / max(len(df_pos), 1), 2),
    }
    with open(OUT / f"summary_{TS}.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── 8. Print final summary ────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("DOWNLOAD COMPLETE")
    print("=" * 70)
    print(f"  EC classes with CO2 Km:  {len(ec_found)}")
    print(f"  Km entries (raw):        {len(df_km):,}")
    print(f"  Km curated (unique UID): {len(df_km_curated):,}")
    print(f"  Km range:                {df_km_curated['km_best_mM'].min():.4f} – "
          f"{df_km_curated['km_best_mM'].max():.2f} mM")
    print(f"  Positive sequences:      {len(df_pos):,}")
    print(f"  Negative sequences:      {len(df_neg):,}")
    print(f"  Neg/Pos ratio:           {summary['ratio_neg_pos']:.1f}×")
    print(f"\n  EC distribution (top 15):")
    for ec, cnt in df_km["ec_number"].value_counts().head(15).items():
        print(f"    {ec}: {cnt}")
    print(f"\n  All files → {OUT}/")


def _write_fasta(df: pd.DataFrame, path: Path, label: int):
    with open(path, "w") as f:
        for _, row in df.iterrows():
            if row["sequence"]:
                f.write(
                    f">{row['uniprot_id']}|{row['ec_number']}|"
                    f"{row['organism'].replace(' ', '_')}|label={label}\n"
                )
                seq = str(row["sequence"])
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def get_credentials(args) -> tuple[str, str]:
    if args.interactive:
        email = input("BRENDA email: ").strip()
        password = getpass.getpass("BRENDA password: ")
        return email, password
    email    = args.email    or os.environ.get("BRENDA_EMAIL")
    password = args.password or os.environ.get("BRENDA_PASSWORD")
    if not email:
        email = input("BRENDA email: ").strip()
    if not password:
        password = getpass.getpass("BRENDA password: ")
    return email, password


def main():
    parser = argparse.ArgumentParser(
        description="Download BRENDA CO2/HCO3- Km data and sequences"
    )
    parser.add_argument("--email",       default=None)
    parser.add_argument("--password",    default=None)
    parser.add_argument("--interactive", action="store_true",
                        help="Prompt for credentials interactively")
    args = parser.parse_args()

    email, password = get_credentials(args)

    client = BRENDAClient(email, password)

    # Phase 1: scan all EC numbers for CO2 Km → positive dataset
    km_rows, seq_rows, ec_found = scan_for_co2_enzymes(client)

    # Phase 2: download negative sequences
    neg_rows = download_negatives(client)

    # Save everything
    save_results(km_rows, seq_rows, neg_rows, ec_found)


if __name__ == "__main__":
    main()
