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
    print(f"  Unique sequences:             {len(seq_rows)}")
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
    df_km.to_csv(OUT / f"brenda_co2_km_{TS}.tsv", sep="\t", index=False)
    print(f"  Km entries:      {len(df_km):>8,}  → brenda_co2_km_{TS}.tsv")

    # ── 2. Curate Km: one row per UniProt ID ──────────────────────────────────
    # Keep only rows with a UniProt ID and a numeric Km
    df_km_uid = df_km[
        df_km["uniprot_id"].str.len() > 0
    ].copy()
    df_km_uid["km_value_mM"] = pd.to_numeric(df_km_uid["km_value_mM"], errors="coerce")
    df_km_uid = df_km_uid[df_km_uid["km_value_mM"] > 0]
    df_km_uid = df_km_uid[df_km_uid["km_value_mM"] <= 1000]  # remove outliers

    # Best Km = minimum wild-type value per sequence
    df_km_curated = (
        df_km_uid.sort_values("km_value_mM")
        .drop_duplicates(subset="uniprot_id", keep="first")
        [["uniprot_id", "ec_number", "km_value_mM", "organism"]]
        .rename(columns={"km_value_mM": "km_best_mM"})
    )
    df_km_curated.to_csv(OUT / f"brenda_km_curated_{TS}.tsv", sep="\t", index=False)
    print(f"  Unique seqs+Km:  {len(df_km_curated):>8,}  → brenda_km_curated_{TS}.tsv")

    # ── 3. Positive sequences ─────────────────────────────────────────────────
    df_pos = pd.DataFrame(seq_rows)
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
