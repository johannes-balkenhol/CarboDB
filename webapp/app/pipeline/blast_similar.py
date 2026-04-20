"""
blast_similar.py — nearest-neighbor lookup via pre-built per-EC BLAST DBs.

Flow:
    user sequence + predicted EC
    -> find EC-specific BLAST DB from manifest.json
    -> write query to temp FASTA
    -> run blastp with outfmt 6 (tabular)
    -> parse top K hits
    -> for each hit: lookup experimental Km + curation tier from SQLite
    -> return enriched list with tier badges
"""

import json
import os
import sqlite3
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, List, Dict


# Module-level cache for the manifest — read once, reuse across requests
_MANIFEST_CACHE: Optional[dict] = None


def _load_manifest(manifest_path: Path) -> dict:
    global _MANIFEST_CACHE
    if _MANIFEST_CACHE is None:
        if not manifest_path.exists():
            _MANIFEST_CACHE = {}
        else:
            with open(manifest_path) as f:
                _MANIFEST_CACHE = json.load(f)
    return _MANIFEST_CACHE


def _tier_of(source: str, reviewed: int, has_exp_km: bool) -> Dict[str, str]:
    """Assign a display tier badge for a sequence, based on curation + Km evidence."""
    source = (source or "").lower()
    if has_exp_km:
        if source == "swissprot" or reviewed:
            return {"tier": "swissprot_experimental",
                    "tier_label": "SwissProt · experimental Km"}
        if source == "brenda":
            return {"tier": "brenda_experimental",
                    "tier_label": "BRENDA · experimental Km"}
        return {"tier": "trembl_experimental",
                "tier_label": f"{source.title() or 'DB'} · experimental Km"}
    # No experimental Km — fall back to curation-only labels
    if source == "swissprot" or reviewed:
        return {"tier": "swissprot_predicted",
                "tier_label": "SwissProt · predicted Km"}
    if source == "brenda":
        return {"tier": "brenda_predicted",
                "tier_label": "BRENDA · predicted Km"}
    return {"tier": "trembl_predicted",
            "tier_label": f"{source.title() or 'DB'} · predicted Km"}


def _lookup_hit_metadata(conn, uniprot_id: str) -> dict:
    """Fetch organism, source, reviewed, km_experimental, km_predicted in one pass."""
    # Prefer CO2/HCO3- substrates if multiple Km entries exist
    row = conn.execute("""
        SELECT s.uniprot_id, s.organism, s.source, s.reviewed, s.length,
               (SELECT km_value_mM * 1000
                FROM km_evidence
                WHERE uniprot_id = s.uniprot_id
                  AND evidence_tier = 1
                ORDER BY
                    CASE substrate
                        WHEN 'CO2' THEN 1
                        WHEN 'HCO3-' THEN 2
                        WHEN 'bicarbonate' THEN 3
                        ELSE 4
                    END
                LIMIT 1) AS km_exp_uM,
               (SELECT substrate
                FROM km_evidence
                WHERE uniprot_id = s.uniprot_id
                  AND evidence_tier = 1
                LIMIT 1) AS km_exp_substrate,
               (SELECT km_pred_mM * 1000
                FROM predictions p
                WHERE p.sequence_id = s.id) AS km_pred_uM
        FROM sequences s
        WHERE s.uniprot_id = ?
    """, (uniprot_id,)).fetchone()
    if row is None:
        return {}
    return {
        "uniprot_id": row[0], "organism": row[1], "source": row[2] or "",
        "reviewed": bool(row[3]), "length": row[4],
        "km_experimental_uM": round(row[5], 2) if row[5] is not None else None,
        "km_exp_substrate": row[6],
        "km_predicted_uM": round(row[7], 2) if row[7] is not None else None,
    }


def run_blast_similar(
    sequence: str,
    ec_predicted: str,
    limit: int = 3,
    db_path: Optional[str] = None,
    manifest_path: Optional[str] = None,
) -> List[Dict]:
    """
    Run blastp of `sequence` against the EC-specific clustered DB, return top-K hits
    ranked by sequence identity.

    Returns:
        List of dicts, each with keys: rank, uniprot_id, organism, ec_number,
        identity_pct, align_length, evalue, bitscore, km_experimental_uM,
        km_predicted_uM, source, reviewed, tier, tier_label.
    """
    # Resolve paths
    root = Path(os.environ.get("CARBODB_ROOT", ".")).resolve()
    sqlite_path = Path(db_path or os.environ.get("DB_PATH", "data/primary/carbodb.sqlite"))
    mfest_path = Path(manifest_path or "data/blast_ec_dbs/manifest.json")

    if not sqlite_path.is_absolute():
        sqlite_path = root / sqlite_path
    if not mfest_path.is_absolute():
        mfest_path = root / mfest_path

    manifest = _load_manifest(mfest_path)
    ec_entry = manifest.get(ec_predicted)
    if not ec_entry:
        # No BLAST DB for this EC — can't run. Return empty.
        return []

    blast_db = Path(ec_entry["db_path"])
    if not blast_db.is_absolute():
        blast_db = root / blast_db

    # Write query sequence to a temp FASTA
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fa", delete=False, prefix="cdb_query_"
    ) as f:
        f.write(f">query\n{sequence}\n")
        query_fa = f.name

    try:
        # Run blastp. outfmt 6 gives tab-separated fields:
        # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
        cmd = [
            "blastp",
            "-query", query_fa,
            "-db", str(blast_db),
            "-outfmt", "6 sseqid pident length evalue bitscore",
            "-max_target_seqs", str(limit * 3),   # over-fetch, filter dupes later
            "-evalue", "1e-5",
            "-num_threads", "2",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            # Log stderr but don't raise — neighbor lookup is non-critical
            print(f"BLAST failed for EC {ec_predicted}: {result.stderr[:200]}")
            return []
        raw = result.stdout.strip()
    finally:
        try:
            os.unlink(query_fa)
        except OSError:
            pass

    if not raw:
        return []

    # Parse hits — header is "uid|EC|source|reviewed|organism" from makeblastdb
    hits = []
    seen_uids = set()
    for line in raw.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        sseqid, pident, align_len, evalue, bitscore = parts[:5]
        # sseqid format: {uid}|{EC}|{source}|{reviewed}|{organism}
        uid = sseqid.split("|", 1)[0].strip()
        if not uid or uid in seen_uids:
            continue
        seen_uids.add(uid)
        hits.append({
            "uniprot_id": uid,
            "identity_pct": round(float(pident), 1),
            "align_length": int(align_len),
            "evalue": float(evalue),
            "bitscore": round(float(bitscore), 1),
        })

    # Enrich each hit with DB metadata (organism, Km, tier)
    conn = sqlite3.connect(str(sqlite_path), timeout=10)
    try:
        enriched = []
        for rank, h in enumerate(hits[:limit], 1):
            meta = _lookup_hit_metadata(conn, h["uniprot_id"])
            if not meta:
                continue
            has_km = meta.get("km_experimental_uM") is not None
            tier_info = _tier_of(meta["source"], int(meta["reviewed"]), has_km)
            enriched.append({
                "rank": rank,
                "uniprot_id": h["uniprot_id"],
                "organism": meta.get("organism"),
                "ec_number": ec_predicted,
                "identity_pct": h["identity_pct"],
                "align_length": h["align_length"],
                "evalue": h["evalue"],
                "bitscore": h["bitscore"],
                "km_experimental_uM": meta.get("km_experimental_uM"),
                "km_exp_substrate": meta.get("km_exp_substrate"),
                "km_predicted_uM": meta.get("km_predicted_uM"),
                "source": meta.get("source"),
                "reviewed": meta.get("reviewed"),
                "tier": tier_info["tier"],
                "tier_label": tier_info["tier_label"],
            })
    finally:
        conn.close()

    return enriched


if __name__ == "__main__":
    # Minimal self-test: spinach RuBisCO large subunit against the 4.1.1.39 DB
    import sys
    seq = ("MSPQTETKAGAGFKAGVKDYRLTYYTPDYVVRDTDILAAFRMTPQPGVPPEECGAAVAAESSTGTWTTVWT"
           "DGLTSLDRYKGRCYDIEPVPGEDNQYIAFVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIPP"
           "AYSKTFQGPPHGIQVERDKLNKYGRPLLGCTIKPKLGLSAKNYGRAVYECLRGGLDFTKDDENVNSQPFMR"
           "WRDRFLFCAEAIYKAQAETGEIKGHYLNATAGTCEEMIKRAVFARELGVPIVMHDYLTGGFTANTSLQYC"
           "RDNGLLLHIHRAMHAVIDRQKNHGMHFRVLAKALRLSGGDHIHSGTVVGKLEGERDITLGFVDLLRDDFIE"
           "KDRSRGIYFTQDWVSLPGVIPVASGGIHVWHMPALTEIFGDDSVLQFGGGTLGHPWGNAPGAVANRVALEA"
           "CVQARNEGRDLAREGNEIIREACKWSPELAAACEVWKEIKFEFPAMDTV")
    ec = sys.argv[1] if len(sys.argv) > 1 else "4.1.1.39"
    hits = run_blast_similar(seq, ec, limit=3)
    print(json.dumps(hits, indent=2))
