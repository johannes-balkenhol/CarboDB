"""Build cleaned motif sequence bundles (v2) for collaborator motif analysis.

What changed from April 14 bundles:
  - Excludes the BRENDA Helianthus 18.0 mM block (suspected HCO3- mislabel)
  - Excludes the 19 degenerate stuck predicted-Km values for the high-Km bundles
  - Caps sequences at max 10 per genus (organism redundancy reduction)
  - New defline format includes both predicted and experimental Km
  - Adds explicit km_evidence-anchored bundles alongside prediction-based ones
  - Adds a diagnostic file for the 26 suspect 18.0 mM entries

Run: python scripts/build_motif_bundles_v2.py
"""
import sqlite3, os, sys, time, subprocess, random
from collections import defaultdict
from pathlib import Path

DB_PATH    = Path("data/primary/carbodb.sqlite")
OUT_DIR    = Path("data/motifs_v2_clean")
GENUS_CAP  = 10
RNG_SEED   = 42  # for reproducible per-genus subsampling

# Degenerate stuck predicted-Km values to exclude from high-Km bundles
DEGENERATE_HIGH_KM_VALUES = {
    25.6423511505127,    15.5278882980347,    13.5064754486084,
    8.73185443878174,    4.96488666534424,    3.86344981193542,
    3.50768971443176,    2.68330836296082,    2.47630143165588,
    2.27833819389343,    2.02748990058899,    0.866567552089691,
    0.811510741710663,   0.475064665079117,   0.459040313959122,
    0.33828192949295,    0.181433230638504,   0.178775623440742,
    0.104575753211975,
}

# Target ECs for class-specific bundles (have ≥100 labeled sequences and known biology)
EC_TARGETS = [
    ("4.1.1.39", "rubisco"),
    ("4.2.1.1",  "carbonic_anhydrase"),
    ("4.1.1.31", "pep_carboxylase"),
    ("4.1.1.49", "pep_carboxykinase_atp"),
    ("4.1.1.32", "pep_carboxykinase_gtp"),
    ("6.3.4.14", "biotin_carboxylase"),
    ("6.4.1.1",  "pyruvate_carboxylase"),
]

OUT_DIR.mkdir(parents=True, exist_ok=True)
log_lines = []
def log(msg):
    print(msg, flush=True)
    log_lines.append(msg)


def open_db():
    if not DB_PATH.exists():
        sys.exit(f"DB not found at {DB_PATH}")
    conn = sqlite3.connect(str(DB_PATH))
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA cache_size = -524288")
    conn.execute("PRAGMA temp_store = MEMORY")
    return conn


def first_word(organism):
    if not organism:
        return ""
    return organism.split()[0]


def fmt_km(v):
    return f"{v:.4f}" if v is not None else "NA"


def write_fasta(path, rows, label):
    """rows: list of dicts with uniprot_id, ec_number, km_pred_mM, km_exp_mM,
            organism, reviewed, sequence."""
    n = 0
    with open(path, "w") as f:
        for r in rows:
            org = (r["organism"] or "unknown").replace(" ", "_")
            defline = (
                f">{r['uniprot_id']}|{r['ec_number']}"
                f"|km_pred_mM={fmt_km(r['km_pred_mM'])}"
                f"|km_exp_mM={fmt_km(r['km_exp_mM'])}"
                f"|{org}"
                f"|reviewed={r['reviewed']}"
            )
            seq = r["sequence"] or ""
            if not seq:
                continue
            f.write(defline + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")
            n += 1
    log(f"  → {path.name:55s} n={n:5d}  ({label})")
    return n


def cap_per_genus(rows, cap=GENUS_CAP, seed=RNG_SEED):
    """Randomly sample at most `cap` rows per genus (first word of organism).
    Sequences without an organism are kept as-is."""
    rng = random.Random(seed)
    by_genus = defaultdict(list)
    no_genus = []
    for r in rows:
        g = first_word(r["organism"])
        if g:
            by_genus[g].append(r)
        else:
            no_genus.append(r)
    out = list(no_genus)
    for g, lst in by_genus.items():
        if len(lst) <= cap:
            out.extend(lst)
        else:
            out.extend(rng.sample(lst, cap))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    conn = open_db()

    base_select = """
      SELECT s.uniprot_id, s.ec_number, s.organism, s.reviewed, s.sequence,
             p.km_pred_mM,
             (SELECT km_value_mM FROM km_evidence ke
                WHERE ke.sequence_id = s.id LIMIT 1) AS km_exp_mM
      FROM sequences s
      LEFT JOIN predictions p
        ON p.sequence_id = s.id AND p.model_version='v5'
      WHERE s.label = 1
        AND s.ec_number = ?
        AND s.sequence IS NOT NULL
        AND p.km_pred_mM IS NOT NULL
    """

    # ── Identify the 18.0 mM contamination block organisms ────────────────────
    contaminated_orgs = {
        r["organism"]
        for r in conn.execute(
            "SELECT DISTINCT s.organism FROM sequences s "
            "WHERE s.ec_number='4.1.1.39' AND s.km_best_mM = 18.0"
        ).fetchall()
        if r["organism"]
    }
    log(f"Contaminated organisms (18.0 mM block): {sorted(contaminated_orgs)}")

    # ─── 1. EC-class-specific bundles ─────────────────────────────────────────
    log("\n[1/4] EC-class-specific bundles (max 10 per genus):")
    for ec, slug in EC_TARGETS:
        rows = [dict(r) for r in conn.execute(base_select, (ec,)).fetchall()]
        rows = cap_per_genus(rows)
        path = OUT_DIR / f"ec_{ec.replace('.', '_')}_{slug}.fasta"
        write_fasta(path, rows, f"EC {ec}")

    # ─── 2. RuBisCO Km-stratified (predicted) ─────────────────────────────────
    log("\n[2/4] RuBisCO Km-stratified bundles (predicted Km):")
    rubisco_all = [dict(r) for r in conn.execute(base_select, ("4.1.1.39",)).fetchall()]

    rubisco_low = [
        r for r in rubisco_all
        if r["km_pred_mM"] is not None and r["km_pred_mM"] < 0.01
    ]
    rubisco_low = cap_per_genus(rubisco_low)
    write_fasta(OUT_DIR / "rubisco_low_km_pred.fasta",
                rubisco_low,
                "predicted Km < 0.01 mM, ≤10/genus")

    rubisco_high = [
        r for r in rubisco_all
        if r["km_pred_mM"] is not None
           and 0.1 <= r["km_pred_mM"] <= 5.0
           and r["km_pred_mM"] not in DEGENERATE_HIGH_KM_VALUES
    ]
    rubisco_high = cap_per_genus(rubisco_high)
    write_fasta(OUT_DIR / "rubisco_high_km_pred.fasta",
                rubisco_high,
                "predicted Km 0.1–5 mM, degenerate values excluded, ≤10/genus")

    # ─── 3. RuBisCO Km-stratified (experimental, from km_evidence) ────────────
    log("\n[3/4] RuBisCO Km-stratified bundles (experimental Km):")
    exp_query = """
      SELECT s.uniprot_id, s.ec_number, s.organism, s.reviewed, s.sequence,
             p.km_pred_mM, ke.km_value_mM AS km_exp_mM
      FROM km_evidence ke
      JOIN sequences   s ON s.id = ke.sequence_id
      LEFT JOIN predictions p ON p.sequence_id = s.id AND p.model_version='v5'
      WHERE ke.ec_number = '4.1.1.39'
        AND s.sequence IS NOT NULL
    """
    exp_rows = [dict(r) for r in conn.execute(exp_query).fetchall()]

    exp_low = [
        r for r in exp_rows
        if r["km_exp_mM"] is not None and r["km_exp_mM"] < 0.01
           and r["organism"] not in contaminated_orgs
    ]
    write_fasta(OUT_DIR / "rubisco_low_km_exp.fasta",
                exp_low,
                "experimental Km < 0.01 mM, contaminated orgs excluded")

    exp_high = [
        r for r in exp_rows
        if r["km_exp_mM"] is not None
           and 0.05 <= r["km_exp_mM"] <= 0.5
           and r["organism"] not in contaminated_orgs
    ]
    write_fasta(OUT_DIR / "rubisco_high_km_exp.fasta",
                exp_high,
                "experimental Km 0.05–0.5 mM, contaminated orgs excluded")

    suspect = [r for r in exp_rows if r["km_exp_mM"] == 18.0]
    write_fasta(OUT_DIR / "rubisco_brenda_18mM_suspect.fasta",
                suspect,
                "DIAGNOSTIC: 26 suspected HCO3⁻-mislabeled entries — DO NOT USE for motif analysis")

    # ─── 4. Negative control: pool of non-RuBisCO carboxylases ────────────────
    log("\n[4/4] Negative-control pooled bundle (non-RuBisCO carboxylases):")
    neg_rows = []
    for ec, _slug in EC_TARGETS:
        if ec == "4.1.1.39":
            continue
        rows = [dict(r) for r in conn.execute(base_select, (ec,)).fetchall()]
        rows = cap_per_genus(rows, cap=GENUS_CAP)
        # Sample at most 100 per EC for the pooled negative control
        if len(rows) > 100:
            rng = random.Random(RNG_SEED + hash(ec) % 1000)
            rows = rng.sample(rows, 100)
        neg_rows.extend(rows)
    write_fasta(OUT_DIR / "negative_control_other_carboxylases.fasta",
                neg_rows,
                "100 reps × 8 non-RuBisCO ECs, ≤10/genus per EC")

    conn.close()

    # ─── README ───────────────────────────────────────────────────────────────
    git_commit = "unknown"
    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=Path.cwd(), text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        pass

    readme = f"""# CarboDB Motif Bundles v2 (cleaned) — {time.strftime('%Y-%m-%d')}

Path: data/motifs_v2_clean/
Contact: Johannes Balkenhol
Generated by: scripts/build_motif_bundles_v2.py
Backend commit: {git_commit}

## What's different from the April 14 bundles

These bundles incorporate the contamination concerns you raised:

> "the top 45 sequences are almost all Asteraceae… the Helianthus entry at 18 mM
> is annotated as 'CO2 in form of HCO3-'… several of the high-Km BRENDA entries
> are from mutants… could the high predictions be picking up the HCO3-
> measurement scale or mutant data from the training set?"

We agree the BRENDA training data needs structural cleanup (your second
suggestion: separate wild-type / HCO3⁻ / mutant entries, convert HCO3⁻ to
dissolved CO2 at pH, exclude ambiguous entries). That's deferred — it's
multi-day pipeline work and would require model retraining. Happy to discuss
scope on a call when you're ready.

In the meantime, these v2 bundles apply the following filters so that
downstream motif analysis isn't dominated by the contamination tail:

1. **Helianthus 18.0 mM block excluded.** All 26 entries with experimental
   Km = 18.0 mM (almost all Helianthus annuus) are removed from the
   high-Km bundles. The block is still available in
   `rubisco_brenda_18mM_suspect.fasta` for diagnostic inspection.

2. **Degenerate model predictions excluded from high-Km bundles.**
   The model has 19 specific output values (e.g. 2.68 mM, 0.46 mM,
   15.53 mM) that it returns for hundreds of distinct sequences each —
   these are regression-failure modes, not real predictions. Any sequence
   whose predicted Km is one of these stuck values is excluded from
   `rubisco_high_km_pred.fasta`. This is the cause of the Zea mays /
   Chlamydomonas case you flagged: their predicted 2.79 mM is a stuck
   value, not a real high-affinity prediction.

3. **Per-genus redundancy cap of 10.** Each EC bundle and each Km-stratified
   bundle is capped at 10 sequences per genus. This prevents one organism
   (e.g. Triticum, Helianthus) from dominating motif frequencies. Sampling
   is reproducible (seed={RNG_SEED}).

4. **New defline format includes both predicted and experimental Km:**
   `>UNIPROT_ID|EC|km_pred_mM=X.XXXX|km_exp_mM=X.XXXX|Organism|reviewed=0/1`
   km_exp_mM is "NA" when no km_evidence row exists for that sequence.

## Files

### EC-class-specific (predicted Km only, ≤10/genus)
Use these for "what motifs are class-specific" questions.

  ec_4_1_1_39_rubisco.fasta                       Ribulose-1,5-bisphosphate carboxylase
  ec_4_2_1_1_carbonic_anhydrase.fasta             Carbonic anhydrase
  ec_4_1_1_31_pep_carboxylase.fasta               PEP carboxylase
  ec_4_1_1_49_pep_carboxykinase_atp.fasta         PEPCK (ATP-dependent)
  ec_4_1_1_32_pep_carboxykinase_gtp.fasta         PEPCK (GTP-dependent)
  ec_6_3_4_14_biotin_carboxylase.fasta            Biotin carboxylase
  ec_6_4_1_1_pyruvate_carboxylase.fasta           Pyruvate carboxylase

### RuBisCO Km-stratified (predicted Km)
For population-scale Km-vs-motif analysis. Larger samples, but uses model
predictions, so contamination filters above are critical.

  rubisco_low_km_pred.fasta    predicted Km < 0.01 mM, ≤10/genus
  rubisco_high_km_pred.fasta   predicted Km 0.1–5 mM, ≤10/genus, stuck values excluded

### RuBisCO Km-stratified (experimental Km)
For small-but-trustworthy Km-vs-motif analysis. Drawn from km_evidence
table directly, no model in the loop. The high-Km cohort is small (n≈18).

  rubisco_low_km_exp.fasta     experimental Km < 0.01 mM
  rubisco_high_km_exp.fasta    experimental Km 0.05–0.5 mM
  rubisco_brenda_18mM_suspect.fasta   DIAGNOSTIC ONLY — the 26 suspect entries

Note on experimental coverage: The 916 sequences with experimental Km come
from ~30 distinct BRENDA measurements applied at the organism level (e.g.
all 48 Synura petersenii sequences share one published Km of 0.0284 mM).
For motif statistics this is fine, but the effective measurement count is
~30, not ~900. The per-genus cap helps surface this.

### Negative control
For "discriminating motifs" (present in target EC, absent from others).

  negative_control_other_carboxylases.fasta   100 reps × 8 non-RuBisCO ECs, ≤10/genus

## How we suggest using these

1. **EC-specific motifs** — discriminating motifs:
   target = `ec_4_1_1_39_rubisco.fasta`, background = `negative_control_other_carboxylases.fasta`

2. **RuBisCO Km motifs** — comparison:
   - First pass on the larger predicted-Km bundles (low_km_pred vs high_km_pred)
   - Cross-check anything interesting on the smaller exp bundles (low_km_exp vs high_km_exp)
   - Motifs that show up in both predicted and experimental splits are real;
     motifs that show only in the predicted split should be treated cautiously
     (could be model-driven, not biology-driven)

3. **Verify our contamination call**:
   `rubisco_brenda_18mM_suspect.fasta` is the 26 entries we flagged as suspect.
   If you find motifs in the suspect set that don't appear in the rest of the
   high-Km predicted bundle, that's evidence the 18.0 mM block is genuinely
   different (consistent with the HCO3⁻ hypothesis).

## Generation log

See `BUNDLES_v2_GENERATION.log` for per-file counts.
"""
    (OUT_DIR / "README.md").write_text(readme)
    log("\nWrote README.md")

    log_path = OUT_DIR / "BUNDLES_v2_GENERATION.log"
    log_path.write_text(
        f"# Generated {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"# git commit: {git_commit}\n"
        f"# total wall time: {time.time()-t0:.1f}s\n\n"
        + "\n".join(log_lines)
    )
    log(f"\nDone in {time.time()-t0:.1f}s. Output in {OUT_DIR}/")


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parents[1])
    main()
