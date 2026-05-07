# Motif Analysis — Status and Next Steps

What the colleague's first analysis found, what the v2 cleanup fixed,
what the **v3 within-Form-I bundle** adds, and the protocol for evaluating
amino-acid substitution effects on Km.

---

## Where things stand

### Analysis 1 — colleague's first run (April 14)

**Bundle:** `data/motifs_v1/` (now superseded). Took the top high-Km
RuBisCO predictions from CarboDB v3 and ran motif/MEME analysis vs the
low-Km set.

**Result:** found a Helianthus 18 mM cluster as the dominant high-Km
"signal." On inspection (May 5) this was contamination:
- 26 distinct UniProt IDs all sharing exactly 18.0 mM
- BRENDA annotation says "CO2 in form of HCO3-"
- Plus mutant entries (D117H, K122R, Met330Leu) treated as wild-type
- Plus 19 distinct XGBoost stuck values (e.g. 15.5278882980347 mM repeated
  across 400 sequences) showing up as "high-Km" in the prediction set

**Conclusion:** the v1 bundle was unreliable. We rebuilt as v2.

### v2 bundle — May 5 cleanup

**Where:** `data/motifs_v2_clean/`, also tarball `motifs_v2_clean.tar.gz`.

**Generator:** `scripts/build_motif_bundles_v2.py` (seed=42).

**Contents (13 FASTAs + README + log):**

| File | Sequences | Notes |
|---|---|---|
| `ec_4_1_1_39_rubisco.fasta` | 32K | All RuBisCO, max 10 per genus |
| `ec_4_1_1_31_pep_carboxylase.fasta` | 6K | PEPC |
| `ec_4_1_1_32_pep_carboxykinase_gtp.fasta` | 4K | PEPCK-GTP |
| `ec_4_1_1_49_pep_carboxykinase_atp.fasta` | 5K | PEPCK-ATP |
| `ec_4_2_1_1_carbonic_anhydrase.fasta` | 9K | Carbonic anhydrase |
| `ec_6_3_4_14_biotin_carboxylase.fasta` | 14K | Biotin carboxylase |
| `ec_6_4_1_1_pyruvate_carboxylase.fasta` | 9K | Pyruvate carboxylase |
| `rubisco_low_km_pred.fasta` | 32,067 | predicted Km < 0.01 mM |
| `rubisco_high_km_pred.fasta` | 2,192 | predicted Km 0.1–5 mM, 19 stuck values excluded |
| `rubisco_low_km_exp.fasta` | 372 | experimental Km < 0.01 mM, BRENDA-tier-A |
| `rubisco_high_km_exp.fasta` | 18 | experimental Km > 0.1 mM, BRENDA-tier-A |
| `rubisco_brenda_18mM_suspect.fasta` | 26 | diagnostic — the contamination set |
| `negative_control_other_carboxylases.fasta` | 600 | non-RuBisCO carboxylases |

**Colleague's analysis 2 result on v2:**

11 amino-acid positions differ significantly between low-Km and high-Km
RuBisCO (Bonferroni-corrected p < 0.05):

| Position | Low-Km | High-Km | p-value | Region |
|---|---|---|---|---|
| 86 | R | H | 1.4×10⁻³⁷ | N-terminal domain |
| 251 | Y | I | 2.5×10⁻⁴ | near active site |
| 255 | E | V | 4.6×10⁻⁶ | near active site |
| 258 | K | R | 1.9×10⁻¹⁷ | near active site |
| 391 | V | T | 1.2×10⁻²⁸ | C-terminal |
| 445 | L | I | sig. | C-terminal cluster |
| 461 | L | C | sig. | C-terminal cluster |
| 467 | Q | E | sig. | C-terminal cluster |
| 91 | A | P | 2.4×10⁻³³ | high-Km enriched |
| 145 | S | V | 7.4×10⁻²² | high-Km enriched |
| 449 | A | S | 5.3×10⁻¹⁵ | high-Km enriched |

**Key observations:**
- Catalytic residues (K201, D203, E204) and specificity residues (A375,
  S376) are **invariant** — they don't differ between Km classes.
- Differences cluster spatially around the active site (positions
  251–258) and at the subunit interface (445–467).
- Panel A of the report shows **low-Km RuBisCO is more conserved overall**
  than high-Km, consistent with stronger purifying selection.

### Cross-check against v5 SHAP feature importance

The Analysis page's per-EC SHAP plot for RuBisCO (EC 4.1.1.39) shows:

| Rank | Feature | Group | Importance |
|---|---|---|---|
| 1 | PF02788 (RuBisCO_large_N) | Pfam | 26.12% |
| 2 | PF00016 (RuBisCO_large) | Pfam | 21.92% |
| 3 | ESM-2 dim 1083 | ESM-2 | 6.40% |
| 4 | dipeptide YK | Dipeptide | 2.82% |
| 5 | ESM-2 dim 1059 | ESM-2 | 2.58% |
| 6 | dipeptide QP | Dipeptide | 2.39% |
| 7 | ESM-2 dim 448 | ESM-2 | 2.37% |
| 8 | dipeptide TD | Dipeptide | 2.09% |
| 9 | AAC Y | Composition | 1.56% |
| 10 | dipeptide WT | Dipeptide | 1.47% |

Group totals: Pfam 49.3%, ESM-2 44.1%, Dipeptide 6.7%.

**Interpretation matched against the motif report:**
- **Pfam carries enzyme identity** (both PF00016 and PF02788 are present
  in *every* RuBisCO regardless of Km — they discriminate carboxylase
  vs not, not Km class). The Pfam **e-value** is a numeric feature too —
  lower e-value means stronger HMM match, which under purifying
  selection (Panel A) correlates with low Km. So Pfam's 49% is partly
  identity and partly conservation-score-as-Km-proxy.
- **ESM-2 carries the actual Km signal.** 44% importance, no interpretable
  per-feature meaning, but the residue-level changes the motif analysis
  finds (positions 251, 255, 258, etc.) are exactly the kind of long-range
  context ESM-2 picks up well. Specifically: hydrophobic↔charged
  substitutions (V255E), polar↔hydrophobic (T391V, S145V), backbone
  flexibility (A91P) all shift embedding coordinates measurably.
- **Two of four top dipeptides match motif-report substitutions:**
  - **YK** ← `I251Y` substitution near active site (if 252 is K, this
    creates a YK pair that low-Km has and high-Km lacks)
  - **QP** ← `A91P` substitution (if 90 is Q, creates QP that high-Km has)
  - TD and WT don't have a clean motif-report explanation; likely Form I
    vs Form II/III phylogenetic background composition signal.

**Caveat:** YK/QP confirmation requires checking what the actual neighboring
residues are at positions 90 and 252 in the alignment. The logic is
consistent but a quick alignment-grep would settle it definitively.

---

## The phylogenetic confound — and v3

The v2 analysis's 11 residue changes have a known confound (the report's
"Concerning Question 1"):

> Are we seeing real Km-associated changes, or just Form I vs Form II/III
> taxonomic differences? The low-Km group is likely enriched in C4 grasses
> and the high-Km group in archaea/bacteria.

This matters because the differences could reflect lineage-specific
sequence drift that incidentally also correlates with Km, rather than
selection on Km per se.

**v3 bundle: within-Form-I split by Km.**

Generator: `scripts/build_motif_bundles_v3_form_split.py` (this PR).

Contents:
- `form_I_low_km_pred.fasta` — predicted Km < 0.01 mM, restricted to Form I
- `form_I_high_km_pred.fasta` — predicted Km 0.1–5 mM, restricted to Form I
  (still excluding the 19 stuck values and the Helianthus contamination)
- `form_II_III_pooled.fasta` — Form II and Form III together for contrast
- `README.md` documenting the Form classification logic

Form classification uses InterPro PANTHER family identifiers (PTHR42704
distinguishes the Form-I large chain, derived via `features_interpro`).
Where InterPro is uninformative we fall back to organism kingdom.

**Test of the v2 finding:**
- If the 251/255/258 cluster, the 86 H→R, the 391 T→V remain significant
  within Form I: the signal is Km-related, not phylogenetic. Strong claim.
- If they vanish: the v2 finding was driven by Form I/II/III divergence.
  Reframe the result.

---

## Strategy: characterizing AA-Km associations and motif localization

For evaluating each of the 11 candidate residues (or any new residues the
v3 analysis surfaces), the protocol is:

### Step 1 — Confirm Km association is Form-I-specific

Repeat Fisher's exact test position-by-position **within Form I only**.
Apply Bonferroni correction over the alignment length (~470 positions).
Report:
- Position
- Low-Km dominant AA
- High-Km dominant AA
- p-value within Form I (compare to v2's whole-RuBisCO p-value)
- Survival outcome: kept / lost / weakened

Likely-survivors based on biochemistry: positions in the active-site
neighborhood (251, 255, 258) and the C-terminal subunit interface
(445, 461, 467). Likely-vanishers: positions far from the active site
that may simply track Form I/II/III divergence.

### Step 2 — Map the surviving positions onto 1RCX

In PyMOL:
```
fetch 1RCX
select cat, resi 201+203+204 and chain A   # catalytic (invariant)
select spec, resi 375+376 and chain A       # specificity (invariant)
select km_assoc, resi 86+251+255+258+391+445+461+467 and chain A
color grey80, all
color magenta, cat
color cyan, spec
color yellow, km_assoc
show sticks, cat | spec | km_assoc
```

Visual question: are the Km-associated residues forming a distinct
"second shell" around the catalytic core? Are they clustered or
distributed?

### Step 3 — Categorize the substitutions

For each surviving position:

| Category | Example | Mechanism |
|---|---|---|
| Charge change | V255E (low-Km) | introduces salt bridge / changes electrostatic complementarity at active site |
| Hydrophobicity change | T391V (low-Km), S145V (high-Km) | shifts local packing / dynamics |
| Backbone flexibility | A91P (high-Km) | proline restricts ϕ angle; could destabilize a transient state |
| Conservative | R258K (low-Km) | minor — both basic; may matter for specific salt bridge geometry |
| Steric | I251Y (low-Km) | adds aromatic ring that may engage substrate |

Group results into a table for the writeup. The mechanistic categories
should hint at *which step* of catalysis is being tuned (substrate
binding vs transition-state stabilization vs product release).

### Step 4 — Are they within a known motif?

Check overlap with:
- The 7 expert motifs in `features_expert_motifs` (RuBisCO K-K,
  RuBisCO G-K, CA H-H, CA His cluster, PEPC R-R, Biotin M-K, Biotin A-M-K)
- Pfam domain spans (PF00016 RuBisCO_large covers most of the catalytic
  domain; PF02788 RuBisCO_large_N covers the N-terminal beta barrel)
- FIMO motif hits in `features_fimo` (data-driven, MEME-derived)
- Loop 6 (residues 332–338 in spinach) — known specificity-determining
  region

For each Km-associated position, report which motif/domain it falls in
(or whether it's in a "novel" un-motif-mapped region).

### Step 5 — Co-occurrence

Are the changes independent or do they always travel together? For the
8 surviving positions (hypothetical), build the 2⁸ pattern table from
the alignment, look at which patterns are populated:

- If `(low-Km AAs) at all 8 positions` and `(high-Km AAs) at all 8`
  together account for >80% of sequences: the changes are *coupled*,
  almost certainly phylogenetic.
- If many intermediate patterns exist with comparable frequencies:
  the changes are *independent* — stronger evidence each is under
  selection on its own.

The distinction matters for any "engineer this single residue and
shift Km" follow-up: only independent positions are worth single-mutation
experiments.

### Step 6 — Tie back to the SHAP picture

Once the surviving Form-I residues are known, return to the SHAP feature
list and check:
- Are the dipeptides explained? (YK, QP, TD, WT — which match a surviving
  position pair?)
- Does the Pfam e-value still correlate with Form-I-only Km after
  controlling for the per-residue features?
- For ESM-2 dim 1083, dim 1059, dim 448 — do these correlate with the
  surviving residues across the alignment? (Approach: sort sequences by
  ESM-2 dim 1083 value, look at the AA distribution at position 251 in
  the top vs bottom decile. This is a coarse interpretability probe.)

A clean writeup for the meeting/paper would have the v3 analysis
answer Concerning Question 1, the structural mapping showing the
"second shell" architecture, and the SHAP cross-reference confirming
the model has internalized the same biology the alignment finds.

---

## What the colleague should run on the v3 bundle

In their preferred environment (the same MAFFT L-INS-i + Fisher's exact
they used for v2):

```bash
# 1. Within-Form-I differential
mafft --localpair --maxiterate 1000 form_I_low_km_pred.fasta > form_I_low.aln
mafft --localpair --maxiterate 1000 form_I_high_km_pred.fasta > form_I_high.aln
# Then reuse the existing Fisher-test pipeline.

# 2. Form II/III sanity check (do the v2 positions show the OPPOSITE
#    or the SAME pattern in Form II/III? If opposite, strong evidence
#    that v2 was phylogenetic.)
mafft --localpair --maxiterate 1000 form_II_III_pooled.fasta > form_II_III.aln
# Inspect the v2 11 positions in this alignment.

# 3. Co-occurrence pattern table
# Custom script — see the README in motifs_v3_form_split/ for an
# example skeleton.
```

Estimated time: ~1-2 days of analysis once the bundles are in hand.

---

## Open questions for the writeup

1. **Why are catalytic residues invariant if RuBisCO Km varies 1000-fold
   between Form I and Form II?** The textbook answer: catalysis is an
   ancient capability that selection won't touch; what selection tunes
   is the dynamics around the catalytic core. The v3 results should
   either confirm this or surface a counterexample.

2. **Is the C-terminal cluster (445–467) a subunit-assembly effect
   rather than a Km effect?** RuBisCO is L₈S₈; mutations near the L–S
   interface could affect assembly stoichiometry rather than per-monomer
   kinetics. Worth a structural-biology consultation.

3. **Does the model really learn the same biology?** The fact that the
   SHAP top features (Pfam e-value, dipeptides YK/QP, ESM-2 dims) match
   the alignment-derived signal suggests yes, but a formal causal test
   (mask the residues at positions 251–258 in input → measure ΔKm
   prediction) would be much stronger evidence.
