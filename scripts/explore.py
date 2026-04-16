import sqlite3, pandas as pd, numpy as np
from scipy import stats

conn = sqlite3.connect('data/primary/carbodb.sqlite')

# ── 1. PARASITIC PLANTS ───────────────────────────────────────────────────
print("=== PARASITIC PLANTS ===")
q = pd.read_sql("""
    SELECT s.organism, p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Cuscuta%' THEN 'Holoparasite'
             WHEN s.organism LIKE 'Orobanche%' THEN 'Holoparasite'
             WHEN s.organism LIKE 'Phelipanche%' THEN 'Holoparasite'
             WHEN s.organism LIKE 'Viscum%' THEN 'Hemiparasite'
             WHEN s.organism LIKE 'Striga%' THEN 'Hemiparasite'
             WHEN s.organism LIKE 'Rhinanthus%' THEN 'Hemiparasite'
        END as ptype
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Cuscuta%' OR s.organism LIKE 'Orobanche%'
         OR s.organism LIKE 'Viscum%' OR s.organism LIKE 'Striga%'
         OR s.organism LIKE 'Rhinanthus%' OR s.organism LIKE 'Phelipanche%')
""", conn)
free = pd.read_sql("""
    SELECT p.km_pred_mM FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Ipomoea%' OR s.organism LIKE 'Nicotiana%'
         OR s.organism LIKE 'Solanum%') LIMIT 300
""", conn)
print(q.groupby(['ptype','organism']).agg(
    n=('km_pred_mM','count'),
    mean_uM=('km_pred_mM', lambda x: round(x.mean()*1000,1))
).to_string())
u,p = stats.mannwhitneyu(q['km_pred_mM'], free['km_pred_mM'])
print(f"Parasites mean: {q.km_pred_mM.mean()*1000:.1f} uM | Free-living: {free.km_pred_mM.mean()*1000:.1f} uM | p={p:.4f}")

# ── 2. CROP vs WILD RELATIVE ──────────────────────────────────────────────
print("\n=== CROP vs WILD RELATIVE Km ===")
crops = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Zea mays%' THEN 'Maize_crop'
             WHEN s.organism LIKE 'Zea%' THEN 'Zea_wild'
             WHEN s.organism LIKE 'Oryza sativa%' THEN 'Rice_crop'
             WHEN s.organism LIKE 'Oryza%' THEN 'Oryza_wild'
             WHEN s.organism LIKE 'Hordeum vulgare%' THEN 'Barley_crop'
             WHEN s.organism LIKE 'Hordeum%' THEN 'Hordeum_wild'
             WHEN s.organism LIKE 'Helianthus annuus%' THEN 'Sunflower_crop'
             WHEN s.organism LIKE 'Helianthus%' THEN 'Helianthus_wild'
             WHEN s.organism LIKE 'Solanum tuberosum%' THEN 'Potato_crop'
             WHEN s.organism LIKE 'Solanum%' THEN 'Solanum_wild'
             WHEN s.organism LIKE 'Triticum aestivum%' THEN 'Wheat_6x'
             WHEN s.organism LIKE 'Triticum durum%' THEN 'Wheat_4x'
             WHEN s.organism LIKE 'Triticum%' THEN 'Triticum_wild'
        END as grp
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Zea%' OR s.organism LIKE 'Oryza%'
         OR s.organism LIKE 'Hordeum%' OR s.organism LIKE 'Helianthus%'
         OR s.organism LIKE 'Solanum%' OR s.organism LIKE 'Triticum%')
""", conn).dropna(subset=['grp'])
print(crops.groupby('grp').agg(
    n=('km_pred_mM','count'),
    mean_uM=('km_pred_mM', lambda x: round(x.mean()*1000,1)),
    min_uM=('km_pred_mM', lambda x: round(x.min()*1000,1)),
    max_uM=('km_pred_mM', lambda x: round(x.max()*1000,1))
).sort_values('mean_uM').to_string())

# ── 3. CONVERGENT C4 ──────────────────────────────────────────────────────
print("\n=== CONVERGENT C4 — same Km across independent lineages? ===")
c4 = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Zea%' OR s.organism LIKE 'Sorghum%'
                  OR s.organism LIKE 'Cynodon%' OR s.organism LIKE 'Chloris%'
                  OR s.organism LIKE 'Sporobolus%' OR s.organism LIKE 'Muhlenbergia%'
                  THEN 'Poaceae_C4'
             WHEN s.organism LIKE 'Amaranthus%' OR s.organism LIKE 'Atriplex%'
                  OR s.organism LIKE 'Suaeda%' THEN 'Caryophyllales_C4'
             WHEN s.organism LIKE 'Portulaca%' THEN 'Portulacaceae_C4'
             WHEN s.organism LIKE 'Flaveria%' THEN 'Asteraceae_C4'
             WHEN s.organism LIKE 'Cleome%' THEN 'Cleomaceae_C4'
        END as lineage
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Zea%' OR s.organism LIKE 'Sorghum%'
         OR s.organism LIKE 'Cynodon%' OR s.organism LIKE 'Chloris%'
         OR s.organism LIKE 'Sporobolus%' OR s.organism LIKE 'Muhlenbergia%'
         OR s.organism LIKE 'Amaranthus%' OR s.organism LIKE 'Atriplex%'
         OR s.organism LIKE 'Portulaca%' OR s.organism LIKE 'Flaveria%'
         OR s.organism LIKE 'Cleome%' OR s.organism LIKE 'Suaeda%')
""", conn).dropna(subset=['lineage'])
print(c4.groupby('lineage').agg(
    n=('km_pred_mM','count'),
    mean_uM=('km_pred_mM', lambda x: round(x.mean()*1000,2)),
    median_uM=('km_pred_mM', lambda x: round(x.median()*1000,2))
).sort_values('mean_uM').to_string())
groups = [c4[c4['lineage']==l]['km_pred_mM'].values for l in c4['lineage'].unique() if len(c4[c4['lineage']==l])>=3]
if len(groups)>=3:
    f,p = stats.kruskal(*groups)
    print(f"Kruskal-Wallis p={p:.4f} — {'CONVERGED (same Km)' if p>0.05 else 'DIFFER (independent solutions)'}")

# ── 4. SYMBIODINIACEAE ────────────────────────────────────────────────────
print("\n=== CORAL SYMBIONT CLADES ===")
sym = pd.read_sql("""
    SELECT s.organism, p.km_pred_mM
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Symbiodinium%' OR s.organism LIKE 'Cladocopium%'
         OR s.organism LIKE 'Durusdinium%' OR s.organism LIKE 'Breviolum%'
         OR s.organism LIKE 'Fugacium%')
""", conn)
if len(sym)>0:
    sym['genus'] = sym['organism'].apply(lambda x: x.split()[0])
    print(sym.groupby('genus').agg(
        n=('km_pred_mM','count'),
        mean_uM=('km_pred_mM', lambda x: round(x.mean()*1000,2))
    ).sort_values('mean_uM').to_string())
    print("Durusdinium=heat tolerant | Cladocopium=bleaching sensitive")
else:
    print("No Symbiodiniaceae in RuBisCO predictions")

# ── 5. PLOIDY ─────────────────────────────────────────────────────────────
print("\n=== PLOIDY EFFECT ===")
ploidy = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Triticum aestivum%' THEN 'Wheat_6x'
             WHEN s.organism LIKE 'Triticum durum%' THEN 'Wheat_4x'
             WHEN s.organism LIKE 'Triticum urartu%' THEN 'Wheat_diploid_ancestor'
             WHEN s.organism LIKE 'Aegilops%' THEN 'Wheat_diploid_ancestor'
             WHEN s.organism LIKE 'Avena sativa%' THEN 'Oat_6x'
             WHEN s.organism LIKE 'Avena strigosa%' THEN 'Oat_diploid'
             WHEN s.organism LIKE 'Gossypium hirsutum%' THEN 'Cotton_4x'
             WHEN s.organism LIKE 'Gossypium arboreum%' THEN 'Cotton_2x'
             WHEN s.organism LIKE 'Gossypium raimondii%' THEN 'Cotton_2x'
             WHEN s.organism LIKE 'Oryza sativa%' THEN 'Rice_2x'
             WHEN s.organism LIKE 'Hordeum vulgare%' THEN 'Barley_2x'
        END as ploidy_group
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.1.1.39' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Triticum%' OR s.organism LIKE 'Aegilops%'
         OR s.organism LIKE 'Avena%' OR s.organism LIKE 'Gossypium%'
         OR s.organism LIKE 'Oryza sativa%' OR s.organism LIKE 'Hordeum vulgare%')
""", conn).dropna(subset=['ploidy_group'])
print(ploidy.groupby('ploidy_group').agg(
    n=('km_pred_mM','count'),
    mean_uM=('km_pred_mM', lambda x: round(x.mean()*1000,1)),
    median_uM=('km_pred_mM', lambda x: round(x.median()*1000,1))
).sort_values('mean_uM').to_string())

# ── 6. CA — temperature optima effect ────────────────────────────────────
print("\n=== CA Km by organism temperature optima ===")
ca_t = pd.read_sql("""
    SELECT p.km_pred_mM,
        CASE WHEN s.organism LIKE 'Sulfolobus%' OR s.organism LIKE 'Thermococcus%'
                  OR s.organism LIKE 'Pyrococcus%' OR s.organism LIKE 'Thermus%'
                  THEN 'Thermophile_70-90C'
             WHEN s.organism LIKE 'Psychrobacter%' OR s.organism LIKE 'Colwellia%'
                  OR s.organism LIKE 'Polaribacter%' THEN 'Psychrophile_0-15C'
             WHEN s.organism LIKE 'Homo%' OR s.organism LIKE 'Mus%'
                  OR s.organism LIKE 'Bos%' THEN 'Mesophile_animal_37C'
             WHEN s.organism LIKE 'Escherichia%' OR s.organism LIKE 'Bacillus subtilis%'
                  THEN 'Mesophile_bact_37C'
             WHEN s.organism LIKE 'Chlamydomonas%' OR s.organism LIKE 'Chlorella%'
                  THEN 'Mesophile_alga_20C'
        END as temp_group
    FROM sequences s JOIN predictions p ON p.sequence_id=s.id
    WHERE s.label=1 AND s.ec_number='4.2.1.1' AND p.km_pred_mM IS NOT NULL
    AND (s.organism LIKE 'Sulfolobus%' OR s.organism LIKE 'Thermococcus%'
         OR s.organism LIKE 'Pyrococcus%' OR s.organism LIKE 'Thermus%'
         OR s.organism LIKE 'Psychrobacter%' OR s.organism LIKE 'Colwellia%'
         OR s.organism LIKE 'Polaribacter%' OR s.organism LIKE 'Homo%'
         OR s.organism LIKE 'Mus%' OR s.organism LIKE 'Bos%'
         OR s.organism LIKE 'Escherichia%' OR s.organism LIKE 'Bacillus subtilis%'
         OR s.organism LIKE 'Chlamydomonas%' OR s.organism LIKE 'Chlorella%')
""", conn).dropna(subset=['temp_group'])
print(ca_t.groupby('temp_group').agg(
    n=('km_pred_mM','count'),
    mean_mM=('km_pred_mM', lambda x: round(x.mean(),3)),
    median_mM=('km_pred_mM', lambda x: round(x.median(),3))
).sort_values('mean_mM').to_string())
print("Expected: psychrophiles lowest Km (more CO2 dissolved in cold)")
print("          thermophiles highest (CO2 less soluble at high T)")

conn.close()
print("\n=== ALL DONE ===")
