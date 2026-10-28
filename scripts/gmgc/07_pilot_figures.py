#!/usr/bin/env python3
"""
Sanity-check figures for the GMGC marine pilot.

Reads:  data/gmgc/results/predictions_joined.parquet
Writes: data/gmgc/figures/*.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT  = Path.home() / "Projects_shared/CarboDB_v3"
PRED  = ROOT / "data/gmgc/results/predictions_joined.parquet"
FIGS  = ROOT / "data/gmgc/figures"
FIGS.mkdir(parents=True, exist_ok=True)

df = pd.read_parquet(PRED)
print(f"Loaded {len(df):,} predictions", file=sys.stderr)

carb = df[df["is_carboxylase"]].copy()
print(f"Carboxylase positives: {len(carb):,}", file=sys.stderr)

# Figure 1: Binary probability histogram
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(df["bin_prob"], bins=50, color="steelblue", edgecolor="white")
ax.axvline(0.5, ls="--", c="red", lw=1, label="threshold 0.5")
ax.set_xlabel("v5 binary probability (carboxylase)")
ax.set_ylabel("Pfam-filtered candidates (count)")
ax.set_title(f"Binary classifier output, n={len(df):,}")
ax.legend()
plt.tight_layout()
plt.savefig(FIGS / "01_binary_prob_hist.png", dpi=120)
plt.close()
print("  fig 01 saved", file=sys.stderr)

# Figure 2: EC class distribution
ec_counts = carb["ec_pred"].value_counts().head(20)
fig, ax = plt.subplots(figsize=(10, 6))
ax.barh(ec_counts.index[::-1], ec_counts.values[::-1], color="seagreen")
ax.set_xlabel("Count")
ax.set_title(f"Top 20 EC classes — predicted carboxylases (n={len(carb):,})")
plt.tight_layout()
plt.savefig(FIGS / "02_ec_class_dist.png", dpi=120)
plt.close()
print("  fig 02 saved", file=sys.stderr)

# Figure 3: Km distribution
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(carb["km_pred_log10"], bins=50, color="darkorange", edgecolor="white")
ax.set_xlabel("log10(Km mM) — predicted")
ax.set_ylabel("Predicted carboxylases (count)")
ax.set_title(f"Predicted Km distribution, n={len(carb):,}")
if len(carb) > 0:
    med = carb["km_pred_log10"].median()
    ax.axvline(med, ls="--", c="red", lw=1,
               label=f"median = {10**med:.3f} mM")
    ax.legend()
plt.tight_layout()
plt.savefig(FIGS / "03_km_distribution.png", dpi=120)
plt.close()
print("  fig 03 saved", file=sys.stderr)

# Figure 4: World map
gps = carb.dropna(subset=["latitude", "longitude"]).copy()
print(f"  {len(gps):,} predicted carboxylases with GPS for map", file=sys.stderr)

if len(gps) > 0:
    site = gps.groupby(["latitude", "longitude"]).agg(
        n_carboxy=("seq_id", "count"),
        median_km_mM=("km_pred_mM", "median"),
    ).reset_index()

    fig, ax = plt.subplots(figsize=(12, 6))
    sc = ax.scatter(site["longitude"], site["latitude"],
                    c=np.log10(site["median_km_mM"]),
                    s=4 * np.sqrt(site["n_carboxy"]),
                    cmap="coolwarm", alpha=0.8, edgecolor="black", lw=0.3)
    cb = plt.colorbar(sc, ax=ax, label="log10(median predicted Km mM)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"GMGC marine carboxylases — {len(gps):,} predictions at {len(site)} sites")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(FIGS / "04_world_map_km.png", dpi=120)
    plt.close()
    print("  fig 04 saved", file=sys.stderr)

# Figure 5: Km vs latitude
if len(gps) > 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(gps["latitude"], gps["km_pred_log10"], s=5, alpha=0.3, c="darkblue")
    ax.set_xlabel("Latitude")
    ax.set_ylabel("log10(Km mM)")
    ax.set_title("Predicted Km vs latitude")
    plt.tight_layout()
    plt.savefig(FIGS / "05_km_vs_latitude.png", dpi=120)
    plt.close()
    print("  fig 05 saved", file=sys.stderr)

# Figure 6: Top genera
if "tax_name" in carb.columns and "tax_rank" in carb.columns:
    genera = carb[carb["tax_rank"] == "genus"]["tax_name"].value_counts().head(30)
    if len(genera) > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.barh(genera.index[::-1], genera.values[::-1], color="purple")
        ax.set_xlabel("Predicted carboxylases (count)")
        ax.set_title("Top 30 genera — predicted marine carboxylases")
        plt.tight_layout()
        plt.savefig(FIGS / "06_top_genera.png", dpi=120)
        plt.close()
        print("  fig 06 saved", file=sys.stderr)

print(f"\nAll figures written to {FIGS}", file=sys.stderr)
