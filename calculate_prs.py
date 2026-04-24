"""
Step 3: Calculate Polygenic Risk Score + Visualize
===================================================
Author: Hafsah Shamsi
GitHub: github.com/HafsahShamsi

PRS = sum of (effect_size * genotype_dosage) across all SNPs
Effect sizes are log(OR) — so we take log of the OR values first.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

os.makedirs("./outputs", exist_ok=True)

GOLD       = "#C9A84C"
DARK       = "#1A1A2E"
MUTED      = "#4A4A6A"
HIGHLIGHT  = "#E05C5C"
BLUE       = "#5C9BE0"

def calculate_prs(geno_file, snp_file):
    print("Loading data...")
    geno = pd.read_csv(geno_file, index_col=0)
    snps = pd.read_csv(snp_file)

    # Align SNPs to genotype columns
    snps = snps[snps["rsid"].isin(geno.columns)].copy()
    snps = snps.set_index("rsid")
    geno = geno[snps.index]

    print(f"  {len(snps)} SNPs matched, {len(geno)} individuals")

    # Effect size: OR -> log(OR) = beta
    snps["beta"] = np.log(snps["effect_size"].clip(0.01))

    # PRS = genotype matrix (n x p) dot beta vector (p,)
    betas = snps["beta"].values
    prs = geno.values @ betas

    prs_df = pd.DataFrame({
        "individual_id": geno.index,
        "prs": prs
    }).set_index("individual_id")

    # Standardize to mean=0, sd=1
    prs_df["prs_z"] = (prs_df["prs"] - prs_df["prs"].mean()) / prs_df["prs"].std()

    # Percentile
    prs_df["percentile"] = prs_df["prs_z"].rank(pct=True) * 100

    print(f"\nPRS Summary:")
    print(f"  Mean:   {prs_df['prs_z'].mean():.4f}")
    print(f"  SD:     {prs_df['prs_z'].std():.4f}")
    print(f"  Min:    {prs_df['prs_z'].min():.4f}")
    print(f"  Max:    {prs_df['prs_z'].max():.4f}")

    prs_df.to_csv("./outputs/ms_prs_scores.csv")
    print(f"\nSaved -> ./outputs/ms_prs_scores.csv")
    return prs_df, snps

def plot_prs_distribution(prs_df):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor(DARK)

    # ── Plot 1: PRS distribution ──
    ax = axes[0]
    ax.set_facecolor(DARK)

    n, bins, patches = ax.hist(prs_df["prs_z"], bins=50,
                                color=BLUE, edgecolor=DARK, linewidth=0.3, alpha=0.85)

    # Colour top 10% red
    top10 = prs_df["prs_z"].quantile(0.90)
    for patch, left in zip(patches, bins[:-1]):
        if left >= top10:
            patch.set_facecolor(HIGHLIGHT)

    ax.axvline(top10, color=GOLD, linestyle="--", linewidth=1.2, label="Top 10%")
    ax.axvline(0, color=MUTED, linestyle=":", linewidth=0.8)

    ax.set_xlabel("Standardised PRS (Z-score)", color=GOLD, fontsize=10)
    ax.set_ylabel("Number of Individuals", color=GOLD, fontsize=10)
    ax.set_title("MS Polygenic Risk Score\nDistribution (n=1000)",
                 color=GOLD, fontsize=11, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)
    ax.legend(facecolor=DARK, edgecolor=GOLD, labelcolor="white", fontsize=9)

    # ── Plot 2: Top 20 SNPs by effect size ──
    ax2 = axes[1]
    ax2.set_facecolor(DARK)

    top_snps = prs_df  # reuse snps from outer scope via closure
    return fig, top10

def plot_top_snps(snps, ax):
    top20 = snps.nlargest(20, "beta").copy()
    top20["label"] = top20.index + "\n" + top20["gene"].fillna("").astype(str).str.split(",").str[0].str.strip()

    colors = [HIGHLIGHT if "HLA" in str(g) else BLUE for g in top20["gene"].fillna("")]

    bars = ax.barh(range(len(top20)), top20["beta"].values,
                   color=colors, edgecolor=GOLD, linewidth=0.5)
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20.index.values, fontsize=7, color="white")
    ax.set_xlabel("Effect Size log(OR)", color=GOLD, fontsize=10)
    ax.set_title("Top 20 MS Risk SNPs\nby Effect Size",
                 color=GOLD, fontsize=11, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)

    hla_patch  = mpatches.Patch(color=HIGHLIGHT, label="HLA region")
    other_patch = mpatches.Patch(color=BLUE,     label="Non-HLA")
    ax.legend(handles=[hla_patch, other_patch],
              facecolor=DARK, edgecolor=GOLD, labelcolor="white", fontsize=8)

if __name__ == "__main__":
    prs_df, snps = calculate_prs(
        "./outputs/simulated_genotypes.csv",
        "./outputs/ms_gwas_snps.csv"
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.patch.set_facecolor(DARK)
    for ax in axes:
        ax.set_facecolor(DARK)

    # Distribution plot
    top10 = prs_df["prs_z"].quantile(0.90)
    n, bins, patches = axes[0].hist(prs_df["prs_z"], bins=50,
                                     color=BLUE, edgecolor=DARK, linewidth=0.3, alpha=0.85)
    for patch, left in zip(patches, bins[:-1]):
        if left >= top10:
            patch.set_facecolor(HIGHLIGHT)
    axes[0].axvline(top10, color=GOLD, linestyle="--", linewidth=1.2, label="Top 10%")
    axes[0].axvline(0, color=MUTED, linestyle=":", linewidth=0.8)
    axes[0].set_xlabel("Standardised PRS (Z-score)", color=GOLD, fontsize=10)
    axes[0].set_ylabel("Number of Individuals", color=GOLD, fontsize=10)
    axes[0].set_title("MS Polygenic Risk Score\nDistribution (n=1000)",
                      color=GOLD, fontsize=11, fontweight="bold")
    axes[0].tick_params(colors="white")
    for spine in axes[0].spines.values():
        spine.set_edgecolor(MUTED)
    axes[0].legend(facecolor=DARK, edgecolor=GOLD, labelcolor="white", fontsize=9)

    # Top SNPs plot
    plot_top_snps(snps, axes[1])

    plt.tight_layout()
    plt.savefig("./outputs/ms_prs_results.png", dpi=180, bbox_inches="tight", facecolor=DARK)
    plt.close()
    print("Plot saved -> ./outputs/ms_prs_results.png")

    print(f"\nTop 10% PRS threshold (Z): {top10:.3f}")
    print(f"Individuals above threshold: {(prs_df['prs_z'] >= top10).sum()}")