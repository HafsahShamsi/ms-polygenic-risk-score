"""
Step 2: Simulate genotype data for PRS demonstration
=====================================================
Author: Hafsah Shamsi
GitHub: github.com/HafsahShamsi

Since real patient genotype data isn't publicly available in a simple format,
we simulate a cohort using population-level allele frequencies from the GWAS hits.
This is standard practice for PRS pipeline demonstration.

Each individual gets 0, 1, or 2 copies of the risk allele per SNP,
drawn from a binomial distribution based on the risk allele frequency.
"""

import pandas as pd
import numpy as np
import os

os.makedirs("./outputs", exist_ok=True)
np.random.seed(42)

def simulate_genotypes(snp_file, n_individuals=1000):
    snps = pd.read_csv(snp_file)
    print(f"Loaded {len(snps)} SNPs")

    # Drop SNPs with no allele frequency info
    snps["risk_allele_freq"] = pd.to_numeric(snps["risk_allele_freq"], errors="coerce")
    
    # Fill missing frequencies with population average ~0.3
    snps["risk_allele_freq"] = snps["risk_allele_freq"].fillna(0.3)
    snps["risk_allele_freq"] = snps["risk_allele_freq"].clip(0.01, 0.99)

    print(f"Simulating genotypes for {n_individuals} individuals x {len(snps)} SNPs...")

    # Each person gets 0/1/2 copies of risk allele (diploid)
    # P(0) = (1-f)^2, P(1) = 2f(1-f), P(2) = f^2
    genotypes = np.array([
        np.random.binomial(2, freq, n_individuals)
        for freq in snps["risk_allele_freq"]
    ]).T  # shape: (n_individuals, n_snps)

    geno_df = pd.DataFrame(
        genotypes,
        columns=snps["rsid"].values
    )
    geno_df.index.name = "individual_id"
    geno_df.index = [f"IND_{i:04d}" for i in range(n_individuals)]

    geno_df.to_csv("./outputs/simulated_genotypes.csv")
    print(f"Saved -> ./outputs/simulated_genotypes.csv")
    print(f"Shape: {geno_df.shape}")
    return geno_df, snps

if __name__ == "__main__":
    geno_df, snps = simulate_genotypes("./outputs/ms_gwas_snps.csv")
    print("\nSample (first 5 individuals, first 5 SNPs):")
    print(geno_df.iloc[:5, :5])