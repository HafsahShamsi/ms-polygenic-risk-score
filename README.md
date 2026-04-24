# MS Polygenic Risk Score

**Author:** Hafsah Shamsi · [github.com/HafsahShamsi](https://github.com/HafsahShamsi)

A end-to-end polygenic risk score (PRS) pipeline for Multiple Sclerosis, built from public GWAS data.

---

## What this does

Polygenic risk scores aggregate the cumulative genetic burden of a complex disease across thousands of variants. For MS, this means summing the weighted contributions of 379 genome-wide significant SNPs (p < 5×10⁻⁸) identified across published GWAS studies.

**PRS = Σ (log(OR) × genotype dosage)** across all risk SNPs

The resulting score places each individual on a continuous risk spectrum — the top 10% of this distribution carries disproportionately higher genetic liability for MS.

---

## Pipeline

| Script | What it does |
|---|---|
| `fetch_ms_gwas.py` | Downloads GWAS Catalog associations, filters for MS, extracts genome-wide significant SNPs with effect sizes |
| `simulate_genotypes.py` | Simulates a cohort of 1000 individuals using population-level allele frequencies (binomial dosage model) |
| `calculate_prs.py` | Computes PRS via dot product of genotype matrix and log(OR) weights, standardises scores, visualises distribution and top SNPs |

---

## Key findings

- **379 SNPs** passed genome-wide significance threshold (p < 5×10⁻⁸)
- **HLA region dominance** — the majority of high-effect SNPs map to HLA-DQA1, HLA-DRB1, HLA-DRA, HLA-DRB9. This reflects the well-established role of MHC class II variation in MS susceptibility
- **Non-HLA hits** include CD58, IL2RA, CLEC16A, TNFRSF1A — immune regulation and T-cell signalling genes, consistent with MS as a T-cell mediated disease
- **Top 10% threshold:** Z-score ≥ 1.327

---

## Outputs

```
outputs/
├── ms_gwas_snps.csv           # 379 MS GWAS hits with effect sizes
├── simulated_genotypes.csv    # 1000 × 379 genotype dosage matrix
├── ms_prs_scores.csv          # PRS Z-scores + percentiles per individual
└── ms_prs_results.png         # Distribution + top SNPs visualisation
```

---

## Limitations & notes

- Genotypes are **simulated** from population allele frequencies — this is a pipeline demonstration, not a clinical tool
- Effect sizes are drawn from GWAS Catalog 2022 release; some OR values reflect heterogeneous study populations
- PRS does not capture gene-gene interactions, rare variants, or environmental factors

---

## Dependencies

pip install requests pandas numpy matplotlib seaborn scipy

---

## Context

This project is part of a broader computational neuroimmunology pipeline:

- [`neuro-deg-scanner`](https://github.com/HafsahShamsi/neuro-deg-scanner) — reusable DEG pipeline for GEO datasets
- [`ms-pathway-enrichment`](https://github.com/HafsahShamsi/ms-pathway-enrichment) — KEGG/GO enrichment on MS DEGs
- `ms-polygenic-risk-score` — this repo
- `drug-target-analysis` — coming next
- `drug-target-analysis` — coming next
