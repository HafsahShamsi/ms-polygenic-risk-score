"""
Step 1: Fetch MS GWAS hits from GWAS Catalog
=============================================
Author: Hafsah Shamsi
GitHub: github.com/HafsahShamsi
"""

import requests
import pandas as pd
import io
import os

os.makedirs("./outputs", exist_ok=True)

def fetch_ms_gwas():
    print("Downloading GWAS Catalog associations TSV...")
    url = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2022/05/23/gwas-catalog-associations_ontology-annotated.tsv"
    r = requests.get(url, stream=True)
    print(f"  Status: {r.status_code}")
    if r.status_code != 200:
        print(f"  Failed: {r.text[:200]}")
        return pd.DataFrame()
    print("  Download complete, parsing...")
    content = r.content.decode("utf-8")
    df = pd.read_csv(io.StringIO(content), sep="\t", low_memory=False)
    print(f"  Total rows: {len(df)}")
    print(f"  Columns: {df.columns.tolist()[:10]}")
    return df

def parse_ms_associations(df):
    if df.empty:
        return pd.DataFrame()
    ms = df[df["DISEASE/TRAIT"].str.contains("multiple sclerosis", case=False, na=False)].copy()
    print(f"\nMS associations: {len(ms)}")
    ms["P-VALUE"] = pd.to_numeric(ms["P-VALUE"], errors="coerce")
    ms = ms[ms["P-VALUE"] < 5e-8]
    print(f"Genome-wide significant (p<5e-8): {len(ms)}")
    cols = {
        "SNPS": "rsid",
        "RISK ALLELE FREQUENCY": "risk_allele_freq",
        "P-VALUE": "pvalue",
        "OR or BETA": "effect_size",
        "STRONGEST SNP-RISK ALLELE": "snp_risk_allele",
        "MAPPED_GENE": "gene",
        "CHR_ID": "chromosome",
        "CHR_POS": "position"
    }
    available = {k: v for k, v in cols.items() if k in ms.columns}
    result = ms[list(available.keys())].rename(columns=available)
    result["effect_size"] = pd.to_numeric(result["effect_size"], errors="coerce")
    result = result.dropna(subset=["rsid", "effect_size"])
    result = result[result.rsid.str.startswith("rs", na=False)]
    result = result.drop_duplicates(subset="rsid")
    result = result.sort_values("pvalue")
    return result

if __name__ == "__main__":
    df = fetch_ms_gwas()
    result = parse_ms_associations(df)
    if result.empty:
        print("No SNPs parsed.")
    else:
        print(f"\nFinal SNP count: {len(result)}")
        print(result.head(10).to_string(index=False))
        result.to_csv("./outputs/ms_gwas_snps.csv", index=False)
        print(f"\nSaved -> ./outputs/ms_gwas_snps.csv")
