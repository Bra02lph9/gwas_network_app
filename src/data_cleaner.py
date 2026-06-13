import pandas as pd


def clean_gwas_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep useful GWAS fields and normalize column names.
    """
    if df.empty:
        return pd.DataFrame()

    clean_rows = []

    for _, row in df.iterrows():
        snp = row.get("riskAlleleName")
        p_value = row.get("pvalue")
        mapped_gene = row.get("mappedGene")
        trait = row.get("efoTraits")

        if not snp or not mapped_gene:
            continue

        clean_rows.append({
            "snp": snp,
            "gene": mapped_gene,
            "disease": str(trait),
            "p_value": p_value,
        })

    return pd.DataFrame(clean_rows)


def filter_significant(df: pd.DataFrame, threshold: float = 5e-8) -> pd.DataFrame:
    df = df.copy()
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    return df[df["p_value"] < threshold]
