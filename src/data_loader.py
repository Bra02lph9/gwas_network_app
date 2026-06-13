from __future__ import annotations

from typing import Dict

import pandas as pd


REQUIRED_COLUMNS = [
    "snp",
    "gene",
    "disease",
    "p_value",
]


def read_uploaded_file(uploaded_file) -> pd.DataFrame:
    """
    Read uploaded CSV/TSV/TXT file from Streamlit.
    """

    filename = uploaded_file.name.lower()

    if filename.endswith(".csv"):
        df = pd.read_csv(uploaded_file)

    elif filename.endswith(".tsv"):
        df = pd.read_csv(uploaded_file, sep="\t")

    elif filename.endswith(".txt"):
        try:
            df = pd.read_csv(uploaded_file, sep="\t")
        except Exception:
            uploaded_file.seek(0)
            df = pd.read_csv(uploaded_file)

    else:
        raise ValueError(
            "Unsupported file format. "
            "Supported formats: CSV, TSV, TXT."
        )

    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
    )

    return df


def validate_network_input(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate and clean network input table.
    """

    missing_columns = [
        column
        for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}\n"
            f"Required columns: {REQUIRED_COLUMNS}"
        )

    df = df.copy()

    df["snp"] = (
        df["snp"]
        .astype(str)
        .str.strip()
    )

    df["gene"] = (
        df["gene"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    df["disease"] = (
        df["disease"]
        .astype(str)
        .str.strip()
    )

    df["p_value"] = pd.to_numeric(
        df["p_value"],
        errors="coerce",
    )

    df = df.dropna(
        subset=[
            "snp",
            "gene",
            "disease",
            "p_value",
        ]
    )

    df = df[
        df["snp"].str.startswith("rs", na=False)
    ]

    df = df[
        df["p_value"] > 0
    ]

    df = df.drop_duplicates()

    return df.reset_index(drop=True)


def filter_by_pvalue(
    df: pd.DataFrame,
    pvalue_threshold: float,
) -> pd.DataFrame:
    """
    Keep only significant associations.
    """

    if df.empty:
        return df

    filtered_df = df[
        df["p_value"] <= pvalue_threshold
    ]

    return filtered_df.reset_index(drop=True)


def keep_top_significant_rows(
    df: pd.DataFrame,
    top_n: int = 500,
) -> pd.DataFrame:
    """
    Keep top N most significant associations.
    """

    if df.empty:
        return df

    return (
        df.sort_values(
            "p_value",
            ascending=True,
        )
        .head(top_n)
        .reset_index(drop=True)
    )


def keep_genes_with_min_snps(
    df: pd.DataFrame,
    min_snps_per_gene: int = 2,
) -> pd.DataFrame:
    """
    Remove weak genes supported by too few SNPs.
    """

    if df.empty:
        return df

    if min_snps_per_gene <= 1:
        return df

    gene_counts = (
        df.groupby("gene")["snp"]
        .nunique()
    )

    valid_genes = gene_counts[
        gene_counts >= min_snps_per_gene
    ].index

    return (
        df[df["gene"].isin(valid_genes)]
        .reset_index(drop=True)
    )


def get_dataset_summary(
    df: pd.DataFrame,
) -> Dict[str, int]:
    """
    Return dataset statistics.
    """

    if df.empty:
        return {
            "rows": 0,
            "unique_snps": 0,
            "unique_genes": 0,
            "unique_diseases": 0,
        }

    return {
        "rows": len(df),
        "unique_snps": df["snp"].nunique(),
        "unique_genes": df["gene"].nunique(),
        "unique_diseases": df["disease"].nunique(),
    }
