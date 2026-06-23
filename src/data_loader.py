from __future__ import annotations

from typing import Dict, Any

import pandas as pd


REQUIRED_COLUMNS = ["snp", "gene", "disease", "p_value"]


def read_uploaded_file(uploaded_file) -> pd.DataFrame:
    """
    Read uploaded CSV, TSV, or TXT file from Streamlit.
    """

    if uploaded_file is None:
        raise ValueError("No file was uploaded.")

    filename = uploaded_file.name.lower()

    try:
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
                "Unsupported file format. Supported formats: CSV, TSV, TXT."
            )

    except Exception as exc:
        raise ValueError(f"Could not read uploaded file: {exc}") from exc

    df.columns = (
        df.columns
        .astype(str)
        .str.strip()
        .str.lower()
    )

    return df


def validate_network_input(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate and clean GWAS network input table.

    Required columns:
    - snp
    - gene
    - disease
    - p_value
    """

    if df is None or df.empty:
        raise ValueError("Input table is empty.")

    missing_columns = [
        column for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}. "
            f"Required columns are: {REQUIRED_COLUMNS}"
        )

    df = df.copy()

    df["snp"] = df["snp"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip().str.upper()
    df["disease"] = df["disease"].astype(str).str.strip()

    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    df = df.replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})

    df = df.dropna(
        subset=["snp", "gene", "disease", "p_value"]
    )

    df = df[df["snp"].str.lower().str.startswith("rs", na=False)]

    df = df[df["p_value"] > 0]
    df = df[df["p_value"] <= 1]

    df = df.drop_duplicates(
        subset=["snp", "gene", "disease", "p_value"]
    )

    df = df.sort_values("p_value", ascending=True)

    return df.reset_index(drop=True)


def filter_by_pvalue(
    df: pd.DataFrame,
    pvalue_threshold: float,
) -> pd.DataFrame:
    """
    Keep associations with p_value <= threshold.
    """

    if df.empty:
        return df

    if pvalue_threshold <= 0:
        raise ValueError("p-value threshold must be greater than 0.")

    return (
        df[df["p_value"] <= pvalue_threshold]
        .reset_index(drop=True)
    )


def keep_top_significant_rows(
    df: pd.DataFrame,
    top_n: int = 500,
) -> pd.DataFrame:
    """
    Keep top N most significant GWAS associations.
    """

    if df.empty:
        return df

    if top_n <= 0:
        raise ValueError("top_n must be greater than 0.")

    return (
        df.sort_values("p_value", ascending=True)
        .head(top_n)
        .reset_index(drop=True)
    )


def keep_genes_with_min_snps(
    df: pd.DataFrame,
    min_snps_per_gene: int = 2,
) -> pd.DataFrame:
    """
    Keep genes supported by at least N unique SNPs.
    """

    if df.empty:
        return df

    if min_snps_per_gene <= 1:
        return df.reset_index(drop=True)

    gene_snp_counts = df.groupby("gene")["snp"].nunique()

    valid_genes = gene_snp_counts[
        gene_snp_counts >= min_snps_per_gene
    ].index

    return (
        df[df["gene"].isin(valid_genes)]
        .reset_index(drop=True)
    )


def get_dataset_summary(df: pd.DataFrame) -> Dict[str, Any]:
    """
    Return useful dataset statistics for Streamlit display.
    """

    if df.empty:
        return {
            "rows": 0,
            "unique_snps": 0,
            "unique_genes": 0,
            "unique_diseases": 0,
            "min_p_value": None,
            "max_p_value": None,
        }

    return {
        "rows": int(len(df)),
        "unique_snps": int(df["snp"].nunique()),
        "unique_genes": int(df["gene"].nunique()),
        "unique_diseases": int(df["disease"].nunique()),
        "min_p_value": float(df["p_value"].min()),
        "max_p_value": float(df["p_value"].max()),
    }
