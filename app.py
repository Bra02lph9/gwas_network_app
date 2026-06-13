from __future__ import annotations

import pandas as pd
import streamlit as st

from src.config import create_directories
from src.gwas_fetcher import fetch_clean_gwas_network_data
from src.data_loader import (
    read_uploaded_file,
    validate_network_input,
    filter_by_pvalue,
)
from src.graph_builder import build_variant_gene_disease_graph
from src.graph_analyzer import get_graph_summary, get_top_hub_genes
from src.graph_3d_visualizer import create_3d_network_figure


st.set_page_config(
    page_title="GWAS 3D Network Explorer",
    layout="wide",
)

create_directories()

st.title("GWAS 3D Variant–Gene–Disease Network Explorer")
st.caption(
    "Interactive 3D visualization of variant–gene–disease networks from GWAS data."
)

st.sidebar.header("Data source")

data_source = st.sidebar.radio(
    "Choose data source",
    ["Upload file", "GWAS Catalog API"],
)

uploaded_file = None

if data_source == "Upload file":
    uploaded_file = st.sidebar.file_uploader(
        "Upload CSV or TSV file",
        type=["csv", "tsv", "txt"],
    )

st.sidebar.header("Filtering settings")

trait = st.sidebar.text_input(
    "Disease / trait",
    value="coronary artery disease",
    disabled=data_source == "Upload file",
)

pvalue_threshold = st.sidebar.number_input(
    "P-value threshold",
    value=5e-8,
    format="%.1e",
)

min_snps_per_gene = st.sidebar.slider(
    "Minimum SNPs per gene",
    min_value=1,
    max_value=20,
    value=1,
)

top_n_rows = st.sidebar.slider(
    "Top significant rows to keep",
    min_value=50,
    max_value=5000,
    value=500,
    step=50,
)

max_pages = st.sidebar.slider(
    "Max API pages",
    min_value=1,
    max_value=20,
    value=5,
    disabled=data_source == "Upload file",
)

max_rows_to_visualize = st.sidebar.slider(
    "Max rows to visualize",
    min_value=50,
    max_value=1000,
    value=300,
    step=50,
)


def keep_genes_with_min_snps(
    df: pd.DataFrame,
    min_snps: int,
) -> pd.DataFrame:
    if df.empty or min_snps <= 1:
        return df

    snp_counts = df.groupby("gene")["snp"].nunique()
    valid_genes = snp_counts[snp_counts >= min_snps].index

    return df[df["gene"].isin(valid_genes)].reset_index(drop=True)


if st.sidebar.button("Build 3D network"):

    try:
        if data_source == "Upload file":
            if uploaded_file is None:
                st.warning("Please upload a CSV or TSV file.")
                st.stop()

            raw_df = read_uploaded_file(uploaded_file)
            clean_df = validate_network_input(raw_df)
            filtered_df = filter_by_pvalue(clean_df, pvalue_threshold)

        else:
            filtered_df = fetch_clean_gwas_network_data(
                trait=trait,
                max_pages=max_pages,
                pvalue_threshold=pvalue_threshold,
                use_demo_if_empty=False,
            )

    except Exception as error:
        st.error(f"Data loading error: {error}")
        st.stop()

    if filtered_df.empty:
        st.warning("No data found after p-value filtering.")
        st.stop()

    filtered_df = (
        filtered_df
        .sort_values("p_value", ascending=True)
        .head(top_n_rows)
        .reset_index(drop=True)
    )

    filtered_df = keep_genes_with_min_snps(
        filtered_df,
        min_snps=min_snps_per_gene,
    )

    if filtered_df.empty:
        st.warning("No data left after the minimum SNPs per gene filter.")
        st.stop()

    filtered_df = filtered_df.head(max_rows_to_visualize).reset_index(drop=True)

    G = build_variant_gene_disease_graph(filtered_df)
    summary = get_graph_summary(G)

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Nodes", summary["nodes"])
    col2.metric("Edges", summary["edges"])
    col3.metric("Variants", summary["variants"])
    col4.metric("Genes", summary["genes"])

    st.subheader("Interactive 3D network")
    fig = create_3d_network_figure(G)
    st.plotly_chart(fig, use_container_width=True)

    st.subheader("Top hub genes")
    top_hub_genes = get_top_hub_genes(G)
    st.dataframe(top_hub_genes, use_container_width=True)

    st.subheader("Filtered associations used in the graph")
    st.dataframe(filtered_df, use_container_width=True)

    csv_data = filtered_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="Download filtered associations as CSV",
        data=csv_data,
        file_name="filtered_gwas_network_data.csv",
        mime="text/csv",
    )
