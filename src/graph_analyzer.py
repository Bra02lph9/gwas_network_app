from __future__ import annotations

from typing import Any

import pandas as pd
import networkx as nx


def get_graph_summary(G: nx.Graph) -> dict[str, int]:
    """
    Return global graph statistics.
    """

    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "variants": sum(
            1 for _, data in G.nodes(data=True)
            if data.get("node_type") == "variant"
        ),
        "genes": sum(
            1 for _, data in G.nodes(data=True)
            if data.get("node_type") == "gene"
        ),
        "diseases": sum(
            1 for _, data in G.nodes(data=True)
            if data.get("node_type") == "disease"
        ),
    }


def calculate_node_metrics(G: nx.Graph) -> pd.DataFrame:
    """
    Calculate centrality metrics for all nodes.
    """

    if G.number_of_nodes() == 0:
        return pd.DataFrame()

    degree = dict(G.degree())
    degree_centrality = nx.degree_centrality(G)

    if G.number_of_nodes() > 2:
        betweenness = nx.betweenness_centrality(G, normalized=True)
    else:
        betweenness = {node: 0.0 for node in G.nodes()}

    try:
        eigenvector = nx.eigenvector_centrality(
            G,
            max_iter=1000,
            tol=1e-06,
        )
    except nx.NetworkXException:
        eigenvector = {node: 0.0 for node in G.nodes()}

    rows: list[dict[str, Any]] = []

    for node, data in G.nodes(data=True):
        node_type = data.get("node_type", "unknown")

        rows.append(
            {
                "node": node,
                "node_type": node_type,
                "degree": degree.get(node, 0),
                "degree_centrality": degree_centrality.get(node, 0.0),
                "betweenness_centrality": betweenness.get(node, 0.0),
                "eigenvector_centrality": eigenvector.get(node, 0.0),
            }
        )

    metrics_df = pd.DataFrame(rows)

    if metrics_df.empty:
        return metrics_df

    metrics_df["hub_score"] = (
        metrics_df["degree_centrality"]
        + metrics_df["betweenness_centrality"]
        + metrics_df["eigenvector_centrality"]
    )

    return metrics_df.sort_values(
        "hub_score",
        ascending=False,
    ).reset_index(drop=True)


def get_top_hub_genes(
    G: nx.Graph,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Return top hub genes.
    """

    metrics_df = calculate_node_metrics(G)

    if metrics_df.empty:
        return pd.DataFrame(
            columns=[
                "gene",
                "degree",
                "degree_centrality",
                "betweenness_centrality",
                "eigenvector_centrality",
                "hub_score",
            ]
        )

    genes_df = metrics_df[
        metrics_df["node_type"] == "gene"
    ].copy()

    genes_df = genes_df.rename(
        columns={
            "node": "gene",
        }
    )

    return genes_df[
        [
            "gene",
            "degree",
            "degree_centrality",
            "betweenness_centrality",
            "eigenvector_centrality",
            "hub_score",
        ]
    ].head(top_n).reset_index(drop=True)


def get_top_hub_variants(
    G: nx.Graph,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Return top hub variants/SNPs.
    """

    metrics_df = calculate_node_metrics(G)

    if metrics_df.empty:
        return pd.DataFrame(
            columns=[
                "snp",
                "degree",
                "degree_centrality",
                "betweenness_centrality",
                "eigenvector_centrality",
                "hub_score",
            ]
        )

    variants_df = metrics_df[
        metrics_df["node_type"] == "variant"
    ].copy()

    variants_df = variants_df.rename(
        columns={
            "node": "snp",
        }
    )

    return variants_df[
        [
            "snp",
            "degree",
            "degree_centrality",
            "betweenness_centrality",
            "eigenvector_centrality",
            "hub_score",
        ]
    ].head(top_n).reset_index(drop=True)


def get_top_connected_diseases(
    G: nx.Graph,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Return top connected diseases/traits.
    """

    metrics_df = calculate_node_metrics(G)

    if metrics_df.empty:
        return pd.DataFrame(
            columns=[
                "disease",
                "degree",
                "degree_centrality",
                "betweenness_centrality",
                "eigenvector_centrality",
                "hub_score",
            ]
        )

    diseases_df = metrics_df[
        metrics_df["node_type"] == "disease"
    ].copy()

    diseases_df = diseases_df.rename(
        columns={
            "node": "disease",
        }
    )

    return diseases_df[
        [
            "disease",
            "degree",
            "degree_centrality",
            "betweenness_centrality",
            "eigenvector_centrality",
            "hub_score",
        ]
    ].head(top_n).reset_index(drop=True)


def build_network_report(
    G: nx.Graph,
) -> str:
    """
    Create a simple text report for the network.
    """

    summary = get_graph_summary(G)
    top_genes = get_top_hub_genes(G, top_n=5)

    top_gene_names = (
        ", ".join(top_genes["gene"].tolist())
        if not top_genes.empty
        else "None"
    )

    return (
        f"Network contains {summary['nodes']} nodes and {summary['edges']} edges. "
        f"It includes {summary['variants']} variants, {summary['genes']} genes, "
        f"and {summary['diseases']} diseases/traits. "
        f"The top hub genes are: {top_gene_names}."
    )
