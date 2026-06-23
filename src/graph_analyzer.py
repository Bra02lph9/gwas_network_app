from __future__ import annotations

from typing import Any

import pandas as pd
import networkx as nx


def get_graph_summary(G: nx.Graph) -> dict[str, Any]:
    """
    Return global graph statistics.
    """

    if G is None or G.number_of_nodes() == 0:
        return {
            "nodes": 0,
            "edges": 0,
            "variants": 0,
            "genes": 0,
            "diseases": 0,
            "connected_components": 0,
            "density": 0.0,
        }

    node_types = nx.get_node_attributes(G, "node_type")

    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "variants": sum(1 for t in node_types.values() if t == "variant"),
        "genes": sum(1 for t in node_types.values() if t == "gene"),
        "diseases": sum(1 for t in node_types.values() if t == "disease"),
        "connected_components": nx.number_connected_components(G),
        "density": round(nx.density(G), 6),
    }


def calculate_node_metrics(G: nx.Graph) -> pd.DataFrame:
    """
    Calculate centrality metrics for all nodes.
    """

    if G is None or G.number_of_nodes() == 0:
        return pd.DataFrame()

    degree = dict(G.degree())
    degree_centrality = nx.degree_centrality(G)
    closeness = nx.closeness_centrality(G)

    if G.number_of_nodes() > 2:
        betweenness = nx.betweenness_centrality(G, normalized=True)
    else:
        betweenness = {node: 0.0 for node in G.nodes()}

    try:
        eigenvector = nx.eigenvector_centrality(
            G,
            max_iter=1000,
            tol=1e-6,
        )
    except nx.NetworkXException:
        eigenvector = {node: 0.0 for node in G.nodes()}

    rows: list[dict[str, Any]] = []

    for node, data in G.nodes(data=True):
        rows.append(
            {
                "node": node,
                "node_type": data.get("node_type", "unknown"),
                "degree": degree.get(node, 0),
                "degree_centrality": round(degree_centrality.get(node, 0.0), 6),
                "betweenness_centrality": round(betweenness.get(node, 0.0), 6),
                "closeness_centrality": round(closeness.get(node, 0.0), 6),
                "eigenvector_centrality": round(eigenvector.get(node, 0.0), 6),
            }
        )

    metrics_df = pd.DataFrame(rows)

    if metrics_df.empty:
        return metrics_df

    metrics_df["hub_score"] = (
        metrics_df["degree_centrality"]
        + metrics_df["betweenness_centrality"]
        + metrics_df["closeness_centrality"]
        + metrics_df["eigenvector_centrality"]
    ).round(6)

    return (
        metrics_df
        .sort_values("hub_score", ascending=False)
        .reset_index(drop=True)
    )


def get_top_hub_genes(
    G: nx.Graph,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Return top hub genes.
    """

    return _get_top_nodes_by_type(
        G=G,
        node_type="gene",
        output_name="gene",
        top_n=top_n,
    )


def get_top_hub_variants(
    G: nx.Graph,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Return top hub variants/SNPs.
    """

    return _get_top_nodes_by_type(
        G=G,
        node_type="variant",
        output_name="snp",
        top_n=top_n,
    )


def get_top_connected_diseases(
    G: nx.Graph,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Return top connected diseases/traits.
    """

    return _get_top_nodes_by_type(
        G=G,
        node_type="disease",
        output_name="disease",
        top_n=top_n,
    )


def _get_top_nodes_by_type(
    G: nx.Graph,
    node_type: str,
    output_name: str,
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Internal helper to extract top nodes by biological type.
    """

    columns = [
        output_name,
        "degree",
        "degree_centrality",
        "betweenness_centrality",
        "closeness_centrality",
        "eigenvector_centrality",
        "hub_score",
    ]

    metrics_df = calculate_node_metrics(G)

    if metrics_df.empty:
        return pd.DataFrame(columns=columns)

    filtered_df = metrics_df[
        metrics_df["node_type"] == node_type
    ].copy()

    if filtered_df.empty:
        return pd.DataFrame(columns=columns)

    filtered_df = filtered_df.rename(columns={"node": output_name})

    return (
        filtered_df[columns]
        .head(top_n)
        .reset_index(drop=True)
    )


def build_network_report(G: nx.Graph) -> str:
    """
    Create a cautious biological text report for the network.
    """

    summary = get_graph_summary(G)
    top_genes = get_top_hub_genes(G, top_n=5)

    if summary["nodes"] == 0:
        return "No network could be generated from the provided data."

    top_gene_names = (
        ", ".join(top_genes["gene"].tolist())
        if not top_genes.empty
        else "None"
    )

    return (
        f"The network contains {summary['nodes']} nodes and "
        f"{summary['edges']} edges. It includes "
        f"{summary['variants']} variants, {summary['genes']} genes, "
        f"and {summary['diseases']} diseases or traits. "
        f"The graph contains {summary['connected_components']} connected "
        f"component(s), with a density of {summary['density']}. "
        f"The top hub genes are: {top_gene_names}. "
        f"These hub genes may represent central elements in the GWAS network, "
        f"but they should not be interpreted as causal genes without additional "
        f"functional validation."
    )
