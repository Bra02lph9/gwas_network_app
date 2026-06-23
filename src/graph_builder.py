from __future__ import annotations

from typing import Dict, Any

import pandas as pd
import networkx as nx


def build_variant_gene_disease_graph(
    df: pd.DataFrame,
) -> nx.Graph:
    """
    Build an undirected Variant → Gene → Disease biological network.

    Biological structure:
        SNP / Variant → Gene → Disease

    Node types:
        - variant
        - gene
        - disease

    Edge types:
        - mapped_to_gene
        - associated_with
    """

    if df is None or df.empty:
        return nx.Graph()

    required_columns = {"snp", "gene", "disease", "p_value"}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise ValueError(
            f"Missing columns for graph construction: {sorted(missing_columns)}"
        )

    G = nx.Graph()

    for _, row in df.iterrows():
        snp = str(row["snp"]).strip()
        gene = str(row["gene"]).strip().upper()
        disease = str(row["disease"]).strip()
        p_value = float(row["p_value"])

        if not snp or not gene or not disease:
            continue

        G.add_node(
            snp,
            node_type="variant",
            label=snp,
        )

        G.add_node(
            gene,
            node_type="gene",
            label=gene,
        )

        G.add_node(
            disease,
            node_type="disease",
            label=disease,
        )

        _add_or_update_edge(
            G=G,
            source=snp,
            target=gene,
            relation="mapped_to_gene",
            p_value=p_value,
        )

        _add_or_update_edge(
            G=G,
            source=gene,
            target=disease,
            relation="associated_with",
            p_value=p_value,
        )

    _update_node_attributes(G)

    return G


def _add_or_update_edge(
    G: nx.Graph,
    source: str,
    target: str,
    relation: str,
    p_value: float,
) -> None:
    """
    Add an edge or update its weight and best p-value.
    """

    if G.has_edge(source, target):
        G[source][target]["weight"] += 1

        current_best = G[source][target].get("best_p_value", p_value)

        if p_value < current_best:
            G[source][target]["best_p_value"] = p_value

    else:
        G.add_edge(
            source,
            target,
            relation=relation,
            weight=1,
            best_p_value=p_value,
        )


def _update_node_attributes(G: nx.Graph) -> None:
    """
    Add degree information to each node.
    """

    degree_dict = dict(G.degree())

    for node in G.nodes:
        G.nodes[node]["degree"] = degree_dict.get(node, 0)


def extract_subgraph_by_degree(
    G: nx.Graph,
    min_degree: int = 2,
) -> nx.Graph:
    """
    Keep nodes having at least min_degree connections.
    """

    if G is None or G.number_of_nodes() == 0:
        return nx.Graph()

    if min_degree <= 0:
        return G.copy()

    keep_nodes = [
        node for node, degree in G.degree()
        if degree >= min_degree
    ]

    subgraph = G.subgraph(keep_nodes).copy()
    _update_node_attributes(subgraph)

    return subgraph


def extract_top_nodes_by_degree(
    G: nx.Graph,
    top_n: int = 300,
) -> nx.Graph:
    """
    Keep only the top connected nodes based on degree.
    """

    if G is None or G.number_of_nodes() == 0:
        return nx.Graph()

    if top_n <= 0:
        raise ValueError("top_n must be greater than 0.")

    degree_dict = dict(G.degree())

    sorted_nodes = sorted(
        degree_dict,
        key=degree_dict.get,
        reverse=True,
    )

    selected_nodes = sorted_nodes[:top_n]

    subgraph = G.subgraph(selected_nodes).copy()
    _update_node_attributes(subgraph)

    return subgraph


def graph_statistics(
    G: nx.Graph,
) -> Dict[str, Any]:
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

    variants = 0
    genes = 0
    diseases = 0

    for _, data in G.nodes(data=True):
        node_type = data.get("node_type")

        if node_type == "variant":
            variants += 1
        elif node_type == "gene":
            genes += 1
        elif node_type == "disease":
            diseases += 1

    return {
        "nodes": G.number_of_nodes(),
        "edges": G.number_of_edges(),
        "variants": variants,
        "genes": genes,
        "diseases": diseases,
        "connected_components": nx.number_connected_components(G),
        "density": round(nx.density(G), 6),
    }
