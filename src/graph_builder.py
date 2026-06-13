from __future__ import annotations

import pandas as pd
import networkx as nx


def build_variant_gene_disease_graph(
    df: pd.DataFrame,
) -> nx.Graph:
    """
    Build Variant → Gene → Disease graph.

    Node types:
        - variant
        - gene
        - disease

    Edge types:
        - mapped_to_gene
        - associated_with
    """

    G = nx.Graph()

    for _, row in df.iterrows():

        snp = str(row["snp"]).strip()
        gene = str(row["gene"]).strip()
        disease = str(row["disease"]).strip()

        p_value = float(row["p_value"])

        if not G.has_node(snp):
            G.add_node(
                snp,
                node_type="variant",
                label=snp,
            )

        if not G.has_node(gene):
            G.add_node(
                gene,
                node_type="gene",
                label=gene,
            )

        if not G.has_node(disease):
            G.add_node(
                disease,
                node_type="disease",
                label=disease,
            )

        # SNP -> Gene

        if G.has_edge(snp, gene):
            G[snp][gene]["weight"] += 1

            if p_value < G[snp][gene]["best_p_value"]:
                G[snp][gene]["best_p_value"] = p_value

        else:
            G.add_edge(
                snp,
                gene,
                relation="mapped_to_gene",
                weight=1,
                best_p_value=p_value,
            )

        # Gene -> Disease

        if G.has_edge(gene, disease):
            G[gene][disease]["weight"] += 1

            if p_value < G[gene][disease]["best_p_value"]:
                G[gene][disease]["best_p_value"] = p_value

        else:
            G.add_edge(
                gene,
                disease,
                relation="associated_with",
                weight=1,
                best_p_value=p_value,
            )

    return G


def extract_subgraph_by_degree(
    G: nx.Graph,
    min_degree: int = 2,
) -> nx.Graph:
    """
    Keep nodes having at least min_degree connections.
    """

    keep_nodes = [
        node
        for node, degree in G.degree()
        if degree >= min_degree
    ]

    return G.subgraph(keep_nodes).copy()


def extract_top_nodes_by_degree(
    G: nx.Graph,
    top_n: int = 300,
) -> nx.Graph:
    """
    Keep only top connected nodes.
    """

    degree_dict = dict(G.degree())

    sorted_nodes = sorted(
        degree_dict,
        key=degree_dict.get,
        reverse=True,
    )

    selected_nodes = sorted_nodes[:top_n]

    return G.subgraph(selected_nodes).copy()


def graph_statistics(
    G: nx.Graph,
) -> dict:
    """
    Quick graph statistics.
    """

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
    }
