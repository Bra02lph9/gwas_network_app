from __future__ import annotations

import math

import networkx as nx
import plotly.graph_objects as go


NODE_COLORS = {
    "variant": "#4E79A7",
    "gene": "#59A14F",
    "disease": "#E15759",
    "unknown": "#BAB0AC",
}

NODE_SIZES = {
    "variant": 5,
    "gene": 10,
    "disease": 16,
    "unknown": 6,
}


def _safe_log_pvalue(p_value: float | None) -> float:
    """
    Convert p-value to -log10(p-value).
    """

    if p_value is None:
        return 0.0

    try:
        p_value = float(p_value)
    except (TypeError, ValueError):
        return 0.0

    if p_value <= 0:
        return 0.0

    return -math.log10(p_value)


def _get_edge_width(best_p_value: float | None) -> float:
    """
    Make stronger associations visually thicker.
    """

    log_p = _safe_log_pvalue(best_p_value)

    if log_p >= 50:
        return 5.0
    if log_p >= 20:
        return 4.0
    if log_p >= 10:
        return 3.0
    if log_p >= 7.3:
        return 2.0

    return 1.0


def create_3d_network_figure(G: nx.Graph) -> go.Figure:
    """
    Create an interactive 3D Plotly graph from a NetworkX graph.
    """

    if G.number_of_nodes() == 0:
        fig = go.Figure()
        fig.update_layout(
            title="Empty network",
            height=750,
        )
        return fig

    k_value = 0.8

    if G.number_of_nodes() > 300:
        k_value = 1.2

    pos = nx.spring_layout(
        G,
        dim=3,
        seed=42,
        k=k_value,
        iterations=80,
    )

    edge_traces = []

    for source, target, data in G.edges(data=True):
        x0, y0, z0 = pos[source]
        x1, y1, z1 = pos[target]

        best_p_value = data.get("best_p_value")
        relation = data.get("relation", "unknown")
        weight = data.get("weight", 1)

        edge_width = _get_edge_width(best_p_value)

        edge_trace = go.Scatter3d(
            x=[x0, x1, None],
            y=[y0, y1, None],
            z=[z0, z1, None],
            mode="lines",
            line=dict(
                width=edge_width,
                color="rgba(160,160,160,0.35)",
            ),
            hoverinfo="text",
            text=(
                f"Relation: {relation}<br>"
                f"Weight: {weight}<br>"
                f"Best p-value: {best_p_value}"
            ),
            showlegend=False,
        )

        edge_traces.append(edge_trace)

    node_x = []
    node_y = []
    node_z = []
    node_text = []
    node_color = []
    node_size = []
    node_symbol = []

    degrees = dict(G.degree())

    for node, data in G.nodes(data=True):
        x, y, z = pos[node]

        node_type = data.get("node_type", "unknown")
        degree = degrees.get(node, 0)

        node_x.append(x)
        node_y.append(y)
        node_z.append(z)

        node_color.append(
            NODE_COLORS.get(node_type, NODE_COLORS["unknown"])
        )

        size = NODE_SIZES.get(node_type, NODE_SIZES["unknown"]) + degree * 0.6
        node_size.append(min(size, 35))

        if node_type == "variant":
            symbol = "circle"
        elif node_type == "gene":
            symbol = "diamond"
        elif node_type == "disease":
            symbol = "square"
        else:
            symbol = "circle"

        node_symbol.append(symbol)

        node_text.append(
            f"<b>{node}</b><br>"
            f"Type: {node_type}<br>"
            f"Degree: {degree}"
        )

    node_trace = go.Scatter3d(
        x=node_x,
        y=node_y,
        z=node_z,
        mode="markers",
        text=node_text,
        hoverinfo="text",
        marker=dict(
            size=node_size,
            color=node_color,
            symbol=node_symbol,
            opacity=0.92,
            line=dict(
                width=0.7,
                color="white",
            ),
        ),
        showlegend=False,
    )

    fig = go.Figure(
        data=edge_traces + [node_trace]
    )

    fig.update_layout(
        title="3D Variant–Gene–Disease Network",
        height=780,
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=45,
        ),
        showlegend=False,
        scene=dict(
            xaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                title="",
            ),
            yaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                title="",
            ),
            zaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                title="",
            ),
        ),
    )

    return fig
