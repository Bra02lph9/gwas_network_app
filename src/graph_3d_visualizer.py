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
    "variant": 6,
    "gene": 11,
    "disease": 16,
    "unknown": 6,
}

NODE_SYMBOLS = {
    "variant": "circle",
    "gene": "diamond",
    "disease": "square",
    "unknown": "circle",
}


def _safe_log_pvalue(p_value: float | None) -> float:
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


def _format_pvalue(p_value: float | None) -> str:
    try:
        return f"{float(p_value):.2e}"
    except (TypeError, ValueError):
        return "NA"


def create_3d_network_figure(G: nx.Graph) -> go.Figure:
    """
    Create an interactive 3D Plotly graph from a NetworkX graph.
    """

    if G is None or G.number_of_nodes() == 0:
        fig = go.Figure()
        fig.update_layout(
            title="Empty network",
            height=750,
        )
        return fig

    k_value = 1.2 if G.number_of_nodes() > 300 else 0.8

    pos = nx.spring_layout(
        G,
        dim=3,
        seed=42,
        k=k_value,
        iterations=100,
        weight="weight",
    )

    edge_traces = []

    for source, target, data in G.edges(data=True):
        x0, y0, z0 = pos[source]
        x1, y1, z1 = pos[target]

        best_p_value = data.get("best_p_value")
        relation = data.get("relation", "unknown")
        weight = data.get("weight", 1)

        edge_trace = go.Scatter3d(
            x=[x0, x1, None],
            y=[y0, y1, None],
            z=[z0, z1, None],
            mode="lines",
            line=dict(
                width=_get_edge_width(best_p_value),
                color="rgba(150,150,150,0.35)",
            ),
            hoverinfo="text",
            text=(
                f"<b>Edge</b><br>"
                f"Source: {source}<br>"
                f"Target: {target}<br>"
                f"Relation: {relation}<br>"
                f"Weight: {weight}<br>"
                f"Best p-value: {_format_pvalue(best_p_value)}"
            ),
            showlegend=False,
        )

        edge_traces.append(edge_trace)

    degrees = dict(G.degree())

    node_x = []
    node_y = []
    node_z = []
    node_text = []
    node_color = []
    node_size = []
    node_symbol = []

    for node, data in G.nodes(data=True):
        x, y, z = pos[node]

        node_type = data.get("node_type", "unknown")
        degree = degrees.get(node, 0)

        base_size = NODE_SIZES.get(node_type, NODE_SIZES["unknown"])
        final_size = min(base_size + degree * 0.7, 38)

        node_x.append(x)
        node_y.append(y)
        node_z.append(z)

        node_color.append(NODE_COLORS.get(node_type, NODE_COLORS["unknown"]))
        node_size.append(final_size)
        node_symbol.append(NODE_SYMBOLS.get(node_type, NODE_SYMBOLS["unknown"]))

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
            opacity=0.93,
            line=dict(
                width=0.8,
                color="white",
            ),
        ),
        showlegend=False,
    )

    legend_traces = [
        go.Scatter3d(
            x=[None],
            y=[None],
            z=[None],
            mode="markers",
            marker=dict(
                size=10,
                color=color,
                symbol=NODE_SYMBOLS[node_type],
            ),
            name=node_type.capitalize(),
        )
        for node_type, color in NODE_COLORS.items()
        if node_type != "unknown"
    ]

    fig = go.Figure(
        data=edge_traces + [node_trace] + legend_traces
    )

    fig.update_layout(
        title="3D Variant–Gene–Disease Network",
        height=780,
        margin=dict(l=0, r=0, b=0, t=45),
        showlegend=True,
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor="rgba(0,0,0,0)",
        ),
        scene=dict(
            xaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                title="",
                zeroline=False,
            ),
            yaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                title="",
                zeroline=False,
            ),
            zaxis=dict(
                showbackground=False,
                showgrid=False,
                showticklabels=False,
                title="",
                zeroline=False,
            ),
        ),
    )

    return fig
