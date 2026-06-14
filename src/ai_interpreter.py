from __future__ import annotations

from typing import Any, Generator

import pandas as pd
import streamlit as st

try:
    from groq import Groq
except ImportError:  # pragma: no cover - import guard
    Groq = None  # type: ignore[assignment]


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

# Default model. Override via st.secrets["GROQ_MODEL"] if you want.
DEFAULT_MODEL = "llama-3.3-70b-versatile"

# Quick / cheap model for short tasks like single-gene blurbs.
DEFAULT_FAST_MODEL = "llama-3.1-8b-instant"

# Hard guardrails so a giant network doesn't blow up the context window.
MAX_HUB_GENES_IN_PROMPT = 15
MAX_HUB_VARIANTS_IN_PROMPT = 10
MAX_DISEASES_IN_PROMPT = 10
MAX_COMMUNITIES_IN_PROMPT = 5

SYSTEM_DISCLAIMER = (
    "You are an AI research assistant interpreting GWAS-derived networks. "
    "Your role is to EXPLAIN and SUMMARISE, not to make clinical claims or "
    "fabricate facts. If you are uncertain about a gene or variant function, "
    "say so and suggest the user verify with NCBI Gene / Ensembl / UniProt. "
    "Never claim a finding is clinically actionable."
)


# ---------------------------------------------------------------------------
# Client
# ---------------------------------------------------------------------------


def get_groq_client(api_key: str | None) -> Groq | None:
    """
    Build a Groq client. Returns None if the key is missing or the library
    isn't installed, so the UI can fall back gracefully.
    """

    if Groq is None:
        return None

    # Priority: explicit arg > session_state (sidebar input) > st.secrets.
    key = api_key or ""

    if not key:
        key = st.session_state.get("GROQ_API_KEY", "") or ""

    if not key:
        try:
            key = st.secrets.get("GROQ_API_KEY", "")  # type: ignore[attr-defined]
        except Exception:
            key = ""

    if not key:
        return None

    return Groq(api_key=key)


def get_model_name(fast: bool = False) -> str:
    """
    Resolve model name from secrets, else use default.
    """

    key = "GROQ_MODEL_FAST" if fast else "GROQ_MODEL"

    value = st.session_state.get(key, "") or ""

    if not value:
        try:
            value = st.secrets.get(key, "")  # type: ignore[attr-defined]
        except Exception:
            value = ""

    return value or (DEFAULT_FAST_MODEL if fast else DEFAULT_MODEL)


# ---------------------------------------------------------------------------
# Prompt context packer
#
# We don't ship the full node list to the LLM (expensive + noisy).
# We pack a compact, structured summary that the model can reason over.
# ---------------------------------------------------------------------------


def _format_metrics_df(
    df: pd.DataFrame,
    name_col: str,
    n: int,
) -> str:
    if df is None or df.empty:
        return "  (none)"

    keep_cols = [
        c for c in [
            name_col,
            "degree",
            "degree_centrality",
            "betweenness_centrality",
            "hub_score",
        ] if c in df.columns
    ]

    head = df[keep_cols].head(n)

    lines = [f"  - {row[name_col]} "
             f"(degree={row.get('degree', '?')}, "
             f"hub_score={round(float(row.get('hub_score', 0)), 3)})"
             for _, row in head.iterrows()]

    return "\n".join(lines) if lines else "  (none)"


def build_network_context(
    G,
    summary: dict[str, int],
    top_hub_genes: pd.DataFrame,
    top_hub_variants: pd.DataFrame | None = None,
    top_diseases: pd.DataFrame | None = None,
    filtered_df: pd.DataFrame | None = None,
) -> str:
    """
    Pack the network into a compact text block for the LLM.
    """

    pval_min = "n/a"
    pval_median = "n/a"
    n_unique_genes = "n/a"
    n_unique_snps = "n/a"
    n_unique_diseases = "n/a"

    if filtered_df is not None and not filtered_df.empty:
        if "p_value" in filtered_df.columns:
            pval_min = f"{filtered_df['p_value'].min():.2e}"
            pval_median = f"{filtered_df['p_value'].median():.2e}"
        if "gene" in filtered_df.columns:
            n_unique_genes = filtered_df["gene"].nunique()
        if "snp" in filtered_df.columns:
            n_unique_snps = filtered_df["snp"].nunique()
        if "disease" in filtered_df.columns:
            n_unique_diseases = filtered_df["disease"].nunique()

    parts: list[str] = []

    parts.append("NETWORK OVERVIEW")
    parts.append(
        f"- Total nodes: {summary.get('nodes', 0)}\n"
        f"- Total edges: {summary.get('edges', 0)}\n"
        f"- Variants (SNPs): {summary.get('variants', 0)}\n"
        f"- Genes: {summary.get('genes', 0)}\n"
        f"- Diseases/traits: {summary.get('diseases', 0)}"
    )

    if filtered_df is not None and not filtered_df.empty:
        parts.append("\nDATA STATISTICS")
        parts.append(
            f"- Rows used to build the network: {len(filtered_df)}\n"
            f"- Unique genes in data: {n_unique_genes}\n"
            f"- Unique SNPs in data: {n_unique_snps}\n"
            f"- Unique diseases in data: {n_unique_diseases}\n"
            f"- P-value range: min={pval_min}, median={pval_median}"
        )

    parts.append("\nTOP HUB GENES (ranked by composite hub score)")
    parts.append(_format_metrics_df(
        top_hub_genes, "gene", MAX_HUB_GENES_IN_PROMPT,
    ))

    if top_hub_variants is not None and not top_hub_variants.empty:
        parts.append("\nTOP HUB VARIANTS (SNPs with most connections)")
        parts.append(_format_metrics_df(
            top_hub_variants, "snp", MAX_HUB_VARIANTS_IN_PROMPT,
        ))

    if top_diseases is not None and not top_diseases.empty:
        parts.append("\nTOP CONNECTED DISEASES/TRAITS")
        parts.append(_format_metrics_df(
            top_diseases, "disease", MAX_DISEASES_IN_PROMPT,
        ))

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Low-level streaming call
# ---------------------------------------------------------------------------


def stream_chat(
    client: Groq,
    messages: list[dict[str, str]],
    *,
    model: str | None = None,
    temperature: float = 0.4,
    max_tokens: int = 1500,
) -> Generator[str, None, None]:
    """
    Yield text chunks from a streaming Groq chat completion.
    """

    chosen_model = model or get_model_name()

    stream = client.chat.completions.create(
        model=chosen_model,
        messages=messages,
        temperature=temperature,
        max_tokens=max_tokens,
        stream=True,
    )

    for chunk in stream:
        try:
            delta = chunk.choices[0].delta
            piece = getattr(delta, "content", None)
        except (IndexError, AttributeError):
            piece = None

        if piece:
            yield piece


# ---------------------------------------------------------------------------
# Public feature: 1) Network interpretation
# ---------------------------------------------------------------------------


def render_network_interpretation(
    client: Groq,
    context: str,
    trait_hint: str = "",
) -> None:
    """
    Stream a plain-English interpretation of the network.
    """

    user_prompt = (
        "Interpret the following GWAS-derived variant–gene–disease network "
        "for a research audience. Cover:\n"
        "1. What the network shape suggests (dense hubs vs. sparse, modular vs. central).\n"
        "2. What the top hub genes and variants imply biologically (pathways, "
        "   pleiotropy) — only what is well-established. If unsure, say so.\n"
        "3. Any caveats about the data (small N, single trait, missing context).\n"
        "Keep it 4–6 short paragraphs. No bullet spam. No markdown headers.\n\n"
        f"{('Reported trait/query: ' + trait_hint + chr(10)) if trait_hint else ''}"
        f"NETWORK CONTEXT:\n{context}"
    )

    messages = [
        {"role": "system", "content": SYSTEM_DISCLAIMER},
        {"role": "user", "content": user_prompt},
    ]

    st.write_stream(stream_chat(client, messages))


# ---------------------------------------------------------------------------
# Public feature: 2) Per-gene explainer
# ---------------------------------------------------------------------------


def render_gene_explanation(
    client: Groq,
    gene_name: str,
    network_context: str,
) -> None:
    """
    Stream a short explainer for a single gene in the context of the network.
    """

    user_prompt = (
        f"Explain the gene **{gene_name}** in 2–3 short paragraphs for a "
        "researcher who is exploring a GWAS network. Cover:\n"
        "- Known biological function (consensus only; if contested, say so).\n"
        "- Why it might appear as a hub in a GWAS network "
        "(pathways, pleiotropy, common variant burden).\n"
        "- One practical suggestion for the user to verify "
        "(e.g. 'look it up on NCBI Gene / Ensembl / OpenTargets').\n\n"
        "If you don't know the gene, say so explicitly. Do not invent.\n\n"
        f"NETWORK CONTEXT (for grounding only):\n{network_context}"
    )

    messages = [
        {"role": "system", "content": SYSTEM_DISCLAIMER},
        {"role": "user", "content": user_prompt},
    ]

    st.write_stream(stream_chat(client, messages, max_tokens=600))


# ---------------------------------------------------------------------------
# Public feature: 3) Chat with the network
# ---------------------------------------------------------------------------


CHAT_SYSTEM_TEMPLATE = (
    SYSTEM_DISCLAIMER
    + "\n\nYou can answer questions about the GWAS network below. The user "
      "may ask about specific genes, variants, diseases, or ask for "
      "comparisons. Always ground your answer in the NETWORK CONTEXT when "
      "possible, and call out anything that the data does NOT support."
)


def answer_network_question(
    client: Groq,
    chat_history: list[dict[str, str]],
    user_question: str,
    network_context: str,
) -> Generator[str, None, None]:
    """
    Stream an answer to a free-form question grounded in the network context.
    `chat_history` is the prior turns in OpenAI chat format.
    """

    messages: list[dict[str, str]] = [
        {"role": "system", "content": CHAT_SYSTEM_TEMPLATE},
        {"role": "system", "content": f"NETWORK CONTEXT:\n{network_context}"},
        *chat_history,
        {"role": "user", "content": user_question},
    ]

    return stream_chat(client, messages, max_tokens=900)


# ---------------------------------------------------------------------------
# Public feature: 4) Auto-generated markdown report
# ---------------------------------------------------------------------------


REPORT_SYSTEM = (
    "You write concise, neutral, research-style summaries of GWAS network "
    "analyses. Markdown only. No clinical claims. No fabricated gene functions."
)


def build_markdown_report(
    client: Groq,
    network_context: str,
    trait_hint: str = "",
) -> str:
    """
    Return a markdown report as a single string. Non-streaming (we need the
    whole thing for the download button).
    """

    user_prompt = (
        "Write a markdown report on the GWAS-derived network below. Structure:\n\n"
        "# GWAS Network Analysis Report\n"
        "## 1. Summary\n"
        "## 2. Top hub genes\n"
        "## 3. Top hub variants\n"
        "## 4. Disease / trait landscape\n"
        "## 5. Biological interpretation\n"
        "## 6. Caveats and next steps\n\n"
        "Keep each section 2–4 sentences. Use bullet lists sparingly. "
        "Be honest about uncertainty.\n\n"
        f"{('Reported trait/query: ' + trait_hint + chr(10)) if trait_hint else ''}"
        f"NETWORK CONTEXT:\n{network_context}"
    )

    chosen_model = get_model_name()
    response = client.chat.completions.create(
        model=chosen_model,
        messages=[
            {"role": "system", "content": REPORT_SYSTEM},
            {"role": "user", "content": user_prompt},
        ],
        temperature=0.3,
        max_tokens=1800,
    )

    return response.choices[0].message.content or ""


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------


def ai_unavailable_message() -> None:
    st.info(
        "🤖 AI interpretation is unavailable. Add a Groq API key in the sidebar "
        "(or set `GROQ_API_KEY` in `.streamlit/secrets.toml`) and rebuild."
    )


def check_ai_available() -> tuple[bool, str]:
    """
    Return (ok, message). `ok` is True only if the groq library is installed
    AND a key is available.
    """

    if Groq is None:
        return False, (
            "The `groq` Python package is not installed. "
            "Run `pip install groq` to enable AI features."
        )

    key = st.session_state.get("GROQ_API_KEY", "") or ""

    if not key:
        try:
            key = st.secrets.get("GROQ_API_KEY", "")  # type: ignore[attr-defined]
        except Exception:
            key = ""

    if not key:
        return False, "No Groq API key found in secrets or sidebar."

    return True, "ok"
