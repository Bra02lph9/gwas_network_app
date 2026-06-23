from __future__ import annotations

from typing import Generator

import pandas as pd
import streamlit as st

try:
    from groq import Groq
except ImportError:
    Groq = None


DEFAULT_MODEL = "llama-3.3-70b-versatile"
DEFAULT_FAST_MODEL = "llama-3.1-8b-instant"

MAX_HUB_GENES_IN_PROMPT = 15
MAX_HUB_VARIANTS_IN_PROMPT = 10
MAX_DISEASES_IN_PROMPT = 10


SYSTEM_DISCLAIMER = (
    "You are an AI research assistant specialized in GWAS network interpretation. "
    "You explain only the provided network context. "
    "Do not invent gene functions, disease mechanisms, SNP effects, or clinical conclusions. "
    "Use cautious scientific language such as 'may suggest', 'could indicate', "
    "'is consistent with', or 'requires validation'. "
    "Never claim causality. Never claim clinical actionability. "
    "If information is missing or uncertain, clearly say that external validation "
    "using NCBI Gene, Ensembl, UniProt, GWAS Catalog, or OpenTargets is required."
)


def get_groq_client(api_key: str | None) -> Groq | None:
    """
    Create Groq client from sidebar, session state, or Streamlit secrets.
    """

    if Groq is None:
        return None

    key = api_key or ""

    if not key:
        key = st.session_state.get("GROQ_API_KEY", "") or ""

    if not key:
        try:
            key = st.secrets.get("GROQ_API_KEY", "")
        except Exception:
            key = ""

    if not key:
        return None

    return Groq(api_key=key)


def get_model_name(fast: bool = False) -> str:
    """
    Resolve Groq model name from session state, secrets, or defaults.
    """

    key = "GROQ_MODEL_FAST" if fast else "GROQ_MODEL"

    value = st.session_state.get(key, "") or ""

    if not value:
        try:
            value = st.secrets.get(key, "")
        except Exception:
            value = ""

    return value or (DEFAULT_FAST_MODEL if fast else DEFAULT_MODEL)


def _format_metrics_df(
    df: pd.DataFrame,
    name_col: str,
    n: int,
) -> str:
    """
    Convert top metrics DataFrame into compact text for the LLM.
    """

    if df is None or df.empty:
        return "  (none)"

    keep_cols = [
        c for c in [
            name_col,
            "degree",
            "degree_centrality",
            "betweenness_centrality",
            "closeness_centrality",
            "eigenvector_centrality",
            "hub_score",
        ] if c in df.columns
    ]

    head = df[keep_cols].head(n)

    lines = [
        f"  - {row[name_col]} "
        f"(degree={row.get('degree', '?')}, "
        f"betweenness={round(float(row.get('betweenness_centrality', 0)), 3)}, "
        f"hub_score={round(float(row.get('hub_score', 0)), 3)})"
        for _, row in head.iterrows()
    ]

    return "\n".join(lines) if lines else "  (none)"


def build_network_context(
    G,
    summary: dict,
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
        f"- Diseases/traits: {summary.get('diseases', 0)}\n"
        f"- Connected components: {summary.get('connected_components', 'n/a')}\n"
        f"- Network density: {summary.get('density', 'n/a')}"
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

    parts.append("\nTOP HUB GENES")
    parts.append(_format_metrics_df(
        top_hub_genes,
        "gene",
        MAX_HUB_GENES_IN_PROMPT,
    ))

    if top_hub_variants is not None and not top_hub_variants.empty:
        parts.append("\nTOP HUB VARIANTS")
        parts.append(_format_metrics_df(
            top_hub_variants,
            "snp",
            MAX_HUB_VARIANTS_IN_PROMPT,
        ))

    if top_diseases is not None and not top_diseases.empty:
        parts.append("\nTOP CONNECTED DISEASES / TRAITS")
        parts.append(_format_metrics_df(
            top_diseases,
            "disease",
            MAX_DISEASES_IN_PROMPT,
        ))

    return "\n".join(parts)


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

    try:
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

    except Exception as exc:
        yield (
            "AI interpretation failed. "
            f"Reason: {exc}. "
            "Please check your Groq API key, model name, or internet connection."
        )


def render_network_interpretation(
    client: Groq,
    context: str,
    trait_hint: str = "",
) -> None:
    """
    Stream a biological interpretation of the network.
    """

    user_prompt = (
        "Interpret the following GWAS-derived variant–gene–disease network "
        "for a research audience. Cover:\n"
        "1. What the network topology suggests.\n"
        "2. What the top hub genes and variants may indicate biologically.\n"
        "3. Potential biological pathways only if they are well supported.\n"
        "4. Main limitations and next validation steps.\n\n"
        "Use 4–6 short paragraphs. "
        "Use cautious scientific language. "
        "Do not claim causality. "
        "Do not make clinical recommendations.\n\n"
        f"{('Reported trait/query: ' + trait_hint + chr(10)) if trait_hint else ''}"
        f"NETWORK CONTEXT:\n{context}"
    )

    messages = [
        {"role": "system", "content": SYSTEM_DISCLAIMER},
        {"role": "user", "content": user_prompt},
    ]

    st.write_stream(stream_chat(client, messages))


def render_gene_explanation(
    client: Groq,
    gene_name: str,
    network_context: str,
) -> None:
    """
    Stream a cautious explanation of a selected gene.
    """

    user_prompt = (
        f"Explain the gene {gene_name} in 2–3 short paragraphs for a researcher "
        "exploring a GWAS network. Cover:\n"
        "1. Known biological function, only if well established.\n"
        "2. Why this gene may appear as a hub in the network.\n"
        "3. What the current network does and does not support.\n"
        "4. One practical validation suggestion.\n\n"
        "If you do not know the gene, say so explicitly. "
        "Do not invent functions or mechanisms. "
        "Do not claim causality.\n\n"
        f"NETWORK CONTEXT:\n{network_context}"
    )

    messages = [
        {"role": "system", "content": SYSTEM_DISCLAIMER},
        {"role": "user", "content": user_prompt},
    ]

    st.write_stream(
        stream_chat(
            client,
            messages,
            max_tokens=700,
        )
    )


CHAT_SYSTEM_TEMPLATE = (
    SYSTEM_DISCLAIMER
    + "\n\nYou answer questions about the GWAS network below. "
    "Always ground your answer in the NETWORK CONTEXT when possible. "
    "Clearly state when the data does not support a conclusion."
)


def answer_network_question(
    client: Groq,
    chat_history: list[dict[str, str]],
    user_question: str,
    network_context: str,
) -> Generator[str, None, None]:
    """
    Answer user questions about the current network.
    """

    messages: list[dict[str, str]] = [
        {"role": "system", "content": CHAT_SYSTEM_TEMPLATE},
        {"role": "system", "content": f"NETWORK CONTEXT:\n{network_context}"},
        *chat_history,
        {"role": "user", "content": user_question},
    ]

    return stream_chat(
        client,
        messages,
        max_tokens=900,
    )


REPORT_SYSTEM = (
    "You write concise, neutral, research-style summaries of GWAS network analyses. "
    "Use markdown. Do not make clinical claims. Do not fabricate gene functions. "
    "Use cautious scientific language."
)


def build_markdown_report(
    client: Groq,
    network_context: str,
    trait_hint: str = "",
) -> str:
    """
    Generate a markdown report using Groq.
    """

    user_prompt = (
        "Write a markdown report on the GWAS-derived network below.\n\n"
        "# GWAS Network Analysis Report\n"
        "## 1. Summary\n"
        "## 2. Top hub genes\n"
        "## 3. Top hub variants\n"
        "## 4. Disease / trait landscape\n"
        "## 5. Biological interpretation\n"
        "## 6. Caveats and next steps\n\n"
        "Keep each section 2–4 sentences. "
        "Use bullet lists sparingly. "
        "Be honest about uncertainty. "
        "Do not claim causality or clinical actionability.\n\n"
        f"{('Reported trait/query: ' + trait_hint + chr(10)) if trait_hint else ''}"
        f"NETWORK CONTEXT:\n{network_context}"
    )

    chosen_model = get_model_name()

    try:
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

    except Exception as exc:
        return (
            "# AI Report Generation Failed\n\n"
            f"Reason: {exc}\n\n"
            "Please check your Groq API key, model name, or internet connection."
        )


def ai_unavailable_message() -> None:
    """
    Display AI unavailable message in Streamlit.
    """

    st.info(
        "AI interpretation is unavailable. Add a Groq API key in the sidebar "
        "or set GROQ_API_KEY in `.streamlit/secrets.toml`."
    )


def check_ai_available() -> tuple[bool, str]:
    """
    Check if Groq package and API key are available.
    """

    if Groq is None:
        return False, (
            "The `groq` Python package is not installed. "
            "Run `pip install groq` to enable AI features."
        )

    key = st.session_state.get("GROQ_API_KEY", "") or ""

    if not key:
        try:
            key = st.secrets.get("GROQ_API_KEY", "")
        except Exception:
            key = ""

    if not key:
        return False, "No Groq API key found in secrets or sidebar."

    return True, "ok"
