from __future__ import annotations

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

<<<<<<< HEAD
from src.config import create_directories, GROQ_API_KEY
=======
from src.config import create_directories, get_groq_api_key
>>>>>>> e15f54a2b5abd53396b25f22282ef45f72f14a7d
from src.gwas_fetcher import fetch_clean_gwas_network_data
from src.data_loader import (
    read_uploaded_file,
    validate_network_input,
    filter_by_pvalue,
    keep_genes_with_min_snps,
)
from src.graph_builder import build_variant_gene_disease_graph
from src.graph_analyzer import (
    get_graph_summary,
    get_top_hub_genes,
    get_top_hub_variants,
    get_top_connected_diseases,
)
from src.graph_3d_visualizer import create_3d_network_figure
from src import ai_interpreter


st.set_page_config(
    page_title="GWAS 3D Network Explorer",
    layout="wide",
    initial_sidebar_state="expanded",
)

create_directories()


<<<<<<< HEAD
CUSTOM_CSS = """
<style>
.block-container {
    padding-top: 2.5rem;
    max-width: 1400px;
}
=======
# ---------------------------------------------------------------------------
# Floating AI assistant
#
# The popover button itself is the floating robot button. We don't try to
# inject a separate toolbar button and then trigger the popover — that is
# fragile across Streamlit versions because the popover DOM contract changes.
# Instead we render the popover inline and style its trigger button to be
# fixed-positioned at the bottom-right of the viewport, so it floats over
# the rest of the content. The Groq API key is configured in `src/config.py`
# (or via the GROQ_API_KEY env var) — never typed into the web app.
# ---------------------------------------------------------------------------

# Inject CSS once to:
#   1. Make the popover trigger button a fixed-position floating robot button.
#   2. Ensure the popover body sits on top of everything.
# We use components.html() because st.markdown() strips <script>/<style>
# tags, but the iframe can reach the parent document via window.parent —
# the standard Streamlit trick.
if "ai_floating_style_injected" not in st.session_state:
    st.session_state["ai_floating_style_injected"] = True
    components.html(
        """
        <script>
        (function () {
            var p = window.parent.document;

            if (p.getElementById('st-ai-floating-style')) return;

            var style = p.createElement('style');
            style.id = 'st-ai-floating-style';
            style.textContent = `
                /* Floating AI popover trigger — fixed bottom-right */
                button[data-testid="stPopoverButton"][kind="primary"],
                button[data-testid="stPopoverButton"] {
                    position: fixed !important;
                    bottom: 24px !important;
                    right: 24px !important;
                    z-index: 99999 !important;
                    width: 60px !important;
                    height: 60px !important;
                    border-radius: 50% !important;
                    border: none !important;
                    cursor: pointer !important;
                    font-size: 28px !important;
                    line-height: 1 !important;
                    padding: 0 !important;
                    background: linear-gradient(135deg,#4E79A7 0%,#59A14F 100%) !important;
                    color: #fff !important;
                    box-shadow: 0 6px 18px rgba(0,0,0,.30) !important;
                    transition: transform .15s ease,
                                box-shadow .15s ease !important;
                }
                button[data-testid="stPopoverButton"]:hover {
                    transform: translateY(-2px) scale(1.05) !important;
                    box-shadow: 0 10px 24px rgba(0,0,0,.40) !important;
                }
                button[data-testid="stPopoverButton"]:active {
                    transform: translateY(0) scale(1) !important;
                }

                /* Popover panel — float above all content */
                div[data-testid="stPopoverBody"] {
                    z-index: 100000 !important;
                    max-width: 480px !important;
                    max-height: 70vh !important;
                    overflow-y: auto !important;
                }
            `;
            p.head.appendChild(style);
        })();
        </script>
        """,
        height=0,
    )

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
>>>>>>> e15f54a2b5abd53396b25f22282ef45f72f14a7d

.main-title {
    font-size: 2.6rem;
    font-weight: 900;
    color: #ffffff;
    letter-spacing: -0.8px;
    margin-bottom: 0.25rem;
    line-height: 1.1;
}

.subtitle {
    color: #94a3b8;
    font-size: 1.05rem;
    margin-bottom: 1.4rem;
}

[data-testid="stSidebar"] {
    border-right: 1px solid rgba(120, 120, 120, 0.15);
}

.stAlert {
    border-radius: 10px;
}
</style>
"""

<<<<<<< HEAD
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)
=======
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

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
>>>>>>> e15f54a2b5abd53396b25f22282ef45f72f14a7d


def initialize_state() -> None:
    defaults = {
        "ai_chat_history": [],
        "ai_gene_cache": {},
        "ai_report_text": "",
        "ai_context": "",
    }

    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


<<<<<<< HEAD
def reset_ai_state() -> None:
    for key in [
        "ai_chat_history",
        "ai_gene_cache",
        "ai_report_text",
        "ai_context",
    ]:
        st.session_state.pop(key, None)


def build_network_from_settings(
    data_source: str,
    uploaded_file,
    trait: str,
    pvalue_threshold: float,
    min_snps_per_gene: int,
    top_n_rows: int,
    max_pages: int,
    max_rows_to_visualize: int,
) -> None:
    if data_source == "Upload file":
        if uploaded_file is None:
            st.warning("Please upload a CSV, TSV, or TXT file.")
            st.stop()
=======
# ---------------------------------------------------------------------------
# Build network button
# ---------------------------------------------------------------------------
>>>>>>> e15f54a2b5abd53396b25f22282ef45f72f14a7d

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

    if filtered_df.empty:
        st.warning(
            "No data found after p-value filtering. "
            "Try a less strict threshold such as 1e-5 or 1e-3."
        )
        st.stop()

    filtered_df = (
        filtered_df.sort_values("p_value", ascending=True)
        .head(top_n_rows)
        .reset_index(drop=True)
    )

    filtered_df = keep_genes_with_min_snps(
        filtered_df,
        min_snps_per_gene=min_snps_per_gene,
    )

    if filtered_df.empty:
        st.warning("No data left after the minimum SNPs per gene filter.")
        st.stop()

    filtered_df = filtered_df.head(max_rows_to_visualize).reset_index(drop=True)

    G = build_variant_gene_disease_graph(filtered_df)

    st.session_state["graph"] = G
    st.session_state["graph_summary"] = get_graph_summary(G)
    st.session_state["top_genes"] = get_top_hub_genes(G, top_n=15)
    st.session_state["top_variants"] = get_top_hub_variants(G, top_n=10)
    st.session_state["top_diseases"] = get_top_connected_diseases(G, top_n=10)
    st.session_state["filtered_df"] = filtered_df
    st.session_state["trait"] = trait

    reset_ai_state()


def render_sidebar() -> dict:
    st.sidebar.title("Control panel")

    st.sidebar.header("Data source")

    data_source = st.sidebar.radio(
        "Choose data source",
        ["Upload file", "GWAS Catalog API"],
    )

    uploaded_file = None

    if data_source == "Upload file":
        uploaded_file = st.sidebar.file_uploader(
            "Upload CSV / TSV / TXT",
            type=["csv", "tsv", "txt"],
        )

    st.sidebar.header("Filters")

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
        "Top significant rows",
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

    build_clicked = st.sidebar.button(
        "Build 3D network",
        type="primary",
        use_container_width=True,
    )

    return {
        "data_source": data_source,
        "uploaded_file": uploaded_file,
        "trait": trait,
        "pvalue_threshold": pvalue_threshold,
        "min_snps_per_gene": min_snps_per_gene,
        "top_n_rows": top_n_rows,
        "max_pages": max_pages,
        "max_rows_to_visualize": max_rows_to_visualize,
        "build_clicked": build_clicked,
    }


def render_header() -> None:
    col_title, col_ai = st.columns([7, 1])

    with col_title:
        st.markdown(
            """
            <div class="main-title">GWAS 3D Network Explorer</div>
            <div class="subtitle">
                Interactive Variant–Gene–Disease network exploration with AI-assisted biological interpretation.
            </div>
            """,
            unsafe_allow_html=True,
        )

    with col_ai:
        st.write("")
        if "graph" in st.session_state:
            if st.button("🤖 AI", type="primary", use_container_width=True):
                render_ai_dialog()


def render_metrics(summary: dict) -> None:
    col1, col2, col3, col4, col5 = st.columns(5)

    col1.metric("Nodes", summary.get("nodes", 0))
    col2.metric("Edges", summary.get("edges", 0))
    col3.metric("Variants", summary.get("variants", 0))
    col4.metric("Genes", summary.get("genes", 0))
    col5.metric("Diseases", summary.get("diseases", 0))


def render_results() -> None:
    G = st.session_state["graph"]
    summary = st.session_state["graph_summary"]
    top_genes_df = st.session_state["top_genes"]
    top_variants_df = st.session_state["top_variants"]
    top_diseases_df = st.session_state["top_diseases"]
    filtered_df = st.session_state["filtered_df"]

    render_metrics(summary)

    st.subheader("Interactive 3D network")
    fig = create_3d_network_figure(G)
    st.plotly_chart(fig, use_container_width=True)

    tab_genes, tab_variants, tab_diseases, tab_data = st.tabs(
        [
            "Top hub genes",
            "Top hub variants",
            "Connected diseases",
            "Filtered data",
        ]
    )

<<<<<<< HEAD
    with tab_genes:
        st.dataframe(top_genes_df, use_container_width=True)

    with tab_variants:
        st.dataframe(top_variants_df, use_container_width=True)

    with tab_diseases:
        st.dataframe(top_diseases_df, use_container_width=True)

    with tab_data:
        st.dataframe(filtered_df, use_container_width=True)

        csv_data = filtered_df.to_csv(index=False).encode("utf-8")

        st.download_button(
            label="Download filtered associations CSV",
            data=csv_data,
            file_name="filtered_gwas_network_data.csv",
            mime="text/csv",
            use_container_width=True,
        )


def get_ai_context() -> str:
    if st.session_state.get("ai_context"):
        return st.session_state["ai_context"]

    context = ai_interpreter.build_network_context(
        st.session_state["graph"],
        st.session_state["graph_summary"],
        st.session_state["top_genes"],
        top_hub_variants=st.session_state["top_variants"],
        top_diseases=st.session_state["top_diseases"],
        filtered_df=st.session_state["filtered_df"],
    )

    st.session_state["ai_context"] = context
    return context


@st.dialog("GWAS AI Assistant", width="large")
def render_ai_dialog() -> None:
    if "graph" not in st.session_state:
        st.info("Build a network first.")
        return

    if not GROQ_API_KEY:
        st.info(
            "AI is disabled because `GROQ_API_KEY` is missing from `.streamlit/secrets.toml`."
        )
        return

    client = ai_interpreter.get_groq_client(GROQ_API_KEY)

    if client is None:
        ai_interpreter.ai_unavailable_message()
        return

    network_context = get_ai_context()
    trait_hint = st.session_state.get("trait", "")

    tab_interpret, tab_gene, tab_chat, tab_report = st.tabs(
        [
            "Interpret network",
            "Explain gene",
            "Ask network",
            "Report",
        ]
    )

    with tab_interpret:
        st.caption(
            "AI-generated interpretation. Verify gene and variant claims before using them in a report."
        )

        if st.button("Generate interpretation", key="generate_ai_interpretation"):
            ai_interpreter.render_network_interpretation(
                client,
                network_context,
                trait_hint=trait_hint,
            )

    with tab_gene:
        top_genes_df = st.session_state["top_genes"]

        gene_options = (
            top_genes_df["gene"].tolist()
            if top_genes_df is not None and not top_genes_df.empty
            else []
        )

        selected_gene = st.selectbox(
            "Select hub gene",
            options=gene_options if gene_options else [""],
            disabled=not gene_options,
        )

        custom_gene = st.text_input(
            "Or type another gene",
            placeholder="APOE",
        )

        gene_to_explain = custom_gene.strip().upper() or selected_gene

        if st.button("Explain gene", key="explain_gene_button"):
            if not gene_to_explain:
                st.warning("Please select or type a gene.")
            else:
                ai_interpreter.render_gene_explanation(
                    client,
                    gene_to_explain,
                    network_context,
                )

    with tab_chat:
        st.session_state.setdefault("ai_chat_history", [])

        for turn in st.session_state["ai_chat_history"]:
            with st.chat_message(turn["role"]):
                st.markdown(turn["content"])

        user_question = st.text_input(
            "Ask about this network",
            placeholder="Example: Which hub gene seems most central?",
            key="ai_question_input",
        )

        if st.button("Ask AI", key="ask_ai_button"):
            if not user_question.strip():
                st.warning("Type a question first.")
            else:
                st.session_state["ai_chat_history"].append(
                    {"role": "user", "content": user_question}
                )

                with st.chat_message("user"):
                    st.markdown(user_question)

                with st.chat_message("assistant"):
                    stream = ai_interpreter.answer_network_question(
                        client,
                        st.session_state["ai_chat_history"][:-1],
                        user_question,
                        network_context,
=======
    # -----------------------------------------------------------------------
    # AI INTERPRETATION SECTION — rendered inside a floating popover that
    # opens when the user clicks the 🤖 robot button. The Groq API key is
    # configured in `src/config.py` (or via the GROQ_API_KEY env var), not
    # typed into the web app.
    # -----------------------------------------------------------------------

    # Render the robot trigger button. CSS above makes it position:fixed
    # so it always floats at the bottom-right corner regardless of scroll.
    with st.popover(
        "🤖",
        help="Open the AI assistant",
        use_container_width=False,
    ):
        st.markdown("### 🤖 AI Assistant")
        st.caption(
            "Powered by Groq. API key is configured in `src/config.py`."
        )

        api_key = get_groq_api_key()
        if not api_key:
            st.info(
                "No Groq API key configured. Set `GROQ_API_KEY` in "
                "`src/config.py` or as an environment variable, then "
                "rebuild the network. Get a free key at "
                "[console.groq.com](https://console.groq.com)."
            )
        else:
            # Build the shared context once per graph run.
            if "ai_context" not in st.session_state:
                st.session_state["ai_context"] = (
                    ai_interpreter.build_network_context(
                        G,
                        summary,
                        top_genes_df,
                        top_hub_variants=top_variants_df,
                        top_diseases=top_diseases_df,
                        filtered_df=filtered_df,
                    )
                )

            network_context = st.session_state["ai_context"]
            client = ai_interpreter.get_groq_client(api_key)

            if client is None:
                ai_interpreter.ai_unavailable_message()
            else:
                tab_interpret, tab_genes, tab_chat, tab_report = st.tabs([
                    "🧠 Interpret",
                    "🧬 Gene",
                    "💬 Chat",
                    "📄 Report",
                ])

                # ---- Tab 1: Interpret the network ---------------------
                with tab_interpret:
                    st.caption(
                        "Plain-English narrative of the network "
                        "structure, hub biology, and caveats. Generated "
                        "by an LLM — verify specific gene/variant "
                        "claims with NCBI Gene, Ensembl, or OpenTargets."
>>>>>>> e15f54a2b5abd53396b25f22282ef45f72f14a7d
                    )
                    answer = st.write_stream(stream)

<<<<<<< HEAD
                st.session_state["ai_chat_history"].append(
                    {"role": "assistant", "content": answer}
                )

        if st.session_state["ai_chat_history"]:
            if st.button("Clear chat", key="clear_ai_chat"):
                st.session_state["ai_chat_history"] = []
                st.rerun()

    with tab_report:
        if st.button("Generate markdown report", key="generate_markdown_report"):
            with st.spinner("Generating report..."):
                st.session_state["ai_report_text"] = (
                    ai_interpreter.build_markdown_report(
                        client,
                        network_context,
                        trait_hint=trait_hint,
                    )
                )

        report_text = st.session_state.get("ai_report_text", "")

        if report_text:
            st.markdown(report_text)

            st.download_button(
                label="Download report",
                data=report_text.encode("utf-8"),
                file_name="gwas_network_report.md",
                mime="text/markdown",
                use_container_width=True,
            )


def main() -> None:
    initialize_state()

    settings = render_sidebar()

    if settings["build_clicked"]:
        try:
            build_network_from_settings(
                **{
                    key: value
                    for key, value in settings.items()
                    if key != "build_clicked"
                }
            )
            st.success("Network successfully generated.")

        except Exception as error:
            st.error(f"Network generation failed: {error}")
            st.stop()

    render_header()

    if "graph" in st.session_state:
        render_results()
    else:
        st.info(
            "Configure the data source and filters in the sidebar, then click **Build 3D network**."
        )


if __name__ == "__main__":
    main()
=======
                    if st.button(
                        "✨ Generate interpretation",
                        key="btn_interpret",
                    ):
                        with st.spinner("Interpreting network…"):
                            st.session_state["ai_interpretation"] = (
                                ai_interpreter.render_network_interpretation(
                                    client,
                                    network_context,
                                    trait_hint=trait_hint,
                                )
                            )

                    if st.session_state.get("ai_interpretation"):
                        st.markdown("#### Interpretation")
                        st.markdown(
                            st.session_state["ai_interpretation"]
                        )

                # ---- Tab 2: Explain a single gene ---------------------
                with tab_genes:
                    st.caption(
                        "Pick a hub gene (or type any gene symbol) and "
                        "the LLM will explain its function and why it "
                        "might be a hub."
                    )

                    if top_genes_df is not None and not top_genes_df.empty:
                        gene_options = top_genes_df["gene"].tolist()
                    else:
                        gene_options = []

                    col_a, col_b = st.columns([2, 1])
                    with col_a:
                        selected_gene = st.selectbox(
                            "Hub gene",
                            options=gene_options if gene_options else [""],
                            disabled=not gene_options,
                        )
                    with col_b:
                        custom_gene = st.text_input(
                            "…or type another",
                            placeholder="e.g. BRCA1",
                        )

                    gene_to_explain = (
                        custom_gene.strip().upper() or selected_gene
                    )

                    if st.button(
                        "🔍 Explain gene",
                        key="btn_explain_gene",
                    ):
                        if not gene_to_explain:
                            st.warning("Pick or type a gene symbol first.")
                        else:
                            with st.spinner(
                                f"Looking up {gene_to_explain}…"
                            ):
                                explanation = (
                                    ai_interpreter.render_gene_explanation(
                                        client,
                                        gene_to_explain,
                                        network_context,
                                    )
                                )
                                st.session_state.setdefault(
                                    "ai_gene_cache", {},
                                )[gene_to_explain] = explanation

                    cache = st.session_state.get("ai_gene_cache", {})
                    if gene_to_explain and gene_to_explain in cache:
                        st.markdown(f"#### {gene_to_explain}")
                        st.markdown(cache[gene_to_explain])

                # ---- Tab 3: Chat with the network ---------------------
                with tab_chat:
                    st.caption(
                        "Ask free-form questions about the network. "
                        "The LLM is grounded in the network stats above "
                        "but may still get specific biology wrong — "
                        "verify before publishing."
                    )

                    st.session_state.setdefault("ai_chat_history", [])

                    for turn in st.session_state["ai_chat_history"]:
                        with st.chat_message(turn["role"]):
                            st.markdown(turn["content"])

                    user_msg = st.chat_input(
                        "e.g. Which gene links the most diseases?",
                        key="ai_chat_input",
                    )

                    if user_msg:
                        st.session_state["ai_chat_history"].append(
                            {"role": "user", "content": user_msg},
                        )
                        with st.chat_message("user"):
                            st.markdown(user_msg)

                        with st.chat_message("assistant"):
                            prior = st.session_state[
                                "ai_chat_history"
                            ][:-1]
                            pieces = []
                            gen = ai_interpreter.answer_network_question(
                                client,
                                prior,
                                user_msg,
                                network_context,
                            )
                            for piece in gen:
                                pieces.append(piece)
                            full_answer = "".join(pieces)
                            st.markdown(full_answer)

                        st.session_state["ai_chat_history"].append(
                            {
                                "role": "assistant",
                                "content": full_answer,
                            },
                        )

                        if (
                            len(st.session_state["ai_chat_history"])
                            > 20
                        ):
                            st.session_state["ai_chat_history"] = (
                                st.session_state["ai_chat_history"][-20:]
                            )

                    if st.session_state["ai_chat_history"]:
                        if st.button(
                            "Clear chat", key="btn_clear_chat"
                        ):
                            st.session_state["ai_chat_history"] = []
                            st.rerun()

                # ---- Tab 4: Generate markdown report ------------------
                with tab_report:
                    st.caption(
                        "Generate a markdown report you can download "
                        "and share. Useful for notetaking or as a "
                        "first draft of a discussion section."
                    )

                    if st.button("📝 Generate report", key="btn_report"):
                        with st.spinner("Writing report…"):
                            try:
                                st.session_state["ai_report_text"] = (
                                    ai_interpreter.build_markdown_report(
                                        client,
                                        network_context,
                                        trait_hint=trait_hint,
                                    )
                                )
                            except Exception as error:
                                st.error(
                                    f"Report generation failed: {error}"
                                )

                    report_text = st.session_state.get(
                        "ai_report_text", ""
                    )
                    if report_text:
                        st.markdown(report_text)
                        st.download_button(
                            label="⬇️ Download report (.md)",
                            data=report_text.encode("utf-8"),
                            file_name="gwas_network_report.md",
                            mime="text/markdown",
                            key="btn_download_report",
                        )


else:
    st.info(
        "👈 Configure your data source and filters in the sidebar, "
        "then click **Build 3D network** to start."
    )
>>>>>>> e15f54a2b5abd53396b25f22282ef45f72f14a7d
