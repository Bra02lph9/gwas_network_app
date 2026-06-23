from __future__ import annotations
import streamlit as st

from src.config import create_directories, GROQ_API_KEY
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


CUSTOM_CSS = """
<style>
.block-container {
    padding-top: 2.5rem;
    max-width: 1400px;
}

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

st.markdown(CUSTOM_CSS, unsafe_allow_html=True)


def initialize_state() -> None:
    defaults = {
        "ai_chat_history": [],
        "ai_gene_cache": {},
        "ai_report_text": "",
        "ai_context": "",
    }

    for key, value in defaults.items():
        st.session_state.setdefault(key, value)


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
                    )
                    answer = st.write_stream(stream)

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
