from __future__ import annotations
from pathlib import Path
import streamlit as st


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"

RESULTS_DIR = BASE_DIR / "results"

TABLES_DIR = RESULTS_DIR / "tables"
GRAPHS_DIR = RESULTS_DIR / "graphs"
REPORTS_DIR = RESULTS_DIR / "reports"
EXPORTS_DIR = RESULTS_DIR / "exports"

CACHE_DIR = BASE_DIR / "cache"

GWAS_API_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"

DEFAULT_PVALUE_THRESHOLD = 5e-8

DEFAULT_MAX_API_PAGES = 5

DEFAULT_TOP_ROWS = 500

DEFAULT_MIN_SNPS_PER_GENE = 2

DEFAULT_LAYOUT_SEED = 42

DEFAULT_LAYOUT_ITERATIONS = 80

DEFAULT_MAX_NODE_SIZE = 35

DEFAULT_MIN_NODE_SIZE = 5

GROQ_API_KEY = st.secrets.get("GROQ_API_KEY", "")
GROQ_MODEL = st.secrets.get("GROQ_MODEL", "llama-3.3-70b-versatile")
GROQ_MODEL_FAST = st.secrets.get("GROQ_MODEL_FAST", "llama-3.1-8b-instant")


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


DIRECTORIES = [
    DATA_DIR,
    RAW_DIR,
    PROCESSED_DIR,
    RESULTS_DIR,
    TABLES_DIR,
    GRAPHS_DIR,
    REPORTS_DIR,
    EXPORTS_DIR,
    CACHE_DIR,
]


def create_directories() -> None:

    for directory in DIRECTORIES:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )
