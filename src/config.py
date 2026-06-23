from __future__ import annotations
import os
from pathlib import Path


# PROJECT PATHS

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


# ---------------------------------------------------------------------------
# AI / GROQ SETTINGS
#
# Paste your Groq API key below (get one free at https://console.groq.com).
# You can also set the GROQ_API_KEY environment variable — it takes
# precedence over the value defined here.
#
# SECURITY: do NOT commit a real key to a public repo. If you push this
# project to GitHub, keep this file local or load the key from an env var.
# ---------------------------------------------------------------------------

GROQ_API_KEY: str = "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"  # e.g. "gsk_..."

GROQ_MODEL: str = "llama-3.3-70b-versatile"

GROQ_MODEL_FAST: str = "llama-3.1-8b-instant"


def get_groq_api_key() -> str:
    """
    Resolve the Groq API key from, in order of priority:
      1. The GROQ_API_KEY environment variable (best for deployment)
      2. The GROQ_API_KEY value defined in this config file
    Returns an empty string if neither is set.
    """

    return os.environ.get("GROQ_API_KEY", "") or GROQ_API_KEY or ""


# GWAS SETTINGS

GWAS_API_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"

DEFAULT_PVALUE_THRESHOLD = 5e-8

DEFAULT_MAX_API_PAGES = 5

DEFAULT_TOP_ROWS = 500

DEFAULT_MIN_SNPS_PER_GENE = 2


# NETWORK SETTINGS

DEFAULT_LAYOUT_SEED = 42

DEFAULT_LAYOUT_ITERATIONS = 80

DEFAULT_MAX_NODE_SIZE = 35

DEFAULT_MIN_NODE_SIZE = 5



# VISUALIZATION SETTINGS

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

# DIRECTORIES

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
    """
    Create all required project directories.
    """

    for directory in DIRECTORIES:
        directory.mkdir(
            parents=True,
            exist_ok=True,
        )
