from __future__ import annotations

from typing import Any
import time

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


GWAS_API_BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"

NETWORK_COLUMNS = [
    "snp",
    "gene",
    "disease",
    "p_value",
    "risk_allele",
    "association_id",
]


def _create_session() -> requests.Session:
    session = requests.Session()

    retries = Retry(
        total=3,
        connect=3,
        read=3,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )

    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    session.headers.update(
        {
            "Accept": "application/json",
            "User-Agent": "GWAS-3D-Network-Explorer/1.0",
        }
    )

    return session


SESSION = _create_session()


def _get_json(
    url: str,
    params: dict[str, Any] | None = None,
    debug: bool = False,
) -> dict[str, Any]:
    try:
        response = SESSION.get(
            url,
            params=params,
            timeout=30,
        )
    except requests.RequestException as exc:
        raise RuntimeError(
            f"GWAS Catalog API connection error: {exc}"
        ) from exc

    if debug:
        print("REQUEST URL:", response.url)
        print("STATUS:", response.status_code)
        print("TEXT:", response.text[:500])

    if response.status_code != 200:
        raise RuntimeError(
            f"GWAS Catalog API error {response.status_code}: "
            f"{response.text[:500]}"
        )

    try:
        return response.json()
    except ValueError as exc:
        raise RuntimeError(
            "GWAS Catalog API returned invalid JSON."
        ) from exc


def _extract_associations_from_response(
    data: dict[str, Any],
) -> list[dict[str, Any]]:
    embedded = data.get("_embedded", {})

    if not isinstance(embedded, dict):
        return []

    for key in [
        "associations",
        "associationDtos",
        "associationDtoes",
        "associationDTOs",
    ]:
        value = embedded.get(key)
        if isinstance(value, list):
            return value

    return []


def fetch_associations_by_trait(
    trait: str,
    max_pages: int = 5,
    page_size: int = 100,
    debug: bool = False,
) -> pd.DataFrame:
    """
    Fetch GWAS Catalog associations using REST API v2.

    The v2 endpoint supports filters such as:
    - efo_trait
    - mapped_gene
    - rs_id
    """

    trait = trait.strip()

    if not trait:
        return pd.DataFrame()

    url = f"{GWAS_API_BASE_URL}/associations"

    records: list[dict[str, Any]] = []

    for page in range(max_pages):
        data = _get_json(
            url,
            params={
                "disease_trait": trait,
                "page": page,
                "size": page_size,
            },
            debug=debug,
        )

        associations = _extract_associations_from_response(data)

        if not associations:
            break

        records.extend(associations)
        time.sleep(0.05)

    return pd.DataFrame(records)


def _split_gene_string(value: str) -> list[str]:
    value = (
        value.replace(" - ", ",")
        .replace(";", ",")
        .replace("/", ",")
        .replace("|", ",")
    )

    return [
        gene.strip().upper()
        for gene in value.split(",")
        if gene.strip()
    ]


def _clean_gene_name(gene: Any) -> str | None:
    if gene is None:
        return None

    gene = str(gene).strip().upper()

    if not gene:
        return None

    if gene.lower() in {"none", "null", "nan", "nr"}:
        return None

    return gene


def _extract_genes_from_any(value: Any) -> list[str]:
    genes: list[str] = []

    if value is None:
        return genes

    if isinstance(value, str):
        genes.extend(_split_gene_string(value))

    elif isinstance(value, list):
        for item in value:
            genes.extend(_extract_genes_from_any(item))

    elif isinstance(value, dict):
        for key in [
            "geneName",
            "ensemblGeneName",
            "geneSymbol",
            "symbol",
            "mappedGene",
            "mappedGenes",
            "authorReportedGene",
            "authorReportedGenes",
            "gene",
            "genes",
        ]:
            if key in value:
                genes.extend(_extract_genes_from_any(value.get(key)))

    clean = sorted(
        {
            gene
            for gene in (_clean_gene_name(g) for g in genes)
            if gene is not None
        }
    )

    return clean


def _extract_snps_from_any(value: Any) -> list[tuple[str, str]]:
    """
    Return list of (rsid, risk_allele).
    """

    snps: list[tuple[str, str]] = []

    if value is None:
        return snps

    if isinstance(value, str):
        rs_id = value.split("-")[0].strip()
        if rs_id.startswith("rs"):
            snps.append((rs_id, value.strip()))

    elif isinstance(value, list):
        for item in value:
            snps.extend(_extract_snps_from_any(item))

    elif isinstance(value, dict):
        rs_id = (
            value.get("rs_id")
            or value.get("rsId")
            or value.get("rsid")
            or value.get("variantId")
            or value.get("variantName")
            or value.get("snp")
        )

        risk_allele = (
            value.get("riskAlleleName")
            or value.get("snp_effect_allele")
            or value.get("effect_allele")
            or value.get("riskAllele")
            or rs_id
        )

        if rs_id:
            rs_id_clean = str(rs_id).split("-")[0].strip()
            if rs_id_clean.startswith("rs"):
                snps.append((rs_id_clean, str(risk_allele).strip()))

        for key in [
            "strongestRiskAlleles",
            "snp_allele",
            "snp_effect_allele",
            "riskAlleles",
            "variants",
        ]:
            if key in value:
                snps.extend(_extract_snps_from_any(value.get(key)))

    unique: dict[str, str] = {}
    for rs_id, risk_allele in snps:
        unique.setdefault(rs_id, risk_allele)

    return sorted(unique.items())


def _get_first_value(
    data: dict[str, Any],
    keys: list[str],
) -> Any:
    for key in keys:
        value = data.get(key)
        if value not in [None, "", [], {}]:
            return value
    return None


def _get_association_pvalue(association: dict[str, Any]) -> Any:
    return _get_first_value(
        association,
        [
            "p_value",
            "pvalue",
            "pValue",
            "p-value",
        ],
    )


def _get_association_id(association: dict[str, Any]) -> Any:
    return _get_first_value(
        association,
        [
            "association_id",
            "associationId",
            "accessionId",
            "id",
        ],
    )


def _get_association_trait(
    association: dict[str, Any],
    fallback_trait: str,
) -> str:
    value = _get_first_value(
        association,
        [
            "disease_trait",
            "diseaseTrait",
            "trait",
            "reportedTrait",
            "mappedTrait",
        ],
    )

    if value:
        return str(value)

    efo_traits = association.get("efoTraits")
    if isinstance(efo_traits, list) and efo_traits:
        first = efo_traits[0]
        if isinstance(first, dict):
            label = (
                first.get("trait")
                or first.get("shortForm")
                or first.get("label")
                or first.get("efoTrait")
            )
            if label:
                return str(label)
        elif isinstance(first, str):
            return first

    return fallback_trait


def extract_gwas_network_table(
    associations_df: pd.DataFrame,
    disease_name: str,
) -> pd.DataFrame:
    """
    Convert GWAS Catalog API v2 associations into clean graph table.
    """

    if associations_df is None or associations_df.empty:
        return pd.DataFrame(columns=NETWORK_COLUMNS)

    rows: list[dict[str, Any]] = []

    for _, association in associations_df.iterrows():
        association_dict = association.to_dict()

        association_id = _get_association_id(association_dict)
        p_value = _get_association_pvalue(association_dict)
        disease = _get_association_trait(association_dict, disease_name)

        genes: list[str] = []
        snps: list[tuple[str, str]] = []

        for key in [
            "mapped_genes",
            "mappedGenes",
            "authorReportedGenes",
            "author_reported_genes",
            "gene",
            "genes",
        ]:
            genes.extend(_extract_genes_from_any(association_dict.get(key)))

        for key in [
            "snp_allele",
            "snp_effect_allele",
            "strongestRiskAlleles",
            "riskAlleles",
            "loci",
            "variant",
            "variants",
            "rs_id",
            "rsId",
            "rsid",
        ]:
            snps.extend(_extract_snps_from_any(association_dict.get(key)))

        genes = sorted(set(genes))

        unique_snps: dict[str, str] = {}
        for rs_id, risk_allele in snps:
            unique_snps.setdefault(rs_id, risk_allele)

        if not genes or not unique_snps or p_value is None:
            continue

        for gene in genes:
            for snp, risk_allele in unique_snps.items():
                rows.append(
                    {
                        "snp": snp,
                        "gene": gene,
                        "disease": disease,
                        "p_value": p_value,
                        "risk_allele": risk_allele,
                        "association_id": association_id,
                    }
                )

    df = pd.DataFrame(rows, columns=NETWORK_COLUMNS)

    if df.empty:
        return df

    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    df = df.dropna(subset=["snp", "gene", "disease", "p_value"])
    df = df[df["p_value"] > 0]
    df = df[df["p_value"] <= 1]
    df = df[df["snp"].astype(str).str.startswith("rs")]
    df = df[df["gene"].astype(str).str.upper() != "UNKNOWN_GENE"]

    return (
        df.drop_duplicates()
        .sort_values("p_value", ascending=True)
        .reset_index(drop=True)
    )


def filter_significant_associations(
    df: pd.DataFrame,
    pvalue_threshold: float = 5e-8,
) -> pd.DataFrame:
    if df.empty:
        return df

    df = df.copy()
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
    df = df.dropna(subset=["p_value"])

    return (
        df[df["p_value"] <= pvalue_threshold]
        .sort_values("p_value", ascending=True)
        .reset_index(drop=True)
    )


def fetch_clean_gwas_network_data(
    trait: str,
    max_pages: int = 5,
    pvalue_threshold: float | None = 5e-8,
    use_demo_if_empty: bool = False,
    debug: bool = False,
) -> pd.DataFrame:
    """
    Main function used by Streamlit app.

    trait -> raw GWAS associations -> clean SNP-Gene-Disease table.
    """

    raw_df = fetch_associations_by_trait(
        trait=trait,
        max_pages=max_pages,
        debug=debug,
    )

    if debug:
        print("RAW API rows:", raw_df.shape)
        print("RAW API columns:", raw_df.columns.tolist())
        if not raw_df.empty:
            print("FIRST RAW ROW:")
            print(raw_df.iloc[0].to_dict())

    clean_df = extract_gwas_network_table(
        associations_df=raw_df,
        disease_name=trait,
    )

    if debug:
        print("CLEAN rows before p-value filter:", clean_df.shape)
        if not clean_df.empty:
            print("Min p-value:", clean_df["p_value"].min())
            print("Max p-value:", clean_df["p_value"].max())
            print("Rows <= 5e-8:", (clean_df["p_value"] <= 5e-8).sum())
            print("Rows <= 1e-5:", (clean_df["p_value"] <= 1e-5).sum())
            print("Rows <= 1e-3:", (clean_df["p_value"] <= 1e-3).sum())

    if pvalue_threshold is not None:
        clean_df = filter_significant_associations(
            clean_df,
            pvalue_threshold=pvalue_threshold,
        )

    if clean_df.empty and use_demo_if_empty:
        clean_df = pd.DataFrame(
            [
                {
                    "snp": "rs1333049",
                    "gene": "CDKN2B-AS1",
                    "disease": trait,
                    "p_value": 1e-20,
                    "risk_allele": "rs1333049-C",
                    "association_id": "demo_1",
                },
                {
                    "snp": "rs4977574",
                    "gene": "CDKN2B-AS1",
                    "disease": trait,
                    "p_value": 2e-18,
                    "risk_allele": "rs4977574-G",
                    "association_id": "demo_2",
                },
            ],
            columns=NETWORK_COLUMNS,
        )

    return clean_df.reset_index(drop=True)


if __name__ == "__main__":
    df = fetch_clean_gwas_network_data(
        trait="coronary artery disease",
        max_pages=5,
        pvalue_threshold=1e-5,
        use_demo_if_empty=False,
        debug=True,
    )

    print("FINAL:", df.shape)
    print(df.head())
