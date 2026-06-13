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

    possible_keys = [
        "associations",
        "associationDtoes",
        "associationDtos",
    ]

    for key in possible_keys:
        value = embedded.get(key)
        if isinstance(value, list):
            return value

    return []


def _trait_matches_association(
    association: dict[str, Any],
    trait: str,
) -> bool:
    trait_lower = trait.lower().strip()

    searchable_parts: list[str] = []

    for key in [
        "trait",
        "diseaseTrait",
        "reportedTrait",
        "mappedTrait",
    ]:
        value = association.get(key)
        if value:
            searchable_parts.append(str(value))

    efo_traits = association.get("efoTraits", [])

    if isinstance(efo_traits, list):
        for item in efo_traits:
            if isinstance(item, dict):
                searchable_parts.extend(
                    str(value)
                    for value in item.values()
                    if value
                )
            elif isinstance(item, str):
                searchable_parts.append(item)

    searchable_text = " ".join(searchable_parts).lower()

    return trait_lower in searchable_text


def fetch_associations_by_trait(
    trait: str,
    max_pages: int = 5,
    page_size: int = 100,
) -> pd.DataFrame:
    """
    Fetch raw GWAS Catalog associations using REST API v2.
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
                "efo_trait": trait,
                "page": page,
                "size": page_size,
            },
        )

        associations = data.get("_embedded", {}).get("associations", [])

        if not associations:
            break

        records.extend(associations)

        time.sleep(0.05)

    return pd.DataFrame(records)


def _split_gene_string(
    value: str,
) -> list[str]:
    value = (
        value.replace(" - ", ",")
        .replace(";", ",")
        .replace("/", ",")
    )

    genes = [
        gene.strip()
        for gene in value.split(",")
        if gene.strip()
    ]

    return genes


def _extract_gene_names_from_locus(
    locus: dict[str, Any],
) -> list[str]:
    genes: list[str] = []

    for key in [
        "authorReportedGenes",
        "mappedGenes",
        "gene",
        "genes",
    ]:
        values = locus.get(key, [])

        if isinstance(values, list):
            for item in values:
                if isinstance(item, dict):
                    gene = (
                        item.get("geneName")
                        or item.get("ensemblGeneName")
                        or item.get("geneSymbol")
                        or item.get("symbol")
                    )

                    if gene:
                        genes.append(str(gene).strip())

                elif isinstance(item, str):
                    genes.extend(
                        _split_gene_string(item)
                    )

        elif isinstance(values, str):
            genes.extend(
                _split_gene_string(values)
            )

    clean_genes = sorted(
        {
            gene.upper()
            for gene in genes
            if gene and gene.lower() not in {"none", "null", "nan"}
        }
    )

    return clean_genes or ["UNKNOWN_GENE"]


def _extract_risk_alleles_from_locus(
    locus: dict[str, Any],
) -> list[str]:
    risk_alleles: list[str] = []

    values = locus.get("strongestRiskAlleles", [])

    if isinstance(values, list):
        for item in values:
            if isinstance(item, dict):
                value = (
                    item.get("riskAlleleName")
                    or item.get("variantName")
                    or item.get("snp")
                )

                if value:
                    risk_alleles.append(str(value).strip())

            elif isinstance(item, str):
                risk_alleles.append(item.strip())

    return risk_alleles


def _extract_snps_from_locus(
    locus: dict[str, Any],
) -> list[str]:
    snps: list[str] = []

    risk_alleles = _extract_risk_alleles_from_locus(locus)

    for risk_allele in risk_alleles:
        snp = str(risk_allele).split("-")[0].strip()

        if snp.startswith("rs"):
            snps.append(snp)

    for key in [
        "variantId",
        "variantName",
        "snp",
        "rsId",
        "rsid",
    ]:
        value = locus.get(key)

        if isinstance(value, str):
            snp = value.split("-")[0].strip()

            if snp.startswith("rs"):
                snps.append(snp)

    return sorted(set(snps))


def _get_association_pvalue(
    association: dict[str, Any],
) -> Any:
    return (
        association.get("pvalue")
        or association.get("pValue")
        or association.get("p_value")
    )


def _get_association_id(
    association: dict[str, Any],
) -> Any:
    return (
        association.get("associationId")
        or association.get("accessionId")
        or association.get("id")
    )


def extract_gwas_network_table(
    associations_df: pd.DataFrame,
    disease_name: str,
) -> pd.DataFrame:
    """
    Convert GWAS Catalog API v2 associations into clean graph table.
    """

    if associations_df.empty:
        return pd.DataFrame(columns=NETWORK_COLUMNS)

    rows: list[dict[str, Any]] = []

    for _, association in associations_df.iterrows():
        association_dict = association.to_dict()

        association_id = association_dict.get("association_id")
        p_value = association_dict.get("p_value")

        mapped_genes = association_dict.get("mapped_genes", [])
        snp_allele = association_dict.get("snp_allele", [])
        snp_effect_allele = association_dict.get("snp_effect_allele", [])

        if not isinstance(mapped_genes, list):
            mapped_genes = []

        if not isinstance(snp_allele, list):
            snp_allele = []

        snps = []

        for item in snp_allele:
            if isinstance(item, dict):
                rs_id = item.get("rs_id")
                effect_allele = item.get("effect_allele")

                if rs_id:
                    risk_allele = (
                        f"{rs_id}-{effect_allele}"
                        if effect_allele
                        else rs_id
                    )
                    snps.append((rs_id, risk_allele))

        if not snps and isinstance(snp_effect_allele, list):
            for allele in snp_effect_allele:
                if isinstance(allele, str):
                    rs_id = allele.split("-")[0]
                    if rs_id.startswith("rs"):
                        snps.append((rs_id, allele))

        for gene in mapped_genes:
            gene = str(gene).strip().upper()

            if not gene:
                continue

            for snp, risk_allele in snps:
                rows.append(
                    {
                        "snp": snp,
                        "gene": gene,
                        "disease": disease_name,
                        "p_value": p_value,
                        "risk_allele": risk_allele,
                        "association_id": association_id,
                    }
                )

    df = pd.DataFrame(rows, columns=NETWORK_COLUMNS)

    if df.empty:
        return df

    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    df = df.dropna(subset=["snp", "gene", "p_value"])
    df = df[df["p_value"] > 0]
    df = df[df["snp"].astype(str).str.startswith("rs")]

    return df.drop_duplicates().reset_index(drop=True)


def filter_significant_associations(
    df: pd.DataFrame,
    pvalue_threshold: float = 5e-8,
) -> pd.DataFrame:
    """
    Keep only associations below the selected p-value threshold.
    """

    if df.empty:
        return df

    df = df.copy()

    df["p_value"] = pd.to_numeric(
        df["p_value"],
        errors="coerce",
    )

    df = df.dropna(
        subset=["p_value"]
    )

    return (
        df[df["p_value"] <= pvalue_threshold]
        .reset_index(drop=True)
    )


def fetch_clean_gwas_network_data(
    trait: str,
    max_pages: int = 5,
    pvalue_threshold: float | None = 5e-8,
    use_demo_if_empty: bool = False,
) -> pd.DataFrame:
    """
    Main function used by the Streamlit app.

    trait -> raw GWAS associations -> clean network table.
    """

    raw_df = fetch_associations_by_trait(
        trait=trait,
        max_pages=max_pages,
    )

    clean_df = extract_gwas_network_table(
        associations_df=raw_df,
        disease_name=trait,
    )

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
        pvalue_threshold=5e-8,
        use_demo_if_empty=False,
    )

    print(df.shape)
    print(df.head())
