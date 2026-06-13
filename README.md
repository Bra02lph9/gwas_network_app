# GWAS 3D Network Explorer

An interactive bioinformatics application for exploring Genome-Wide Association Study (GWAS) results through dynamic 3D biological networks.

The application transforms GWAS associations into an intuitive graph structure linking:

- Genetic variants (SNPs)
- Genes
- Diseases / Traits

allowing researchers and students to visually explore genetic architectures and identify important biological entities.

---

# Features

## Data Sources

### Upload Your Own Data

Upload CSV or TSV files containing:

| Column | Description |
|----------|-------------|
| snp | Variant ID (rsID) |
| gene | Associated gene |
| disease | Disease or trait |
| p_value | Association p-value |

Example:

```csv
snp,gene,disease,p_value
rs7412,APOE,Coronary Artery Disease,1e-12
rs10455872,LPA,Coronary Artery Disease,3e-30
rs7903146,TCF7L2,Type 2 Diabetes,1e-50
```

### GWAS Catalog Integration

Fetch association data directly from the GWAS Catalog API.

Supported examples:

- Coronary Artery Disease
- Type 2 Diabetes Mellitus
- Body Mass Index
- Alzheimer's Disease

---

# Network Architecture

The application builds a biological graph:

```text
Variant (SNP)
      │
      ▼
Gene
      │
      ▼
Disease
```

Example:

```text
rs429358
     │
     ▼
APOE
     │
     ▼
Coronary Artery Disease
```

---

# Interactive 3D Visualization

Features:

- Fully interactive 3D network
- Zoom
- Rotate
- Pan
- Hover tooltips
- Dynamic node sizing
- Node coloring by biological type

Node colors:

| Node Type | Color |
|------------|--------|
| Variant | Blue |
| Gene | Green |
| Disease | Red |

---

# Network Analysis

Current analytical features include:

### Graph Summary

Displays:

- Total nodes
- Total edges
- Number of variants
- Number of genes
- Number of diseases

### Hub Gene Detection

Identifies highly connected genes using node degree.

Example:

| Gene | Degree |
|--------|----------|
| APOE | 18 |
| LPA | 15 |
| SORT1 | 12 |

These genes may represent important biological hubs.

---

# Filtering

Users can filter associations by:

### P-value Threshold

Default:

```text
5 × 10⁻⁸
```

Standard genome-wide significance threshold.

### Maximum Associations

Limit network size for performance and visualization clarity.


---

# Installation

Clone repository:

```bash
git clone link._.

cd gwas-network-explorer
```

Create virtual environment:

```bash
python -m venv .venv
```

Activate environment:

Windows:

```bash
.venv\Scripts\activate
```

Linux / macOS:

```bash
source .venv/bin/activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

# Run Application

```bash
streamlit run app.py
```

Application opens automatically in:

```text
http://localhost:8501
```

---

# Requirements

Main dependencies:

```text
streamlit
pandas
networkx
plotly
requests
urllib3
```

---

# Future Development

Planned features:

## Biological Analysis

- Hub SNP detection
- Disease comparison
- Shared gene discovery
- Network communities
- Centrality metrics

## Functional Genomics

- STRING Protein Interaction Integration
- Pathway Enrichment Analysis
- Gene Ontology Enrichment
- Reactome Integration

## AI Integration

- Automatic biological interpretation
- Gene function summaries
- Disease mechanism explanations
- Network insight generation

## Export Options

- CSV
- Excel
- GraphML
- HTML Interactive Networks
- PDF Reports

---

# Scientific Goal

The objective of this project is to transform GWAS association data into interpretable biological networks, helping researchers identify:

- Important disease-associated genes
- Shared genetic mechanisms
- Highly connected hub genes
- Potential biological targets for further investigation

---

Bioinformatics | Genomics | Network Biology | AI Applications in Life Sciences