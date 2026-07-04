# Endometriosis Transcriptomics: Integrative Analysis of Gene Expression, Immune Remodeling, and Fibrotic Pathways

> A reproducible systems biology pipeline integrating transcriptomic datasets to characterize the molecular landscape of endometriosis.

![Python](https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.3+-276DC3?logo=r&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)
![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Transcriptomics-orange)
![Status](https://img.shields.io/badge/Status-Complete-success)

---

## Overview

Endometriosis is a chronic inflammatory disorder characterized by the growth of endometrial-like tissue outside the uterus. Although it affects nearly 10% of reproductive-age women, the molecular mechanisms driving lesion establishment, immune dysregulation, and fibrotic remodeling remain incompletely understood.

This repository presents a **fully reproducible bioinformatics workflow** that integrates two independent Affymetrix microarray datasets to investigate transcriptomic alterations associated with endometriosis.

The analysis combines:

- Differential gene expression analysis
- Cross-study batch correction
- Immune cell signature profiling
- Functional enrichment analysis

to characterize the biological processes underlying the endometriotic microenvironment.

---

## Research Question

> **Which genes, immune cell signatures, and biological pathways characterize the inflammatory and fibrotic microenvironment of ectopic endometrial lesions?**

---

# Highlights

- Integrated **2 independent GEO transcriptomic datasets**
- Analyzed **42 patient samples**
- Harmonized **22,277 shared probes**
- Identified **1,490 differentially expressed genes**
- Quantified **830 immune cell signatures**
- Identified **526 significantly altered immune signatures**
- GO Biological Process and KEGG pathway enrichment
- Fully reproducible **Python + R** workflow

---

# Datasets

| Dataset | Samples | Tissue Types | Platform |
|:---------|---------:|-------------|----------|
| **GSE25628** | 22 | Normal, Eutopic, Ectopic | Affymetrix GPL571 |
| **GSE7305** | 20 | Normal, Diseased | Affymetrix GPL570 |

**Total samples:** **42**

---

# Workflow

```text
                 GEO Datasets
                      │
                      ▼
        Phase 0 ─ Data Loading & QC
                      │
                      ▼
   Phase 1 ─ Batch Correction (ComBat)
                      │
                      ▼
 Phase 2 ─ Differential Expression (limma)
                      │
                      ▼
      Phase 3A ─ Probe Mapping
              + ssGSEA Profiling
                      │
                      ▼
 Phase 3B ─ Immune Signature Statistics
                      │
                      ▼
 Phase 3C ─ GO & KEGG Enrichment
                      │
                      ▼
          Biological Interpretation
```

---

# Methods

| Step | Method |
|:------|:-------|
| Data acquisition | GEOparse |
| Quality control | Boxplots, Density plots, MDS |
| Batch correction | ComBat |
| Differential expression | limma + empirical Bayes |
| Multiple testing | Benjamini–Hochberg FDR |
| Immune profiling | ssGSEA |
| Statistical testing | Kruskal–Wallis + post-hoc analysis |
| Functional enrichment | clusterProfiler |

---

# Repository Structure

```text
endometriosis-transcriptomics/

├── README.md
├── scripts/
│   ├── phase0_data_loading.py
│   ├── phase1_harmonization.py
│   ├── phase2_deg_analysis.R
│   ├── phase3a_probe_to_gene_mapping.R
│   ├── phase3a_immune_deconvolution.py
│   ├── phase3b_differential_immune_infiltration.py
│   └── phase3c_pathway_enrichment.R
│
├── data/
│   ├── raw/
│   └── processed/
│
├── figures/
│
└── results/
```

---

# Pipeline

## Phase 0 — Data Loading & Quality Control

- Download GEO datasets
- Sample annotation
- Log transformation
- Quality assessment
  - Boxplots
  - Density plots
  - Mean expression
  - Multidimensional scaling (MDS)

---

## Phase 1 — Cross-study Harmonization

Datasets were aligned using shared probe identifiers before batch correction with **ComBat**.

**Results**

- 22,277 shared probes retained
- Reduced batch-associated variance
- Preserved biological group separation

---

## Phase 2 — Differential Expression

Differential expression analysis was performed using the **limma** package with empirical Bayes moderation.

### Significant Differentially Expressed Genes

| Comparison | DEGs | Up | Down |
|:------------|----:|---:|-----:|
| Ectopic vs Normal | 1490 | 795 | 695 |
| Diseased vs Normal | 1893 | 991 | 902 |
| Eutopic vs Normal | 216 | 91 | 125 |
| Ectopic vs Eutopic | 145 | 137 | 8 |
| Ectopic vs Diseased | 0 | — | — |

**Filtering criteria**

- FDR < 0.05
- |log₂ Fold Change| > 1

---

## Phase 3A — Immune Cell Signature Profiling

Probe-level expression values were mapped to gene symbols before immune cell enrichment analysis using **ssGSEA**.

### Results

| Metric | Value |
|:-------|------:|
| Shared probes | 22,277 |
| Unique genes | 13,631 |
| Immune signatures | 830 |

---

## Phase 3B — Statistical Analysis

Immune signatures were compared across groups using:

- Kruskal–Wallis tests
- Benjamini–Hochberg correction
- Pairwise post-hoc testing

### Result

**526 immune cell signatures were significantly altered (FDR < 0.05).**

---

## Phase 3C — Functional Enrichment

GO Biological Process and KEGG enrichment analyses were performed using **clusterProfiler**.

Major enriched biological themes included:

- Cell cycle regulation
- Mitotic progression
- Extracellular matrix organization
- Focal adhesion
- Complement activation
- Coagulation cascades

---

# Representative Results

## PCA Before and After Batch Correction


```markdown
![PCA](figures/harmonization/03_comparison_before_after.pdf)
```

---

## Immune Landscape

```markdown
![Immune Heatmap](figures/immune_profiling/authentic_immune_landscape_heatmap.pdf)
```

---

# Biological Insights

The integrated analysis identified several transcriptomic characteristics consistently associated with ectopic endometrial lesions.

### Cell Proliferation

Genes involved in mitosis and cell-cycle progression were significantly upregulated, consistent with increased proliferative activity.

---

### Fibrotic Remodeling

Extracellular matrix organization, focal adhesion pathways, and fibroblast-associated signatures suggest active tissue remodeling and fibrosis.

---

### Chronic Inflammation

Complement activation and coagulation pathways, together with widespread immune alterations, indicate persistent inflammatory signaling.

---

### Immune Dysregulation

Macrophage-associated and immune regulatory signatures were enriched, suggesting a complex immune microenvironment involving both inflammatory activation and immune tolerance.

---

### Disease Progression

The transcriptomic profiles followed the pattern:

```text
Normal
     ↓
Eutopic
     ↓
Ectopic ≈ Diseased
```

No significant differential expression was detected between ectopic and diseased tissues under the selected statistical thresholds, indicating substantial transcriptomic similarity after cross-study harmonization. This observation is consistent with previous studies describing shared molecular characteristics among advanced endometriotic lesions, although additional validation in independent cohorts is warranted.

---

# Limitations

- Bulk microarray data rather than RNA sequencing
- Moderate sample size
- Computational estimation of immune signatures
- Cross-sectional observational design
- Findings represent associations rather than causal mechanisms

---

# Future Directions

Potential extensions include:

- Validation using independent RNA-seq cohorts
- Integration with single-cell RNA sequencing
- Protein–protein interaction network analysis
- Gene regulatory network inference
- Machine learning models for lesion classification
- Comparative analyses with adenomyosis and PCOS

---

# Reproducibility

## Clone the repository

```bash
git clone https://github.com/Theyasna/endometriosis-transcriptomics.git

cd endometriosis-transcriptomics
```

---

## Run the pipeline

```bash
python scripts/phase0_data_loading.py

python scripts/phase1_harmonization.py

Rscript scripts/phase2_deg_analysis.R

Rscript scripts/phase3a_probe_to_gene_mapping.R

python scripts/phase3a_immune_deconvolution.py

python scripts/phase3b_differential_immune_infiltration.py

Rscript scripts/phase3c_pathway_enrichment.R
```

---

# Tools & Packages

### Python

- pandas
- numpy
- scipy
- scikit-learn
- matplotlib
- GEOparse
- GSEApy
- pyComBat
- statsmodels

### R

- limma
- clusterProfiler
- ggplot2
- AnnotationDbi
- org.Hs.eg.db
- hgu133plus2.db



---

# Author

**Yasna Dehestan**

*Bioinformatics • Computational Biology • Transcriptomics • Women's Health Genomics*

📧 theyasnad@gmail.com

🔗 https://github.com/Theyasna

---

# License

Released under the **MIT License**.
