````markdown
# Endometriosis Transcriptomics: Integrative Analysis of Gene Expression, Immune Remodeling, and Fibrotic Pathways

> A fully reproducible systems biology pipeline integrating multiple transcriptomic datasets to investigate the molecular mechanisms underlying endometriosis.

![Python](https://img.shields.io/badge/Python-3.10+-blue)
![R](https://img.shields.io/badge/R-4.3+-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-orange)
![Transcriptomics](https://img.shields.io/badge/Focus-Transcriptomics-red)
![Reproducible](https://img.shields.io/badge/Reproducible-Research-success)

---

# Overview

Endometriosis is a chronic inflammatory disease characterized by the growth of endometrial-like tissue outside the uterus, affecting approximately 10% of reproductive-age women worldwide. Although its clinical manifestations are well recognized, the molecular mechanisms that promote lesion establishment, immune evasion, and fibrotic remodeling remain incompletely understood.

This project presents a **fully reproducible transcriptomic analysis pipeline** integrating two independent GEO microarray datasets to investigate molecular alterations associated with endometriosis.

Rather than focusing solely on differential gene expression, the pipeline combines:

- Differential expression analysis
- Cross-study batch correction
- Immune cell signature profiling
- Functional enrichment analysis

to characterize the biological processes shaping the endometriotic microenvironment.

---

# Research Question

> **Which genes, immune cell programs, and biological pathways characterize the inflammatory and fibrotic microenvironment of ectopic endometrial lesions?**

---

# Project Highlights

- Integrated **2 independent GEO microarray datasets**
- Processed **42 patient samples**
- Harmonized **22,277 shared probes** using ComBat
- Identified **1,490 differentially expressed genes** in ectopic lesions
- Quantified **830 immune cell signatures**
- Detected **526 significantly altered immune programs**
- Performed GO Biological Process and KEGG enrichment analyses
- Fully reproducible **Python + R** workflow

---

# Datasets

| Dataset | Samples | Tissue Types | Platform |
|----------|---------:|-------------|----------|
| GSE25628 | 22 | Normal, Eutopic, Ectopic | Affymetrix |
| GSE7305 | 20 | Normal, Diseased | Affymetrix |

**Total samples:** 42

---

# Analysis Pipeline

```text
Raw GEO Data
      │
      ▼
Phase 0
Quality Control
      │
      ▼
Phase 1
Cross-study Harmonization
(ComBat)
      │
      ▼
Phase 2
Differential Expression
(limma)
      │
      ▼
Phase 3A
Probe Mapping
+
ssGSEA Immune Profiling
      │
      ▼
Phase 3B
Statistical Analysis
      │
      ▼
Phase 3C
GO / KEGG Enrichment
      │
      ▼
Biological Interpretation
```

---

# Methods at a Glance

| Step | Method |
|------|--------|
| Data download | GEOparse |
| Quality control | Boxplots, density plots, MDS |
| Batch correction | ComBat |
| Differential expression | limma + empirical Bayes |
| Multiple testing | Benjamini–Hochberg FDR |
| Immune profiling | ssGSEA |
| Statistical testing | Kruskal–Wallis + post-hoc tests |
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

# Pipeline Details

## Phase 0 — Data Loading & Quality Control

- Download GEO datasets
- Log transformation
- Sample annotation
- Quality assessment using:
  - Boxplots
  - Density plots
  - Mean expression
  - MDS visualization

---

## Phase 1 — Cross-study Harmonization

- Identification of common probes
- Alignment of datasets
- ComBat batch correction
- PCA validation before and after correction

**Result**

- 22,277 shared probes retained
- Batch-associated variance substantially reduced while preserving biological group differences

---

## Phase 2 — Differential Expression

Differential expression analysis was performed using **limma** with empirical Bayes moderation.

### Significant Differentially Expressed Genes

| Comparison | DEGs | Up | Down |
|------------|------:|----:|------:|
| Ectopic vs Normal | 1,490 | 795 | 695 |
| Diseased vs Normal | 1,893 | 991 | 902 |
| Eutopic vs Normal | 216 | 91 | 125 |
| Ectopic vs Eutopic | 145 | 137 | 8 |
| Ectopic vs Diseased | 0 | — | — |

Criteria:

- FDR < 0.05
- |log₂ Fold Change| > 1

---

## Phase 3A — Immune Cell Signature Profiling

Probe-level expression was converted to gene-level expression before performing ssGSEA using MSigDB immune cell signatures.

Results:

- 22,277 probes
- 13,631 unique genes
- 830 immune signatures quantified

---

## Phase 3B — Statistical Analysis

Immune signatures were compared across groups using:

- Kruskal–Wallis tests
- Multiple testing correction
- Pairwise post-hoc comparisons

Result:

**526 immune signatures were significantly altered (FDR < 0.05).**

---

## Phase 3C — Functional Enrichment

GO Biological Process and KEGG pathway enrichment analyses were performed using clusterProfiler.

Major enriched biological themes included:

- Cell cycle regulation
- Mitotic progression
- ECM organization
- Focal adhesion
- Complement activation
- Coagulation pathways

---

# Representative Results

## PCA Before and After Batch Correction

```markdown
![PCA](figures/harmonization/03_comparison_before_after.pdf)
```

## Immune Cell Signature Heatmap

```markdown
![Immune Heatmap](figures/immune_profiling/authentic_immune_landscape_heatmap.pdf)
```

---

# Biological Insights

The integrated analysis revealed several consistent molecular features associated with ectopic endometrial lesions.

## Increased Cellular Proliferation

Genes involved in mitosis and cell-cycle progression were strongly upregulated, consistent with enhanced proliferative activity within ectopic tissue.

---

## Fibrotic Remodeling

Extracellular matrix organization, focal adhesion pathways, and fibroblast-associated signatures were enriched, supporting extensive tissue remodeling.

---

## Chronic Inflammation

Complement activation, coagulation pathways, and widespread immune-cell alterations indicate persistent inflammatory activity.

---

## Immune Dysregulation

Multiple macrophage-associated and immune regulatory signatures were enriched, suggesting a complex immune microenvironment involving both inflammatory activation and immune tolerance.

---

## Disease Progression

The transcriptomic analyses suggest a molecular progression:

```text
Normal
     ↓
Eutopic
     ↓
Ectopic ≈ Diseased
```

No significant differential expression was detected between ectopic and diseased tissues under the selected statistical thresholds, indicating substantial transcriptomic similarity after cross-study harmonization. This observation is consistent with previous reports describing shared molecular characteristics among advanced endometriotic lesions, although further validation in independent cohorts is warranted.

---

# Limitations

- Based on bulk microarray data rather than RNA sequencing
- Moderate sample size
- Immune infiltration estimated computationally using ssGSEA
- Cross-sectional observational study
- Findings represent associations rather than causal mechanisms

---

# Future Directions

Potential extensions of this project include:

- Validation using independent RNA-seq cohorts
- Integration with single-cell RNA sequencing
- Protein-protein interaction network analysis
- Gene regulatory network inference
- Machine learning models for lesion classification
- Comparison with adenomyosis and PCOS transcriptomic datasets

---

# Reproducibility

Clone the repository:

```bash
git clone https://github.com/Theyasna/endometriosis-transcriptomics.git

cd endometriosis-transcriptomics
```

Run the pipeline sequentially:

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

# Author

**Yasna Dehestan**

Bioinformatics • Transcriptomics • Systems Biology • Women's Health Genomics

GitHub: https://github.com/Theyasna

Email: theyasnad@gmail.com

---

# License

This project is released under the MIT License.
````
