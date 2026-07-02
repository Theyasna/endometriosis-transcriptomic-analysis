# Endometriosis Transcriptomics: A Systems Biology Pipeline

## Research Question

Which genes, cell types, and biological pathways drive the inflammatory and fibrotic microenvironment of ectopic endometrial lesions?

---

## Quick Summary

This repository contains a **complete, reproducible pipeline** for analyzing transcriptomic data in endometriosis. We integrated two independent Affymetrix microarray datasets (GSE25628, GSE7305), performed batch correction, identified 1,490 differentially expressed genes in ectopic tissue, deconvoluted the immune landscape using MSigDB signatures, and identified enriched biological pathways including cell cycle dysregulation and chronic inflammation.

**Key finding:** Ectopic tissue is characterized by elevated adventitial fibroblast signatures, complement cascade activation, and dysregulated mitotic control—suggesting that lesion survival depends on active fibrotic remodeling and immune tolerance.

---

## Pipeline Overview

```
Phase 0: Load data        (Python)
         ↓
Phase 1: Harmonize        (Python + ComBat)
         ↓
Phase 2: DEG Analysis     (R + limma)
         ↓
Phase 3a: Gene Mapping    (R) → Immune Deconvolution (Python)
         ↓
Phase 3b: Immune Stats    (Python + Kruskal-Wallis)
         ↓
Phase 3c: Pathways        (R + clusterProfiler)
```

---

## What Each Phase Does

### Phase 0: Data Loading & QC
**Language:** Python  
**Script:** `scripts/phase0_data_loading.py`  
**Runtime:** 30-45 min (includes GEO download)

Downloads GSE25628 (n=22) and GSE7305 (n=20) from NCBI GEO, applies QC checks, detects and corrects log scale differences.

**Outputs:**
- `data/processed/GSE25628_expression_matrix.csv` (22,277 probes × 22 samples)
- `data/processed/GSE7305_expression_matrix.csv` (54,675 probes × 20 samples)  
- Sample info CSVs + QC plots (boxplot, density, MDS)

---

### Phase 1: Harmonization & Batch Correction
**Language:** Python  
**Script:** `scripts/phase1_harmonization.py`  
**Runtime:** 5-10 min

Aligns probes across datasets (22,277 common probes = 100% retention), applies **ComBat batch correction** while preserving biological signal via group covariates, validates effect via PCA.

**Key Results:**
- Batch effect removed: PC1 variance 46.3% → 32.2% ✓
- Merged matrix: 22,277 probes × 42 samples  
- Final composition: 16 normal, 8 ectopic, 8 eutopic, 10 diseased

**Outputs:**
- `data/processed/harmonized_expression_matrix.csv`
- `data/processed/harmonized_sample_info.csv`
- `figures/harmonization/02_PCA_before_after_correction.pdf`

---

### Phase 2: Differential Expression Analysis
**Language:** R  
**Script:** `scripts/phase2_deg_analysis.R`  
**Runtime:** 10-15 min

Fits linear models using **limma** (gold standard for microarray), applies empirical Bayes moderation, computes 4 contrasts, applies Benjamini-Hochberg FDR correction (α=0.05).

**Contrasts & Results:**
| Comparison | DEGs | Up | Down |
|-----------|------|-----|------|
| Ectopic vs. Normal | **1,490** | 795 | 695 |
| Diseased vs. Normal | **1,893** | 991 | 902 |
| Eutopic vs. Normal | **216** | 91 | 125 |
| Ectopic vs. Eutopic | **145** | 137 | 8 |

**Outputs:**
- Volcano plots (4×) + MA plots + heatmap  
- DEG tables (full results + significant genes)

---

### Phase 3a: Gene Mapping & Immune Deconvolution
**Languages:** R (probe mapping) + Python (ssGSEA)  
**Scripts:** 
- `scripts/phase3a_probe_to_gene_mapping.R` (5 min)
- `scripts/phase3a_immune_deconvolution.py` (10-15 min)

**Step 1 (R):** Maps 22,277 Affymetrix probes → 13,631 unique HGNC genes using `hgu133plus2.db`. Collapses duplicate probes by keeping highest-variance transcript.

**Step 2 (Python):** Applies single-sample GSEA (ssGSEA) using MSigDB C8 cell-type signatures (830 signatures). Generates immune landscape heatmap.

**Key Finding:** 526 of 830 immune signatures show significant group differences (FDR < 0.05), suggesting major immune remodeling in endometriosis.

**Outputs:**
- `results/immune_profiling/authentic_immune_enrichment_scores.csv`
- `figures/immune_profiling/authentic_immune_landscape_heatmap.pdf`

---

### Phase 3b: Statistical Validation of Immune Changes
**Language:** Python  
**Script:** `scripts/phase3b_differential_immune_infiltration.py`  
**Runtime:** 5 min

Applies **Kruskal-Wallis omnibus test** (appropriate for non-normal enrichment scores) + Mann-Whitney U pairwise tests with FDR correction.

**Key Results:**
- 526 signatures significant at FDR < 0.05  
- Top 6 plotted with individual sample data points  
- Post-hoc pairwise comparisons included

**Outputs:**
- `results/immune_profiling/differential_immune_infiltration_stats.csv`
- `figures/immune_profiling/differential_immune_boxplots.pdf`

---

### Phase 3c: Pathway Enrichment Analysis
**Language:** R  
**Script:** `scripts/phase3c_pathway_enrichment.R`  
**Runtime:** 10 min

Maps DEGs to Entrez IDs, performs Over-Representation Analysis (ORA) for:
- **GO Biological Processes** (clusterProfiler)  
- **KEGG Pathways** (clusterProfiler)  
- FDR correction (α=0.05)

**Top Enriched Pathways:**
- Cell cycle progression (p=1e-8) — Ectopic lesions proliferate aggressively
- Complement cascades (p=1e-7) — Chronic inflammatory niche
- ECM-receptor interaction (p=1e-6) — Fibrotic tissue remodeling
- Focal adhesion (p=1e-5) — Abnormal cell implantation

**Outputs:**
- `figures/pathway_enrichment/01_go_enrichment_dotplot.pdf`  
- `figures/pathway_enrichment/02_kegg_enrichment_barplot.pdf`

---

## Installation & Setup

### Requirements

**Python 3.10+**
```bash
pip install pandas numpy scipy scikit-learn matplotlib seaborn GEOparse gseapy
```

**R 4.0+**
```R
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", 
                        "hgu133plus2.db", "AnnotationDbi"))

install.packages(c("ggplot2", "ggrepel", "dplyr"))
```

### Project Structure
```
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
├── data/
│   ├── raw/              ← GEO downloads (auto-created)
│   └── processed/        ← Phase output (auto-created)
├── figures/              ← Publication-quality PDFs
└── results/              ← CSV result tables
```

---

## Quick Start

```bash
# Clone repo
git clone https://github.com/Theyasna/endometriosis-transcriptomics.git
cd endometriosis-transcriptomics

# Phase 0: Download & QC (30-45 min, ~2 GB)
python scripts/phase0_data_loading.py

# Phase 1: Harmonization (5-10 min)
python scripts/phase1_harmonization.py

# Phase 2: DEG Analysis (10-15 min)
Rscript scripts/phase2_deg_analysis.R

# Phase 3a: Gene Mapping + Immune Deconvolution (15-20 min)
Rscript scripts/phase3a_probe_to_gene_mapping.R
python scripts/phase3a_immune_deconvolution.py

# Phase 3b: Immune Statistics (5 min)
python scripts/phase3b_differential_immune_infiltration.py

# Phase 3c: Pathway Enrichment (10 min)
Rscript scripts/phase3c_pathway_enrichment.R
```

**Total runtime:** ~90 minutes (first run includes data download)

---

## Key Findings

### Molecular Signature of Ectopic Tissue

1. **Proliferative Gene Expression** 
   - Upregulation of cell cycle, mitosis pathways  
   - Enables survival in hostile environment outside uterus

2. **Fibrotic Remodeling** 
   - Elevated adventitial fibroblast signatures  
   - Active ECM remodeling

3. **Chronic Inflammation** 
   - Complement cascade & coagulation enriched
   - 526/830 immune signatures altered

4. **Immune Tolerance** 
   - Both M1 and M2 macrophages elevated
   - Suggests lesion immune equilibrium

---

## Limitations & Caveats 

### 1. Affymetrix Probes, Not RNA-seq 
- Probe-to-gene mapping adds inference step
- All data from microarray (no RNA-seq validation)


### 2. Small Sample Size 
- 42 total samples; ectopic/eutopic groups n=8 each


### 3. Immune Deconvolution Limits 
- ssGSEA estimates relative, not absolute, abundance
- MSigDB signatures not validated on endometrial tissue


### 4. No Severity Correlation 
- Cannot link genes to disease stage
- All ectopic samples treated equally

---


## Contact

**Yasna Dehestan**  
Bioinformatics | Women's Health Genomics  
Email: theyasnad@gmail.com  
GitHub: [@Theyasna](https://github.com/Theyasna)

---

## License

MIT License
