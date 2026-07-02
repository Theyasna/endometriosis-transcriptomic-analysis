# Endometriosis Transcriptomics: A Systems Biology Pipeline

## Research Question
Which genes, immune cell types, and biological pathways drive the inflammatory and fibrotic microenvironment of ectopic endometrial lesions?

---

## Project Summary
This repository contains a **complete, reproducible bioinformatics pipeline** integrating two independent Affymetrix microarray datasets (GSE25628 and GSE7305, total n=42 samples). 


---

## Pipeline Architecture


```
Phase 0: Data Loading & QC (Python)
↓
Phase 1: Batch Correction & Harmonization (Python + ComBat)
↓
Phase 2: Differential Expression (R + limma)
↓
Phase 3a: Probe-to-Gene Mapping (R) + ssGSEA Immune Deconvolution (Python)
↓
Phase 3b: Statistical Testing of Immune Changes (Python)
↓
Phase 3c: GO & KEGG Pathway Enrichment (R + clusterProfiler)
```


---

## Phase Details

### Phase 0: Data Loading & Quality Control
- **Script**: `scripts/phase0_data_loading_fixed.py`
- Downloads and processes GSE25628 (22 samples) and GSE7305 (20 samples)
- Performs log transformation, group assignment, and comprehensive QC (boxplots, density, MDS, mean expression)

### Phase 1: Harmonization & Batch Correction
- **Script**: `scripts/phase1_harmonization_fixed.py`
- Aligns 22,277 common probes (100% retention from GSE25628)
- Applies **ComBat** batch correction while preserving group effects
- Validates correction with PCA (PC1 variance reduced from 46.3% to 32.2%)

### Phase 2: Differential Expression Analysis
- **Script**: `scripts/phase2_deg_fixed.R`
- Uses **limma + eBayes** for robust statistical testing
- Four contrasts with Benjamini-Hochberg FDR correction

**Significant DEGs (FDR < 0.05, |logFC| > 1)**:

| Comparison              | Total | Up   | Down |
|-------------------------|-------|------|------|
| Ectopic vs Normal       | 1,490 | 795  | 695  |
| Diseased vs Normal      | 1,893 | 991  | 902  |
| Eutopic vs Normal       | 216   | 91   | 125  |
| Ectopic vs Eutopic      | 145   | 137  | 8    |

### Phase 3: Immune Profiling & Pathway Analysis
- **3a**: Probe-to-gene mapping (22,277 probes → 13,039 unique genes) + ssGSEA using MSigDB C8 cell-type signatures
- **3b**: Kruskal-Wallis + post-hoc testing → **526 significant immune signatures** (FDR < 0.05)
- **3c**: GO Biological Process and KEGG pathway enrichment

**Notable Enriched Pathways**:
- Cell cycle & mitotic regulation
- ECM-receptor interaction & focal adhesion (fibrosis)
- Complement & coagulation cascades (inflammation)

---

## Installation & Running the Pipeline

### Requirements
**Python**:
```bash
pip install pandas numpy scipy scikit-learn matplotlib seaborn GEOparse gseapy pycombat
```
**R (Bioconductor)**:
```
RBiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db", "hgu133plus2.db", "AnnotationDbi"))
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
   -Significant upregulation of cell cycle and mitotic pathways
   -Supports aggressive growth and survival of ectopic lesions outside the uterus

2. **Fibrotic Remodeling** 
   - Elevated signatures of adventitial fibroblasts and smooth muscle cells
   - Active extracellular matrix (ECM) remodeling and tissue restructuring

3. **Chronic Inflammation** 
   - Enrichment of complement cascade and coagulation pathways
   - 526 out of 830 immune cell-type signatures show significant alteration (FDR < 0.05)

4. **Immune Tolerance** 
   - Evidence of both pro- and anti-inflammatory signals (M1/M2 macrophage signatures)
   - Suggests establishment of immune tolerance that allows lesion persistence

---

## Limitations & Caveats 

1.Microarray-based (probe-level data)
2.Moderate sample sizes per subgroup (n=8 for ectopic/eutopic)
3.Computational immune estimates (ssGSEA)
4.Observational study — requires experimental validation

---


## Contact

**Yasna Dehestan**  
Bioinformatics | Women's Health Genomics  
Email: theyasnad@gmail.com  
GitHub: [@Theyasna](https://github.com/Theyasna)

---

## License

MIT License
