# Endometriosis Transcriptomics: Per-Dataset Analysis of Gene Expression, Immune Signatures, and Functional Pathways

> A reproducible bioinformatics pipeline integrating two independent transcriptomic datasets to characterize the molecular landscape of endometriosis, with cross-cohort replication of lesion-vs-control findings.

![Python](https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white)
![R](https://img.shields.io/badge/R-4.3+-276DC3?logo=r&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-green)
![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Transcriptomics-orange)

---

## Highlights

- 2 independent GEO datasets, 2 different Affymetrix platforms
- 42 samples across ectopic lesion, eutopic, and control tissue
- Per-dataset statistical design to avoid cross-platform confounding
- 3,604 DEGs identified in ectopic endometrial lesions
- Lesion-vs-control DEGs **replicated across both cohorts** (608-gene overlap)
- 4,872 immune signatures profiled (MSigDB C7 ImmuneSigDB)
- Fully reproducible Python + R workflow

---

## Research Question

> Which genes, immune signatures, and biological pathways characterize endometriotic lesion tissue relative to control endometrium, and do these signatures replicate across independent patient cohorts?

---

## Datasets

| Dataset | Samples | Groups | Platform |
|:---|---:|:---|:---|
| **GSE25628** | 22 | Ectopic lesion (8), Eutopic (8), Normal (6) | Affymetrix GPL571 |
| **GSE7305** | 20 | Ovarian lesion (10), Matched control (10) | Affymetrix GPL570 |

**Total samples:** 42

Both datasets contain a **lesion vs. control-endometrium comparison**, run on independent cohorts and different platforms — the basis for the cross-cohort replication in this project. Group definitions were verified against the original GEO records and source publications (Crispi et al. 2013 for GSE25628; Hever et al. 2007 for GSE7305). See [docs/dataset_notes.md](docs/dataset_notes.md) for full detail, including the important note that GSE7305's "diseased" group is ovarian endometriosis **lesion** tissue, not eutopic endometrium.

---

## Workflow

```
GEO Datasets → Phase 0 (QC) → Phase 1 (Harmonization)
    → Phase 2 (Per-Dataset DEG) → Phase 3A (Gene Mapping + ssGSEA)
    → Phase 3B (Immune Statistics) → Phase 3C (GO + KEGG Enrichment)
    → Cross-Cohort Interpretation
```

---

## Methods

| Step | Method |
|:---|:---|
| Data acquisition | GEOparse |
| Quality control | Boxplots, density plots, MDS |
| Cross-study handling | Per-dataset analysis (see [docs/methodology.md](docs/methodology.md)) |
| Differential expression | limma + empirical Bayes |
| Multiple testing | Benjamini–Hochberg FDR |
| Immune profiling | ssGSEA, MSigDB C7 ImmuneSigDB |
| Statistical testing | Kruskal–Wallis + Mann-Whitney post-hoc |
| Functional enrichment | clusterProfiler (GO-BP, GO-CC, KEGG) |

The two datasets use different Affymetrix platforms and cannot be pooled at the sample level without confounding platform with biology. Analyses were therefore run independently within each dataset, and the two lesion-vs-control results were compared at the gene-list level for replication. See [docs/methodology.md](docs/methodology.md).

---

## Results

### Differential Expression

Per-dataset analysis (FDR < 0.05, |log₂FC| > 1):

| Comparison | Dataset | DEGs | Up | Down |
|:---|:---|---:|---:|---:|
| Ectopic lesion vs Normal | GSE25628 | 3,604 | 1,150 | 1,742 |
| Eutopic vs Normal | GSE25628 | 2,015 | 521 | 1,112 |
| Ectopic vs Eutopic | GSE25628 | 21 | 19 | 0 |
| Ovarian lesion vs Control | GSE7305 | 1,603 | 729 | 552 |

**Cross-cohort replication:** the two lesion-vs-control comparisons (ectopic-vs-normal in GSE25628, ovarian-lesion-vs-control in GSE7305) share **608 DEGs (37.9% of the smaller set)** — replication across two independent cohorts, two platforms, and two slightly different control definitions.

*Not computed, by design:* comparisons pairing GSE25628-only groups against GSE7305-only groups. These would confound biology with platform/cohort and cannot be validly computed from these two datasets.

---

### Immune Signature Profiling

ssGSEA against 4,872 MSigDB C7 ImmuneSigDB gene sets, analyzed per dataset:

| Dataset | Significant signatures (FDR < 0.05) |
|:---|:---|
| GSE25628 | 1,716 / 4,872 (35.2%) |
| GSE7305 | 3,542 / 4,872 (72.7%) |

The dominant immune axis in both datasets tracked group membership (PC1 = 60.9% and 75.6% of variance in GSE25628 and GSE7305 respectively), consistent with a coordinated immune/inflammatory shift rather than independent per-signature effects. See [docs/validation.md](docs/validation.md) for the PCA diagnostic and its caveats.

---

### Pathway Enrichment

| Comparison | GO BP | GO CC | KEGG (total / flagged) |
|:---|---:|---:|---:|
| GSE25628 Ectopic vs Normal | 883 | 124 | 94 / 18 |
| GSE25628 Eutopic vs Normal | 219 | 64 | 48 / 4 |
| GSE25628 Ectopic vs Eutopic | 0 | 3 | 0 / — |
| GSE7305 Lesion vs Control | 1,060 | 117 | 52 / 8 |

Several infection-associated KEGG pathways (e.g. HTLV-1 infection, Malaria, *S. aureus* infection) were flagged as likely gene-overlap artifacts arising from shared complement and inflammatory gene membership rather than evidence of infectious processes. A `likely_gene_overlap_artifact` column is included in each KEGG result table. See [docs/validation.md](docs/validation.md) for the cross-dataset recurrence evidence.

---

## Biological Insights

**Cell cycle and proliferation:** Mitotic progression genes were significantly upregulated in ectopic lesions, consistent with enhanced proliferative activity.

**ECM remodeling and fibrosis:** Focal adhesion, ECM-receptor interaction, and integrin signaling enriched across contrasts, suggesting active tissue remodeling.

**Immune dysregulation:** A coordinated immune signature shift, together with complement and coagulation pathway enrichment, was observed in both cohorts independently — the most consistent cross-dataset theme.

**Eutopic endometrium is transcriptomically altered:** 2,015 DEGs vs. healthy endometrium suggest eutopic tissue from patients is not equivalent to normal endometrium.

**Lesion and eutopic tissue are molecularly similar within GSE25628:** only 21 DEGs between ectopic and eutopic tissue, suggesting a largely shared transcriptional program.

---

## What This Project Does Not Support

- Any comparison pairing GSE25628-only groups against GSE7305-only groups (would confound biology with platform/cohort)
- Interpretation of immune-signature hit rates as thousands of independent discoveries (a single dominant axis drives the signal in each dataset)
- Menstrual-cycle-controlled claims for GSE7305 beyond what the paired within-patient design supports (see docs)

---

## Limitations

- Bulk microarray data rather than RNA-seq
- Moderate per-group sample sizes (n = 6–10)
- Different platforms prevent sample-level pooling; replication is at the gene-list level
- GSE7305 uses a paired within-patient design that the current analysis treats as unpaired (a conservative choice; see [docs/methodology.md](docs/methodology.md))
- Immune signature results reflect correlated gene sets driven by a small number of dominant axes
- ssGSEA findings are computational estimates, not validated experimentally

---

## Repository Structure

```
endometriosis-transcriptomics/
├── README.md
├── docs/
│   ├── methodology.md
│   ├── dataset_notes.md
│   └── validation.md
├── scripts/
│   ├── phase0_data_loading.py
│   ├── phase1_harmonization.py
│   ├── phase2_deg_analysis.R
│   ├── phase3a1_probe_to_gene_mapping.R
│   ├── phase3a2_immune_deconvolution.py
│   ├── phase3b_differential_immune_infiltration.py
│   └── phase3c_pathway_enrichment.R
├── data/
│   ├── raw/
│   └── processed/
├── figures/
└── results/
```

---

## Reproducibility

```bash
git clone https://github.com/Theyasna/endometriosis-transcriptomics.git
cd endometriosis-transcriptomics

python scripts/phase0_data_loading.py
python scripts/phase1_harmonization.py
Rscript scripts/phase2_deg_analysis.R
Rscript scripts/phase3a1_probe_to_gene_mapping.R
python scripts/phase3a2_immune_deconvolution.py
python scripts/phase3b_differential_immune_infiltration.py
Rscript scripts/phase3c_pathway_enrichment.R
```

---

## Tools & Packages

**Python:** pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, GEOparse, GSEApy, statsmodels

**R:** limma, clusterProfiler, ggplot2, ggrepel, AnnotationDbi, org.Hs.eg.db, hgu133plus2.db

---

## Author

**Yasna Dehestan**
*Bioinformatics · Computational Biology · Transcriptomics · Women's Health Genomics*

📧 theyasnad@gmail.com · 🔗 [github.com/Theyasna](https://github.com/Theyasna)

---

## License

MIT License
