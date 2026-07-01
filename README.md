## Project Overview

This repository contains a robust computational pipeline for the integrative analysis of transcriptomic data in endometriosis. The project performs **end-to-end processing**, from raw microarray harmonization to immune microenvironment deconvolution and functional pathway enrichment. The goal is to identify the transcriptomic signatures and cellular populations driving the inflammatory and remodeling processes in ectopic endometrial lesions.

---

## Workflow

The pipeline is divided into **three modular phases**:

### **1. Data Harmonization**
Integration of independent datasets (GSE7305 & GSE25628). PCA-based batch correction was implemented to mitigate technical artifacts while preserving biological variance across clinical cohorts (**Normal, Eutopic, Ectopic, Diseased**).

### **2. Differential Expression Analysis**
Identification of significant transcriptional shifts. Probes were mapped to Entrez IDs to facilitate functional interpretation.

### **3. Immune Microenvironment Deconvolution**
Application of ssGSEA using MSigDB C8 (Cell-Type) signatures to identify specific immune and stromal cell infiltration patterns associated with disease progression.

### **4. Functional Enrichment Analysis**
Over-Representation Analysis (ORA) using Gene Ontology (GO: Biological Processes) and KEGG pathways to link gene lists to underlying mechanisms, such as **cell cycle dysregulation** and **ECM-receptor interaction**.

---

## Key Findings
* **Cellular Drivers**: Analysis revealed a significant enrichment of *Adventitial Fibroblast* and *smooth muscle cell* signatures in ectopic tissue, suggesting active tissue remodeling. 
  ![Immune Landscape](figures/immune_profiling/differential_immune_boxplots.pdf)

* **Inflammatory Signature**: Enriched *Complement and coagulation cascades* (KEGG) correlate with the inflammatory niche of endometriosis lesions. 

* **Proliferative Capacity**: GO analysis highlighted persistent upregulation of *mitotic nuclear division* and *cell cycle progression* pathways, a molecular hallmark of ectopic lesion survival.
  ![Cell Cycle Pathways](figures/pathway_enrichment/01_go_enrichment_dotplot.pdf)



R (v4.x): limma, AnnotationDbi, clusterProfiler, org.Hs.eg.db, enrichplot

Python (v3.x): gseapy, statsmodels, scipy, seaborn, pandas
