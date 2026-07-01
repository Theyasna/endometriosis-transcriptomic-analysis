#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 3a: AUTHENTIC IMMUNE CELL DECONVOLUTION (ssGSEA) via MSigDB C8
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy import Msigdb

# ============================================================================
# 1. CONFIGURATION & PATHS
# ============================================================================
PROJECT_DIR = r"C:\Users\Yasna\OneDrive\Belgeler\endometriosis-transcriptomic-analysis"

GENE_EXPR_PATH = os.path.join(PROJECT_DIR, "data", "processed", "gene_expression_matrix.csv")
SAMPLE_INFO_PATH = os.path.join(PROJECT_DIR, "data", "processed", "harmonized_sample_info.csv")

RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "immune_profiling")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "immune_profiling")

# Ensure output directories exist
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================================
# 2. LOAD & SYNCHRONIZE DATA
# ============================================================================
print("Loading gene-level expression matrix...")
df_genes = pd.read_csv(GENE_EXPR_PATH, index_col=0)
df_sample = pd.read_csv(SAMPLE_INFO_PATH, index_col=0)

# Synchronize matrix columns to match the sample metadata rows exactly
df_genes = df_genes[df_sample.index]
print(f"✓ Data Synced: {df_genes.shape[0]} unique genes across {df_genes.shape[1]} samples.")

# ============================================================================
# 3. FETCH MSigDB C8 CELL SIGNATURES
# ============================================================================
print("\nConnecting to MSigDB to fetch C8 Cell-Type Signatures...")
try:
    msig = Msigdb()
    c8_genesets = msig.get_gmt(category='c8')
    print(f"✓ Successfully loaded {len(c8_genesets)} cell-type gene sets from MSigDB.")
except Exception as e:
    print(f"⚠ Online MSigDB download standard category failed: {e}")
    print("Attempting alternative fallback category name string...")
    msig = Msigdb()
    c8_genesets = msig.get_gmt(category='c8.all')
    print(f"✓ Loaded {len(c8_genesets)} fallback gene sets successfully.")

# ============================================================================
# 4. RUN ssGSEA DECONVOLUTION
# ============================================================================
print("\nRunning single-sample GSEA (ssGSEA)...")
ssgsea_run = gp.ssgsea(
    data=df_genes,
    gene_sets=c8_genesets,          # Pass downloaded dictionary directly
    outdir=None,                    # Keep results in memory for processing
    sample_norm_method='rank',      # Standard normalization for enrichment scores
    permutation_num=0,              # 0 is standard and fast for ssGSEA
    no_plot=True                    # We generate our own custom heatmap below
)

# Crucial Fix: GSEApy 1.3.0 returns a unpivoted long table in .res2d
# We pivot it into a true wide matrix: Rows = Cell Types (Term), Columns = Samples (Name)
print("Reshaping flat GSEApy output into wide matrix format...")
raw_results = ssgsea_run.res2d
es_matrix = raw_results.pivot(index='Term', columns='Name', values='ES')

# Save the pivoted enrichment matrix
scores_out_path = os.path.join(RESULTS_DIR, "authentic_immune_enrichment_scores.csv")
es_matrix.to_csv(scores_out_path)
print(f"✓ Matrix pivoted. Enrichment scores saved to:\n  {scores_out_path}")

# ============================================================================
# 5. VISUALIZE IMMUNE MICROENVIRONMENT CHANGES
# ============================================================================
print("\nGenerating publication-quality visualization...")

# Sort samples by clinical group to group them explicitly on the heatmap
df_sample_sorted = df_sample.sort_values(by="group")
es_matrix_sorted = es_matrix[df_sample_sorted.index]

# Filter down to the top 20 most highly variable signatures to capture major shifts
top_signatures = es_matrix_sorted.var(axis=1).nlargest(20).index
plot_data = es_matrix_sorted.loc[top_signatures]

# Robust cleaning for signature names (stripping complex technical prefixes safely)
cleaned_labels = []
for idx in plot_data.index:
    term_name = idx.split("___")[-1] if "___" in idx else idx
    cleaned_labels.append(term_name.replace("_", " ").title())
plot_data.index = cleaned_labels

# Map clinical groups to distinct, publication-grade color bars
group_colors = df_sample_sorted['group'].map({
    'normal': '#2ecc71',   # Emerald Green
    'eutopic': '#3498db',  # Marine Blue
    'ectopic': '#e74c3c',  # Crimson Red
    'diseased': '#f1c40f'  # Muted Gold
})

# Convert mapping to a clean matching list aligned to our matrix columns
col_colors = group_colors.reindex(plot_data.columns).tolist()

# Generate the Clustergram Heatmap
plt.figure(figsize=(14, 11))
g = sns.clustermap(
    plot_data,
    cmap="RdYlBu_r",       # Classic red-to-blue expression spectrum
    center=0,
    col_cluster=False,     # Maintain our strict clinical group order
    row_cluster=True,      # Cluster signatures that move together
    col_colors=col_colors,
    linewidths=0.3,
    figsize=(14, 11),
    cbar_kws={'label': 'Enrichment Score (ssGSEA)'}
)

# Adjust label alignments and font configurations
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=9)
plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)

# Superimpose the main chart title
g.fig.suptitle("Authentic Immune Landscape of Endometriosis Microenvironment\n(MSigDB C8 Signatures)", 
               y=1.02, fontsize=14, fontweight='bold')

# Save out to high resolution PDF
heatmap_path = os.path.join(FIGURES_DIR, "authentic_immune_landscape_heatmap.pdf")
g.savefig(heatmap_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Heatmap visualization successfully saved to:\n  {heatmap_path}")
print("\n*** Phase 3a execution successfully concluded! ***")