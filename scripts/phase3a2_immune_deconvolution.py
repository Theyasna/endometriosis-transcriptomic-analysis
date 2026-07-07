#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 3a: IMMUNE SIGNATURE DECONVOLUTION (ssGSEA) via MSigDB C7 ImmuneSigDB
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy import Msigdb

# ============================================================================
# 1. CONFIGURATION - DYNAMIC PATH DETECTION
# ============================================================================
parser = argparse.ArgumentParser(description="Phase 3a: Immune Deconvolution")
parser.add_argument("--project_root", default=os.getcwd(), 
                    help="Project root directory")
args = parser.parse_args()

PROJECT_DIR = args.project_root

# Auto-detect if running from scripts folder
if not os.path.exists(os.path.join(PROJECT_DIR, "data")):
    parent_dir = os.path.dirname(PROJECT_DIR)
    if os.path.exists(os.path.join(parent_dir, "data")):
        PROJECT_DIR = parent_dir

print(f"Project directory: {PROJECT_DIR}")

GENE_EXPR_PATH = os.path.join(PROJECT_DIR, "data", "processed", "gene_expression_matrix.csv")
SAMPLE_INFO_PATH = os.path.join(PROJECT_DIR, "data", "processed", "harmonized_sample_info.csv")

RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "immune_profiling")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "immune_profiling")

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
print(f"Data Synced: {df_genes.shape[0]} unique genes across {df_genes.shape[1]} samples.")

# ============================================================================
# 3. FETCH IMMUNE-SPECIFIC GENE SETS (FIXED)
# ============================================================================
# FIX: the original script fetched MSigDB category 'c8' (a pan-tissue,
# cross-organ single-cell reference atlas -- fetal lung, olfactory
# epithelium, developing heart, midbrain, etc.) and labeled the result
# "authentic immune deconvolution." Cross-checking the 707 signatures
# actually returned, only ~18% matched immune-cell keywords; the
# remaining ~82% were unrelated cell types. That mismatch produced
# statistically significant hits (e.g. "Fetal Heart Smooth Muscle
# Cells") that have no immunological interpretation.
#
# C7 (ImmuneSigDB / immunologic signature gene sets) is MSigDB's actual
# immune-focused collection, built from published immune cell expression
# and perturbation studies. This is the correct category for a claim
# about "immune landscape" or "immune infiltration."
print("\nConnecting to MSigDB to fetch C7 immunologic signature gene sets...")
try:
    msig = Msigdb()
    immune_genesets = msig.get_gmt(category='c7.immunesigdb')
    print(f"Successfully loaded {len(immune_genesets)} immune-specific gene sets (C7 ImmuneSigDB).")
except Exception as e:
    print(f"WARNING: C7 ImmuneSigDB fetch failed: {e}")
    print("Attempting fallback: full C7 category...")
    msig = Msigdb()
    immune_genesets = msig.get_gmt(category='c7')
    print(f"Loaded {len(immune_genesets)} fallback immune gene sets (C7).")

# ============================================================================
# 4. RUN ssGSEA DECONVOLUTION
# ============================================================================
print("\nRunning single-sample GSEA (ssGSEA)...")
ssgsea_run = gp.ssgsea(
    data=df_genes,
    gene_sets=immune_genesets,      # C7 ImmuneSigDB gene sets
    outdir=None,                    # Keep results in memory for processing
    sample_norm_method='rank',      # Standard normalization for enrichment scores
    permutation_num=0,              # 0 is standard and fast for ssGSEA
    no_plot=True                    # We generate our own custom heatmap below
)

# Note: GSEApy 1.3.0 returns a unpivoted long table in .res2d
# We pivot it into a true wide matrix: Rows = Cell Types (Term), Columns = Samples (Name)
print("Reshaping flat GSEApy output into wide matrix format...")
raw_results = ssgsea_run.res2d
es_matrix = raw_results.pivot(index='Term', columns='Name', values='ES')

# Save the pivoted enrichment matrix
# FIX: renamed from "authentic_immune_..." (C8, mislabeled) to reflect
# that this is now genuinely immune-signature-based (C7 ImmuneSigDB)
scores_out_path = os.path.join(RESULTS_DIR, "immunesigdb_enrichment_scores.csv")
es_matrix.to_csv(scores_out_path)
print(f"Matrix pivoted. Enrichment scores saved to:\n  {scores_out_path}")

# ============================================================================
# 5. VISUALIZE IMMUNE MICROENVIRONMENT CHANGES
# ============================================================================
print("\nGenerating immune signature heatmap...")

# Sort samples by clinical group to group them explicitly on the heatmap
df_sample_sorted = df_sample.sort_values(by="group")
es_matrix_sorted = es_matrix[df_sample_sorted.index]

# Filter down to the top 20 most highly variable signatures to capture major shifts
top_signatures = es_matrix_sorted.var(axis=1).nlargest(20).index
plot_data = es_matrix_sorted.loc[top_signatures]

# Clean signature names for display (strip technical prefixes)
cleaned_labels = []
for idx in plot_data.index:
    term_name = idx.split("___")[-1] if "___" in idx else idx
    cleaned_labels.append(term_name.replace("_", " ").title())
plot_data.index = cleaned_labels

# Map clinical groups to color bars
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
g.fig.suptitle("Immune Signature Landscape of Endometriosis Microenvironment\n(MSigDB C7 ImmuneSigDB)",
               y=1.02, fontsize=14, fontweight='bold')

# Save out to high resolution PDF
heatmap_path = os.path.join(FIGURES_DIR, "immunesigdb_landscape_heatmap.pdf")
g.savefig(heatmap_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Heatmap visualization successfully saved to:\n  {heatmap_path}")
print("\nPhase 3a complete.")