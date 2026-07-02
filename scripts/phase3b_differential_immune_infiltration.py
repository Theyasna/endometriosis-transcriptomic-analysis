#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 3b: DIFFERENTIAL IMMUNE INFILTRATION STATISTICAL TESTING
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ============================================================================
# 1. CONFIGURATION - DYNAMIC PATH DETECTION
# ============================================================================
parser = argparse.ArgumentParser(description="Phase 3b: Differential Immune Stats")
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

ENRICHMENT_SCORES_PATH = os.path.join(PROJECT_DIR, "results", "immune_profiling", 
                                      "authentic_immune_enrichment_scores.csv")
SAMPLE_INFO_PATH = os.path.join(PROJECT_DIR, "data", "processed", "harmonized_sample_info.csv")

RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "immune_profiling")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "immune_profiling")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================================
# 2. LOAD & PREPARE DATA
# ============================================================================
print("Loading enrichment scores and sample metadata...")
es_matrix = pd.read_csv(ENRICHMENT_SCORES_PATH, index_col=0)
df_sample = pd.read_csv(SAMPLE_INFO_PATH, index_col=0)

# Ensure sample alignment across files
common_samples = df_sample.index.intersection(es_matrix.columns)
es_matrix = es_matrix[common_samples]
df_sample = df_sample.loc[common_samples]

groups_list = ['normal', 'eutopic', 'ectopic', 'diseased']

# ============================================================================
# 3. HIGH-THROUGHPUT STATISTICAL TESTING (830 SIGNATURES)
# ============================================================================
print("\nRunning Kruskal-Wallis Omnibus Tests & Pairwise Post-Hocs...")
stats_results = []

for cell_type in es_matrix.index:
    # Separate scores by clinical cohort
    group_data = {}
    is_valid = True
    for g in groups_list:
        samples_in_g = df_sample[df_sample['group'] == g].index
        scores = es_matrix.loc[cell_type, samples_in_g].values
        if len(scores) < 2:  # Ensure statistical group sizing is valid
            is_valid = False
        group_data[g] = scores
        
    if not is_valid:
        continue
        
    # Calculate Omnibus Significance across all 4 environments
    try:
        kw_stat, kw_p = stats.kruskal(
            group_data['normal'], 
            group_data['eutopic'], 
            group_data['ectopic'], 
            group_data['diseased']
        )
    except ValueError:
        kw_stat, kw_p = 0.0, 1.0

    # Helper function for pairwise post-hoc tests (Mann-Whitney U)
    def run_pairwise(g1, g2):
        try:
            return stats.mannwhitneyu(group_data[g1], group_data[g2], alternative='two-sided')[1]
        except ValueError:
            return 1.0

    # Compile data row
    row = {
        'Signature_Name': cell_type,
        'Clean_Name': cell_type.split("___")[-1].replace("_", " ").title() if "___" in cell_type else cell_type.replace("_", " ").title(),
        'Kruskal_Wallis_Stat': kw_stat,
        'Omnibus_p_value': kw_p,
        'Mean_Normal': np.mean(group_data['normal']),
        'Mean_Eutopic': np.mean(group_data['eutopic']),
        'Mean_Ectopic': np.mean(group_data['ectopic']),
        'Mean_Diseased': np.mean(group_data['diseased']),
        'p_Ectopic_vs_Normal': run_pairwise('ectopic', 'normal'),
        'p_Eutopic_vs_Normal': run_pairwise('eutopic', 'normal'),
        'p_Ectopic_vs_Eutopic': run_pairwise('ectopic', 'eutopic'),
        'p_Diseased_vs_Normal': run_pairwise('diseased', 'normal')
    }
    stats_results.append(row)

df_stats = pd.DataFrame(stats_results)

# Apply False Discovery Rate (FDR) adjustment to control for multi-testing overhead
df_stats['Omnibus_FDR'] = multipletests(df_stats['Omnibus_p_value'], method='fdr_bh')[1]
df_stats = df_stats.sort_values(by='Omnibus_p_value')

# Save comprehensive statistical analysis sheet
stats_out_path = os.path.join(RESULTS_DIR, "differential_immune_infiltration_stats.csv")
df_stats.to_csv(stats_out_path, index=False)
print(f"✓ Statistical profiling complete. Full table saved to:\n  {stats_out_path}")

# Filter out significantly shifting signatures
sig_signatures = df_stats[df_stats['Omnibus_FDR'] < 0.05]
print(f"✓ Found {len(sig_signatures)} signatures showing significant cohort variance at FDR < 0.05.")

# ============================================================================
# 4. GRID VISUALIZATION OF TOP TARGET IMMUNE POPULATIONS
# ============================================================================
print("\nPlotting top 6 highly significant immune signature transitions...")

# Select the top 6 moving signatures by raw p-value
top_6_targets = df_stats.head(6)

fig, axes = plt.subplots(2, 3, figsize=(18, 11), sharex=False)
axes = axes.flatten()

# Custom design color mapping matching your heatmap environment bar
palette_colors = {
    'normal': '#2ecc71',   # Emerald Green
    'eutopic': '#3498db',  # Marine Blue
    'ectopic': '#e74c3c',  # Crimson Red
    'diseased': '#f1c40f'  # Muted Gold
}

for i, (_, row) in enumerate(top_6_targets.iterrows()):
    ax = axes[i]
    sig_id = row['Signature_Name']
    clean_label = row['Clean_Name']
    
    # Isolate vector scores across cohorts for visualization
    plot_rows = []
    for g in groups_list:
        samples_in_g = df_sample[df_sample['group'] == g].index
        for s in samples_in_g:
            plot_rows.append({'Group': g, 'Score': es_matrix.loc[sig_id, s]})
    df_plot_sig = pd.DataFrame(plot_rows)
    
    # Generate hybrid Box+Strip charts for high publication clarity
    sns.boxplot(
        data=df_plot_sig, x='Group', y='Score', hue='Group', order=groups_list,
        palette=palette_colors, ax=ax, width=0.5, fliersize=0, legend=False
    )
    sns.stripplot(
        data=df_plot_sig, x='Group', y='Score', order=groups_list,
        color='black', alpha=0.4, size=4, jitter=0.2, ax=ax
    )
    
    # Refine specific text features
    ax.set_title(f"{clean_label}\nOmnibus FDR: {row['Omnibus_FDR']:.2e}", fontsize=11, fontweight='bold')
    ax.set_ylabel("Enrichment Score (ssGSEA)", fontsize=9)
    ax.set_xlabel("", fontsize=9)
    ax.tick_params(axis='x', labelsize=10)
    ax.grid(axis='y', linestyle='--', alpha=0.5)

# Adjust overall layout spacing
plt.tight_layout()
fig.subplots_adjust(top=0.90)
fig.suptitle("Top Significant Immune Profiling Modifications Across Endometriosis Microenvironments", 
             fontsize=15, fontweight='bold')

boxplots_out_path = os.path.join(FIGURES_DIR, "differential_immune_boxplots.pdf")
plt.savefig(boxplots_out_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Publication-grade grid boxplots saved to:\n  {boxplots_out_path}")
print("\n*** Phase 3b statistical validation concluded cleanly! ***")