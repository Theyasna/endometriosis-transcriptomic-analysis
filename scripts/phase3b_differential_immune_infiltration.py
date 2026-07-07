#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 3b: DIFFERENTIAL IMMUNE INFILTRATION (PER-DATASET)
FIXED: No R syntax in Python, dynamic grid, deterministic group order
"""

import os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ============================================================================
# 1. CONFIGURATION
# ============================================================================
parser = argparse.ArgumentParser(description="Phase 3b: Differential Immune Stats")
parser.add_argument("--project_root", default=os.getcwd())
args = parser.parse_args()

PROJECT_DIR = args.project_root
if not os.path.exists(os.path.join(PROJECT_DIR, "data")):
    parent_dir = os.path.dirname(PROJECT_DIR)
    if os.path.exists(os.path.join(parent_dir, "data")):
        PROJECT_DIR = parent_dir

print(f"Project directory: {PROJECT_DIR}")

ENRICHMENT_SCORES_PATH = os.path.join(
    PROJECT_DIR, "results", "immune_profiling", "immunesigdb_enrichment_scores.csv"
)
SAMPLE_INFO_PATH = os.path.join(PROJECT_DIR, "data", "processed", "harmonized_sample_info.csv")

RESULTS_DIR = os.path.join(PROJECT_DIR, "results", "immune_profiling")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "immune_profiling")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================================
# 2. LOAD DATA
# ============================================================================
print("Loading enrichment scores and sample metadata...")
es_matrix = pd.read_csv(ENRICHMENT_SCORES_PATH, index_col=0)
df_sample = pd.read_csv(SAMPLE_INFO_PATH, index_col=0)

common_samples = df_sample.index.intersection(es_matrix.columns)
es_matrix = es_matrix[common_samples]
df_sample = df_sample.loc[common_samples]

# ============================================================================
# 3. PER-DATASET STATISTICAL TESTING
# ============================================================================
print("\nRunning per-dataset Kruskal-Wallis tests...")

datasets = {
    "GSE25628": {
        "samples": df_sample[df_sample['batch'] == "GSE25628"].index,
        "groups": ["ectopic", "eutopic", "normal"]
    },
    "GSE7305": {
        "samples": df_sample[df_sample['batch'] == "GSE7305"].index,
        "groups": ["diseased", "normal"]
    }
}

all_stats = []

for ds_name, ds_info in datasets.items():
    print(f"\n  Processing: {ds_name}")

    ds_samples = ds_info["samples"]
    ds_es = es_matrix[ds_samples]
    ds_sample_info = df_sample.loc[ds_samples]

    # Use explicit group list for deterministic order
    groups_available = [g for g in ds_info["groups"] if g in ds_sample_info['group'].unique()]
    print(f"    Groups: {groups_available}")

    if len(groups_available) < 2:
        print(f"    Skipping {ds_name} — fewer than 2 groups")
        continue

    for cell_type in ds_es.index:
        group_data = {}
        for g in groups_available:
            samples_in_g = ds_sample_info[ds_sample_info['group'] == g].index
            scores = ds_es.loc[cell_type, samples_in_g].values
            if len(scores) >= 2:
                group_data[g] = scores

        if len(group_data) < 2:
            continue

        try:
            # FIX: iterate over group_data (groups that actually passed the
            # >=2 sample filter) rather than groups_available. If a group had
            # <2 samples it was dropped from group_data, and the old
            # `for g in groups_available` would raise KeyError -- silently
            # swallowed by the except below and recorded as p=1.0. Harmless on
            # the current data (no group drops below 2) but a real trap if a
            # smaller dataset is ever used.
            kw_stat, kw_p = stats.kruskal(*[group_data[g] for g in group_data])
        except ValueError:
            kw_stat, kw_p = 0.0, 1.0

        def run_pairwise(g1, g2):
            try:
                return stats.mannwhitneyu(group_data[g1], group_data[g2],
                                         alternative='two-sided')[1]
            except Exception:
                return 1.0

        row = {
            'Dataset': ds_name,
            'Signature': cell_type,
            'Clean_Name': cell_type.split("___")[-1].replace("_", " ").title()
                          if "___" in cell_type else cell_type.replace("_", " ").title(),
            'Kruskal_Wallis_Stat': kw_stat,
            'Omnibus_p_value': kw_p,
        }

        # FIX: key off group_data (groups that passed the >=2 filter),
        # consistent with the kruskal call above, so column names and the
        # test never reference a group that was dropped for small n.
        present_groups = list(group_data.keys())
        for g in present_groups:
            row[f'Mean_{g.title()}'] = np.mean(group_data[g])

        for i, g1 in enumerate(present_groups):
            for g2 in present_groups[i+1:]:
                row[f'p_{g1}_vs_{g2}'] = run_pairwise(g1, g2)

        all_stats.append(row)

df_stats = pd.DataFrame(all_stats)

if df_stats.empty:
    print("\nERROR: No valid comparisons. Check your sample metadata.")
    exit(1)

# Apply FDR per dataset
for ds_name in datasets.keys():
    mask = df_stats['Dataset'] == ds_name
    if mask.any():
        pvals = df_stats.loc[mask, 'Omnibus_p_value'].values
        # Handle NaNs in p-values
        pvals = np.nan_to_num(pvals, nan=1.0)
        df_stats.loc[mask, 'Omnibus_FDR'] = multipletests(pvals, method='fdr_bh')[1]

# ============================================================================
# 4. SAVE RESULTS
# ============================================================================
stats_out_path = os.path.join(RESULTS_DIR, "differential_immune_infiltration_stats.csv")
df_stats.to_csv(stats_out_path, index=False)
print(f"\nStatistical results saved to: {stats_out_path}")

for ds_name in datasets.keys():
    mask = df_stats['Dataset'] == ds_name
    if mask.any():
        # Use pandas sum with skipna (Python syntax, not R)
        n_sig = df_stats.loc[mask, 'Omnibus_FDR'].lt(0.05).sum()
        n_total = mask.sum()
        print(f"  {ds_name}: {n_sig}/{n_total} significant (FDR < 0.05)")

# ============================================================================
# 5. BOXPLOTS (Top 6 per dataset)
# ============================================================================
print("\nPlotting top significant signatures...")

top_sigs = []
for ds_name in datasets.keys():
    mask = df_stats['Dataset'] == ds_name
    ds_sigs = df_stats[mask].sort_values('Omnibus_p_value').head(6)
    if not ds_sigs.empty:
        top_sigs.append(ds_sigs)

if top_sigs:
    top_combined = pd.concat(top_sigs)
    n_plots = len(top_combined)

    # Dynamic grid: up to 6 per dataset, max 12 total
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols
    # Cap at 2 rows for readability
    n_rows = min(n_rows, 2)
    n_cols = min(n_plots, 3)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
    if n_rows == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    palette_colors = {
        'normal': '#2ecc71',
        'eutopic': '#3498db',
        'ectopic': '#e74c3c',
        'diseased': '#f1c40f'
    }

    for i, (_, row) in enumerate(top_combined.iterrows()):
        if i >= len(axes):
            break

        ax = axes[i]
        sig_id = row['Signature']
        ds_name = row['Dataset']

        ds_samples = df_sample[df_sample['batch'] == ds_name].index
        plot_es = es_matrix.loc[sig_id, ds_samples]
        plot_df = pd.DataFrame({
            'Group': df_sample.loc[ds_samples, 'group'],
            'Score': plot_es.values
        })

        groups_order = [g for g in ['normal', 'ectopic', 'eutopic', 'diseased']
                       if g in plot_df['Group'].unique()]

        sns.boxplot(data=plot_df, x='Group', y='Score', hue='Group',
                   order=groups_order, palette=palette_colors, ax=ax,
                   width=0.5, fliersize=0, legend=False)
        sns.stripplot(data=plot_df, x='Group', y='Score', order=groups_order,
                     color='black', alpha=0.4, size=4, jitter=0.2, ax=ax)

        fdr_val = row['Omnibus_FDR']
        ax.set_title(f"{row['Clean_Name']}\n{ds_name} - FDR: {fdr_val:.2e}",
                    fontsize=10, fontweight='bold')
        ax.set_ylabel("Enrichment Score")
        ax.set_xlabel("")
        ax.tick_params(axis='x', labelsize=9)
        ax.grid(axis='y', linestyle='--', alpha=0.5)

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "differential_immune_boxplots.pdf"),
               dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Boxplots saved to: {FIGURES_DIR}/differential_immune_boxplots.pdf")

print("\nPhase 3b complete.")