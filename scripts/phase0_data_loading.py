#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Project: Endometriosis Transcriptomics
Datasets: GSE25628 and GSE7305
Python script for downloading, processing, and QC with color coding.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS
import GEOparse

# ============================================================================
# Helper functions
# ============================================================================

def create_project_dirs(project_root):
    """Create standard directories under project_root."""
    data_raw = os.path.join(project_root, "data", "raw")
    data_processed = os.path.join(project_root, "data", "processed")
    figures = os.path.join(project_root, "figures")
    for d in [data_raw, data_processed, figures]:
        os.makedirs(d, exist_ok=True)
    return data_raw, data_processed, figures

def get_group_column(gse):
    """Auto-detect which column in phenotype data contains group information."""
    pheno = gse.phenotype_data
    candidates = ['source_name_ch1', 'characteristics_ch1', 'title']
    for col in candidates:
        if col in pheno.columns:
            return col
    raise ValueError("Could not find a column with group labels.")

def assign_groups(series, patterns):
    """
    Assign groups based on regex patterns.
    patterns: dict {group_name: list_of_compiled_regex_patterns}
    """
    groups = []
    for val in series:
        val_str = str(val) if pd.notna(val) else ""
        assigned = False
        for group, pat_list in patterns.items():
            for pat in pat_list:
                if pat.search(val_str):
                    groups.append(group)
                    assigned = True
                    break
            if assigned:
                break
        if not assigned:
            groups.append("unknown")
    return groups

def safe_log2_transform(expr, offset=0.5):
    """Apply log2(x + offset) only if max > 50 (suggesting raw scale)."""
    if expr.max().max() > 50:
        print("Detected raw expression values; applying log2(x + {})".format(offset))
        return np.log2(expr + offset)
    else:
        print("Expression values appear already log-transformed; skipping.")
        return expr

def plot_boxplot_individual(expr_df, sample_colors, group_colors, title, outfile):
    """Individual sample boxplot."""
    # Ensure sample_colors is a list in the same order as expr_df.columns
    color_list = [sample_colors[col] for col in expr_df.columns]
    fig, ax = plt.subplots(figsize=(14, 6))
    positions = np.arange(expr_df.shape[1])
    bp = ax.boxplot(expr_df.values, positions=positions, patch_artist=True,
                    showfliers=False, widths=0.6)
    for patch, color in zip(bp['boxes'], color_list):
        patch.set_facecolor(color)
    ax.set_xticks(positions)
    ax.set_xticklabels(expr_df.columns, rotation=90, fontsize=6)
    ax.set_title(title)
    ax.set_ylabel("Expression (log2 scale)")
    ax.set_xlabel("Samples")
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=col, label=lab) for lab, col in group_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_boxplot_grouped(expr_df, sample_info, group_colors, title, outfile):
    """Grouped boxplot (one box per group)."""
    groups = sorted(sample_info['group'].unique())
    # Remove 'unknown' if present
    groups = [g for g in groups if g != 'unknown']
    if not groups:
        print("No valid groups to plot.")
        return
    data_to_plot = []
    colors = []
    for g in groups:
        cols = sample_info[sample_info['group'] == g].index
        if len(cols) == 0:
            continue
        data_to_plot.append(expr_df[cols].values.flatten())
        colors.append(group_colors.get(g, 'gray'))
    if not data_to_plot:
        print("No data to plot for grouped boxplot.")
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    bp = ax.boxplot(data_to_plot, patch_artist=True, notch=False,
                    showfliers=True, widths=0.6)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    ax.set_xticklabels(groups)
    ax.set_title(title)
    ax.set_ylabel("Expression (log2 scale)")
    ax.set_xlabel("Group")
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_density(expr_df, sample_colors, group_colors, title, outfile):
    """Density plots for each sample, colored by group."""
    fig, ax = plt.subplots(figsize=(10, 6))
    for col in expr_df.columns:
        data = expr_df[col].dropna()
        if len(data) > 0:
            data.plot.kde(ax=ax, color=sample_colors[col], alpha=0.4, lw=0.8, legend=False)
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=col, label=lab) for lab, col in group_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_title(title)
    ax.set_xlabel("log2 Expression")
    ax.set_ylabel("Density")
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_mds(expr_df, sample_colors, group_colors, title, outfile):
    """MDS plot using sklearn's MDS."""
    # Transpose so samples are rows
    data = expr_df.T
    dist = pdist(data, metric='euclidean')
    dist_sq = squareform(dist)
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    coords = mds.fit_transform(dist_sq)
    # Build color list in the same order as data.index (sample IDs)
    color_list = [sample_colors[sample] for sample in data.index]
    fig, ax = plt.subplots(figsize=(10, 8))
    for i, (x, y) in enumerate(coords):
        ax.scatter(x, y, color=color_list[i], s=80, alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel("MDS Dimension 1")
    ax.set_ylabel("MDS Dimension 2")
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=col, label=lab) for lab, col in group_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_mean_barplot(expr_df, sample_colors, group_colors, title, outfile):
    """Barplot of mean expression per sample."""
    means = expr_df.mean(axis=0)
    # Build color list in the same order as means.index
    color_list = [sample_colors[sample] for sample in means.index]
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(means)), means.values, color=color_list, edgecolor='black')
    ax.set_xticks(range(len(means)))
    ax.set_xticklabels(means.index, rotation=90, fontsize=6)
    ax.set_title(title)
    ax.set_ylabel("Mean log2 Expression")
    ax.set_xlabel("Samples")
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=col, label=lab) for lab, col in group_colors.items()]
    ax.legend(handles=legend_elements, loc='upper right')
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

# ============================================================================
# Main processing function for one GEO series
# ============================================================================

def process_geo(geo_id, project_root, patterns, group_colors, log_transform=False, offset=0.5):
    """
    Download, process, and generate QC plots for a GEO series.
    """
    data_raw, data_processed, figures = create_project_dirs(project_root)
    
    print("=" * 70)
    print(f"Processing {geo_id}")
    print("=" * 70)
    
    # 1. Download GEO data (include full data tables)
    print(f"Downloading {geo_id} from GEO...")
    try:
        gse = GEOparse.get_GEO(geo=geo_id, destdir=data_raw, include_data=True)
    except Exception as e:
        print(f"Error downloading {geo_id}: {e}")
        return
    
    # 2. Extract expression matrix from each sample (GSM)
    print("Extracting expression data from samples...")
    sample_ids = list(gse.gsms.keys())
    if not sample_ids:
        print("ERROR: No samples found in GSE.")
        return
    
    # Get the first sample to identify column names
    first_sample = gse.gsms[sample_ids[0]]
    sample_table = first_sample.table
    # Identify the probe ID column and value column
    id_col = None
    value_col = None
    for col in sample_table.columns:
        if 'id' in col.lower() or 'probe' in col.lower():
            id_col = col
        if 'value' in col.lower():
            value_col = col
    if id_col is None or value_col is None:
        # Fallback: assume first column is ID, second is value
        id_col = sample_table.columns[0]
        value_col = sample_table.columns[1]
    print(f"Using ID column: '{id_col}', value column: '{value_col}'")
    
    # Build expression DataFrame: probes as index, samples as columns
    expr_dict = {}
    for gsm_id in sample_ids:
        sample = gse.gsms[gsm_id]
        if id_col not in sample.table.columns or value_col not in sample.table.columns:
            print(f"Warning: sample {gsm_id} missing expected columns. Skipping.")
            continue
        ser = sample.table.set_index(id_col)[value_col]
        expr_dict[gsm_id] = ser
    
    expr_df = pd.DataFrame(expr_dict)
    print(f"Expression matrix: {expr_df.shape[0]} probes x {expr_df.shape[1]} samples")
    
    # 3. Metadata
    pheno = gse.phenotype_data
    if pheno.empty:
        print("WARNING: No phenotype data found. Using sample names as grouping fallback.")
        pheno = pd.DataFrame({'title': sample_ids}, index=sample_ids)
    
    # 4. Detect grouping column
    group_col = get_group_column(gse)
    print(f"Using column '{group_col}' for grouping.")
    group_values = pheno[group_col].astype(str)
    print("Unique values in grouping column:")
    print(group_values.unique())
    
    # 5. Assign groups using the provided patterns
    groups = assign_groups(group_values, patterns)
    sample_info = pd.DataFrame({
        'sample_id': expr_df.columns,
        'source': group_values,
        'group': groups
    }, index=expr_df.columns)
    
    # Fallback: if all unknown, try using 'title'
    if (sample_info['group'] == 'unknown').all():
        print("WARNING: All samples assigned to 'unknown' from main column.")
        print("Trying to extract groups from 'title' column...")
        if 'title' in pheno.columns:
            title_vals = pheno['title'].astype(str)
            groups_fallback = assign_groups(title_vals, patterns)
            sample_info['group'] = groups_fallback
            print("Fallback group assignment:")
            print(sample_info['group'].value_counts())
    
    # Check again
    if (sample_info['group'] == 'unknown').all():
        print("ERROR: Still all unknown. Please adjust patterns based on the printed unique values.")
        return
    
    print("\nFinal group assignment:")
    print(sample_info['group'].value_counts())
    
    # 6. Apply log transformation if requested (but check if already log)
    if log_transform:
        expr_processed = safe_log2_transform(expr_df, offset)
    else:
        expr_processed = expr_df
    
    # 7. Save processed data
    expr_out = os.path.join(data_processed, f"{geo_id}_expression_matrix.csv")
    sample_out = os.path.join(data_processed, f"{geo_id}_sample_info.csv")
    expr_processed.to_csv(expr_out)
    sample_info.to_csv(sample_out)
    print(f"Saved expression matrix to {expr_out}")
    print(f"Saved sample info to {sample_out}")
    
    # 8. Create color vector (Series indexed by sample ID)
    sample_colors = sample_info['group'].map(group_colors)
    sample_colors = sample_colors.fillna('gray')
    
    # 9. Generate QC plots
    print("\nGenerating QC plots...")
    base_title = f"{geo_id}: Endometriosis Dataset"
    
    plot_boxplot_individual(expr_processed, sample_colors, group_colors,
                            f"{base_title} - Expression by Sample",
                            os.path.join(figures, f"{geo_id}_boxplot_individual.pdf"))
    plot_boxplot_grouped(expr_processed, sample_info, group_colors,
                         f"{base_title} - Expression by Group",
                         os.path.join(figures, f"{geo_id}_boxplot.pdf"))
    plot_density(expr_processed, sample_colors, group_colors,
                 f"{base_title} - Density",
                 os.path.join(figures, f"{geo_id}_density.pdf"))
    plot_mds(expr_processed, sample_colors, group_colors,
             f"{base_title} - MDS",
             os.path.join(figures, f"{geo_id}_MDS.pdf"))
    plot_mean_barplot(expr_processed, sample_colors, group_colors,
                      f"{base_title} - Mean Expression",
                      os.path.join(figures, f"{geo_id}_mean_expression.pdf"))
    
    print(f"\n✅ {geo_id} processing complete. All plots saved in {figures}.\n")

# ============================================================================
# Define patterns and colors for each dataset
# ============================================================================

# For GSE25628: groups = normal, ectopic, eutopic
patterns_25628 = {
    'normal': [re.compile(r'normal', re.I)],
    'ectopic': [re.compile(r'ectopic', re.I)],
    'eutopic': [re.compile(r'eutopic', re.I)]
}
group_colors_25628 = {
    'normal': 'blue',
    'ectopic': 'red',
    'eutopic': 'orange'
}

# For GSE7305: groups = normal, diseased
patterns_7305 = {
    'normal': [re.compile(r'normal', re.I)],
    'diseased': [re.compile(r'disease|endometrium/ovary-disease', re.I)]
}
group_colors_7305 = {
    'normal': 'blue',
    'diseased': 'red'
}

# ============================================================================
# Run for both datasets (adjust project_root as needed)
# ============================================================================

if __name__ == "__main__":
    # Set your project root directory. By default, use current working directory.
    project_root = os.getcwd()
    print(f"Project root: {project_root}")
    
    # Process GSE25628 (already log2, so log_transform=False)
    process_geo("GSE25628", project_root, patterns_25628, group_colors_25628,
                log_transform=False)
    
    # Process GSE7305 (raw values, so log_transform=True)
    process_geo("GSE7305", project_root, patterns_7305, group_colors_7305,
                log_transform=True, offset=0.5)
    
    print("All datasets processed successfully.")