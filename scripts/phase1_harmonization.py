#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 1: HARMONIZATION & BATCH CORRECTION (Fully Fixed)
"""

import os
import argparse
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Try to import pycombat
try:
    from pycombat import Combat
    print("✓ pycombat imported successfully")
except ImportError:
    print("ERROR: pycombat not installed. Install with: pip install pycombat")
    exit(1)

# ============================================================================
# CONFIGURATION
# ============================================================================
parser = argparse.ArgumentParser(description="Phase 1: Harmonization")
parser.add_argument("--project_root", default=os.getcwd(), 
                    help="Project root directory")
args = parser.parse_args()

PROJECT_DIR = args.project_root
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "harmonization")

os.makedirs(FIGURES_DIR, exist_ok=True)

print("\n" + "="*70)
print("PHASE 1: HARMONIZATION & BATCH CORRECTION")
print("="*70)
print(f"Project directory: {PROJECT_DIR}")

# ============================================================================
# STEP 1-4: LOAD, ALIGN, MERGE
# ============================================================================
print("\nSTEP 1: Loading expression matrices...")

# Check if files exist
expr_25628_path = os.path.join(DATA_PROCESSED, "GSE25628_expression_matrix.csv")
expr_7305_path = os.path.join(DATA_PROCESSED, "GSE7305_expression_matrix.csv")
sample_25628_path = os.path.join(DATA_PROCESSED, "GSE25628_sample_info.csv")
sample_7305_path = os.path.join(DATA_PROCESSED, "GSE7305_sample_info.csv")

if not all(os.path.exists(f) for f in [expr_25628_path, expr_7305_path, 
                                        sample_25628_path, sample_7305_path]):
    print("ERROR: Missing files from Phase 0. Please run phase0 first.")
    exit(1)

expr_25628 = pd.read_csv(expr_25628_path, index_col=0)
sample_25628 = pd.read_csv(sample_25628_path, index_col=0)

expr_7305 = pd.read_csv(expr_7305_path, index_col=0)
sample_7305 = pd.read_csv(sample_7305_path, index_col=0)

print(f"✓ GSE25628: {expr_25628.shape[0]} probes × {expr_25628.shape[1]} samples")
print(f"✓ GSE7305:  {expr_7305.shape[0]} probes × {expr_7305.shape[1]} samples")

# Align probes
common_probes = expr_25628.index.intersection(expr_7305.index)
expr_25628_aligned = expr_25628.loc[common_probes].copy()
expr_7305_aligned = expr_7305.loc[common_probes].copy()

print(f"✓ Aligned to {len(common_probes)} common probes")

# Batch labels & merge
sample_25628['batch'] = 'GSE25628'
sample_7305['batch'] = 'GSE7305'
sample_combined = pd.concat([sample_25628, sample_7305])
expr_merged = pd.concat([expr_25628_aligned, expr_7305_aligned], axis=1)

print(f"✓ Merged: {expr_merged.shape[0]} probes × {expr_merged.shape[1]} samples")
assert all(expr_merged.columns == sample_combined.index), "Sample order mismatch!"

# ============================================================================
# STEP 5: PCA BEFORE
# ============================================================================
print("\nSTEP 5: Visualizing BEFORE correction...")

def plot_pca(expr_data, sample_info, title, filename, color_by='batch'):
    """Plot PCA colored by batch or group."""
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(expr_data.T)
    
    pca_df = pd.DataFrame(pca_coords, columns=['PC1', 'PC2'], index=expr_data.columns)
    pca_df['batch'] = sample_info['batch'].values
    pca_df['group'] = sample_info['group'].values
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot by batch
    for batch in pca_df['batch'].unique():
        subset = pca_df[pca_df['batch'] == batch]
        axes[0].scatter(subset['PC1'], subset['PC2'], label=batch, alpha=0.7, s=80)
    axes[0].set_title(f"{title} - Colored by Batch")
    axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot by group
    for group in pca_df['group'].unique():
        subset = pca_df[pca_df['group'] == group]
        axes[1].scatter(subset['PC1'], subset['PC2'], label=group, alpha=0.7, s=80)
    axes[1].set_title(f"{title} - Colored by Group")
    axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    axes[1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    return pca_df, pca.explained_variance_ratio_

# Before correction
pca_before, var_before = plot_pca(
    expr_merged, sample_combined, 
    "PCA BEFORE Batch Correction",
    os.path.join(FIGURES_DIR, "01_PCA_before_correction.pdf")
)
print(f"✓ PCA before correction saved (PC1: {var_before[0]:.1%}, PC2: {var_before[1]:.1%})")

# ============================================================================
# STEP 6: COMBAT
# ============================================================================
print("\nSTEP 6: Applying ComBat...")

# Prepare batch vector (0 for GSE25628, 1 for GSE7305)
batch_vector = np.where(sample_combined['batch'] == 'GSE25628', 0, 1)

print(f"  Batch distribution: GSE25628={sum(batch_vector==0)}, GSE7305={sum(batch_vector==1)}")

# Apply ComBat
try:
    combat_model = Combat()
    expr_corrected_array = combat_model.fit_transform(
        expr_merged.values.T, batch_vector
    ).T
    
    expr_corrected = pd.DataFrame(
        expr_corrected_array,
        index=expr_merged.index,
        columns=expr_merged.columns
    )
    print("✓ ComBat completed successfully")
except Exception as e:
    print(f"ERROR in ComBat: {e}")
    print("Trying alternative method...")
    # Alternative: Use combat with different parameters
    combat_model = Combat()
    expr_corrected_array = combat_model.fit_transform(
        expr_merged.values.T, 
        batch_vector,
        mean_only=False
    ).T
    expr_corrected = pd.DataFrame(
        expr_corrected_array,
        index=expr_merged.index,
        columns=expr_merged.columns
    )
    print("✓ ComBat completed with alternative parameters")

# ============================================================================
# STEP 7: PCA AFTER
# ============================================================================
print("\nSTEP 7: Visualizing AFTER correction...")

pca_after, var_after = plot_pca(
    expr_corrected, sample_combined,
    "PCA AFTER Batch Correction",
    os.path.join(FIGURES_DIR, "02_PCA_after_correction.pdf")
)
print(f"✓ PCA after correction saved (PC1: {var_after[0]:.1%}, PC2: {var_after[1]:.1%})")

# ============================================================================
# STEP 8: COMPARISON PLOT
# ============================================================================
print("\nSTEP 8: Creating comparison plots...")

# Create combined before/after plot
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Before - by batch
pca_before_df = pca_before
for batch in pca_before_df['batch'].unique():
    subset = pca_before_df[pca_before_df['batch'] == batch]
    axes[0, 0].scatter(subset['PC1'], subset['PC2'], label=batch, alpha=0.7, s=80)
axes[0, 0].set_title("BEFORE: By Batch")
axes[0, 0].set_xlabel("PC1")
axes[0, 0].set_ylabel("PC2")
axes[0, 0].legend()
axes[0, 0].grid(True, alpha=0.3)

# Before - by group
for group in pca_before_df['group'].unique():
    subset = pca_before_df[pca_before_df['group'] == group]
    axes[0, 1].scatter(subset['PC1'], subset['PC2'], label=group, alpha=0.7, s=80)
axes[0, 1].set_title("BEFORE: By Group")
axes[0, 1].set_xlabel("PC1")
axes[0, 1].set_ylabel("PC2")
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# After - by batch
pca_after_df = pca_after
for batch in pca_after_df['batch'].unique():
    subset = pca_after_df[pca_after_df['batch'] == batch]
    axes[1, 0].scatter(subset['PC1'], subset['PC2'], label=batch, alpha=0.7, s=80)
axes[1, 0].set_title("AFTER: By Batch")
axes[1, 0].set_xlabel("PC1")
axes[1, 0].set_ylabel("PC2")
axes[1, 0].legend()
axes[1, 0].grid(True, alpha=0.3)

# After - by group
for group in pca_after_df['group'].unique():
    subset = pca_after_df[pca_after_df['group'] == group]
    axes[1, 1].scatter(subset['PC1'], subset['PC2'], label=group, alpha=0.7, s=80)
axes[1, 1].set_title("AFTER: By Group")
axes[1, 1].set_xlabel("PC1")
axes[1, 1].set_ylabel("PC2")
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.suptitle("Batch Correction Comparison", fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "03_comparison_before_after.pdf"), dpi=300, bbox_inches='tight')
plt.close()
print("✓ Comparison plot saved")

# ============================================================================
# STEP 9: SAVE
# ============================================================================
print("\nSTEP 9: Saving harmonized data...")

# Save corrected expression
expr_corrected.to_csv(os.path.join(DATA_PROCESSED, "harmonized_expression_matrix.csv"))
print(f"✓ Saved: {os.path.join(DATA_PROCESSED, 'harmonized_expression_matrix.csv')}")

# Save sample info
sample_combined.to_csv(os.path.join(DATA_PROCESSED, "harmonized_sample_info.csv"))
print(f"✓ Saved: {os.path.join(DATA_PROCESSED, 'harmonized_sample_info.csv')}")

# ============================================================================
# STEP 10: QC SUMMARY
# ============================================================================
print("\nSTEP 10: Generating QC summary...")

summary_file = os.path.join(FIGURES_DIR, "harmonization_summary.txt")
with open(summary_file, 'w') as f:
    f.write("="*70 + "\n")
    f.write("HARMONIZATION QC SUMMARY\n")
    f.write("="*70 + "\n\n")
    
    f.write("DATASET OVERVIEW:\n")
    f.write(f"  - GSE25628: {expr_25628.shape[1]} samples, {expr_25628.shape[0]} probes\n")
    f.write(f"  - GSE7305:  {expr_7305.shape[1]} samples, {expr_7305.shape[0]} probes\n")
    f.write(f"  - Common probes: {len(common_probes)}\n")
    f.write(f"  - Total samples after merge: {expr_merged.shape[1]}\n\n")
    
    f.write("BATCH CORRECTION:\n")
    f.write(f"  - Method: ComBat\n")
    f.write(f"  - Before correction - PC1: {var_before[0]:.2%}, PC2: {var_before[1]:.2%}\n")
    f.write(f"  - After correction - PC1: {var_after[0]:.2%}, PC2: {var_after[1]:.2%}\n")
    f.write(f"  - PC1 change: {(var_after[0] - var_before[0]):.2%}\n\n")
    
    f.write("GROUP DISTRIBUTION:\n")
    for group, count in sample_combined['group'].value_counts().items():
        f.write(f"  - {group}: {count} samples\n")
    f.write("\n")
    
    f.write("BATCH DISTRIBUTION:\n")
    for batch, count in sample_combined['batch'].value_counts().items():
        f.write(f"  - {batch}: {count} samples\n")

print(f"✓ QC summary saved to: {summary_file}")

# ============================================================================
# COMPLETE
# ============================================================================
print("\n" + "="*70)
print("✅ PHASE 1 COMPLETE")
print("="*70)
print(f"Output files saved to: {DATA_PROCESSED}")
print(f"Figures saved to: {FIGURES_DIR}")
print("\n✓✓✓ Ready for Phase 2: Differential Expression Analysis ✓✓✓")