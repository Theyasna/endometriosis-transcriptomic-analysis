#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 1: HARMONIZATION & BATCH CORRECTION

Workflow:
  1. Load GSE25628 + GSE7305 expression matrices
  2. Align probes (find common genes)
  3. Merge into single matrix
  4. Apply ComBat batch correction
  5. Visualize batch effect removal
  6. Save harmonized dataset

Expected runtime: 5-10 minutes
Output: Harmonized expression matrix (42 samples × ~22k genes)
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy import stats

# Try to import pycombat (ComBat batch correction)
try:
    from pycombat import Combat
    print("✓ pycombat imported successfully")
except ImportError:
    print("ERROR: pycombat not installed. Install with:")
    print("  pip install pycombat")
    exit(1)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Update this to your actual project path
PROJECT_DIR = "C:/Users/Yasna/OneDrive/Belgeler/endometriosis-transcriptomic-analysis"

# Define directories
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "harmonization")

# Create harmonization figures directory
os.makedirs(FIGURES_DIR, exist_ok=True)

print("\n" + "="*70)
print("PHASE 1: HARMONIZATION & BATCH CORRECTION")
print("="*70)
print(f"Project directory: {PROJECT_DIR}")
print(f"Data directory: {DATA_PROCESSED}")
print(f"Output figures: {FIGURES_DIR}\n")

# ============================================================================
# STEP 1: LOAD EXPRESSION MATRICES
# ============================================================================
print("STEP 1: Loading expression matrices...")

# Load GSE25628 (RMA-normalized, log2 scale)
expr_25628_file = os.path.join(DATA_PROCESSED, "GSE25628_expression_matrix.csv")
sample_25628_file = os.path.join(DATA_PROCESSED, "GSE25628_sample_info.csv")

if not os.path.exists(expr_25628_file):
    print(f"ERROR: GSE25628 file not found at {expr_25628_file}")
    print("Please ensure Phase 0 (data loading) completed successfully.")
    exit(1)

expr_25628 = pd.read_csv(expr_25628_file, index_col=0)
sample_25628 = pd.read_csv(sample_25628_file, index_col=0)

print(f"✓ GSE25628 loaded: {expr_25628.shape[0]} probes × {expr_25628.shape[1]} samples")
print(f"  Sample groups: {sample_25628['group'].value_counts().to_dict()}")

# Load GSE7305 (RMA-normalized, log2 scale — transformed in Phase 0)
expr_7305_file = os.path.join(DATA_PROCESSED, "GSE7305_expression_matrix.csv")
sample_7305_file = os.path.join(DATA_PROCESSED, "GSE7305_sample_info.csv")

if not os.path.exists(expr_7305_file):
    print(f"ERROR: GSE7305 file not found at {expr_7305_file}")
    exit(1)

expr_7305 = pd.read_csv(expr_7305_file, index_col=0)
sample_7305 = pd.read_csv(sample_7305_file, index_col=0)

print(f"✓ GSE7305 loaded: {expr_7305.shape[0]} probes × {expr_7305.shape[1]} samples")
print(f"  Sample groups: {sample_7305['group'].value_counts().to_dict()}")

# ============================================================================
# STEP 2: ALIGN PROBES (FIND COMMON GENES)
# ============================================================================
print("\nSTEP 2: Aligning probes...")

# Both datasets use Affymetrix probes, so row indices are probe IDs
# Find probes present in BOTH datasets
common_probes = expr_25628.index.intersection(expr_7305.index)

print(f"  GSE25628 probes: {expr_25628.shape[0]}")
print(f"  GSE7305 probes: {expr_7305.shape[0]}")
print(f"  Common probes: {len(common_probes)}")
print(f"  Retention: {len(common_probes) / expr_25628.shape[0] * 100:.1f}% of GSE25628")

# Filter both to common probes only
expr_25628_aligned = expr_25628.loc[common_probes].copy()
expr_7305_aligned = expr_7305.loc[common_probes].copy()

print(f"✓ Both datasets aligned to {len(common_probes)} common probes")

# ============================================================================
# STEP 3: CREATE BATCH LABELS
# ============================================================================
print("\nSTEP 3: Creating batch labels...")

# Add batch information to sample info
sample_25628['batch'] = 'GSE25628'
sample_7305['batch'] = 'GSE7305'

# Combine sample info
sample_combined = pd.concat([sample_25628, sample_7305], axis=0)

print(f"✓ Batch labels assigned:")
print(f"  GSE25628: {(sample_combined['batch'] == 'GSE25628').sum()} samples")
print(f"  GSE7305: {(sample_combined['batch'] == 'GSE7305').sum()} samples")

# ============================================================================
# STEP 4: MERGE EXPRESSION MATRICES
# ============================================================================
print("\nSTEP 4: Merging expression matrices...")

# Concatenate horizontally (samples as columns)
expr_merged = pd.concat([expr_25628_aligned, expr_7305_aligned], axis=1)

print(f"✓ Merged matrix: {expr_merged.shape[0]} probes × {expr_merged.shape[1]} samples")

# Verify sample order matches
assert all(expr_merged.columns == sample_combined.index), "Sample order mismatch!"
print("✓ Sample order verified")

# ============================================================================
# STEP 5: VISUALIZE BATCH EFFECT (BEFORE CORRECTION)
# ============================================================================
print("\nSTEP 5: Visualizing batch effects (BEFORE correction)...")

# Use PCA to show batch effect
print("  Computing PCA...")
pca = PCA(n_components=2, random_state=42)
pca_coords = pca.fit_transform(expr_merged.T)  # Transpose: samples × probes

pca_df = pd.DataFrame(
    pca_coords,
    columns=['PC1', 'PC2'],
    index=expr_merged.columns
)
pca_df['batch'] = sample_combined['batch'].values
pca_df['group'] = sample_combined['group'].values

# Color by batch
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Colored by batch
batch_colors = {'GSE25628': 'red', 'GSE7305': 'blue'}
for batch in batch_colors:
    mask = pca_df['batch'] == batch
    axes[0].scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], 
                   label=batch, alpha=0.7, s=100, 
                   color=batch_colors[batch])
axes[0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
axes[0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
axes[0].set_title("BEFORE Batch Correction (colored by dataset)")
axes[0].legend()
axes[0].grid(alpha=0.3)

# Plot 2: Colored by group
group_colors = {'normal': 'blue', 'ectopic': 'red', 'diseased': 'orange'}
for group in group_colors:
    mask = pca_df['group'] == group
    if mask.sum() > 0:
        axes[1].scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], 
                       label=group, alpha=0.7, s=100, 
                       color=group_colors[group])
axes[1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
axes[1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
axes[1].set_title("BEFORE Batch Correction (colored by group)")
axes[1].legend()
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "01_PCA_before_correction.pdf"))
plt.close()

print("✓ PCA plot saved: 01_PCA_before_correction.pdf")
print(f"  Variance explained: PC1={pca.explained_variance_ratio_[0]*100:.1f}%, PC2={pca.explained_variance_ratio_[1]*100:.1f}%")

# ============================================================================
# STEP 6: APPLY COMBAT BATCH CORRECTION
# ============================================================================
print("\nSTEP 6: Applying ComBat batch correction...")
print("  (This may take 1-2 minutes for ~22k probes)")

# Prepare batch vector (numeric: 0 = GSE25628, 1 = GSE7305)
batch_vector = np.where(sample_combined['batch'] == 'GSE25628', 0, 1)

# ComBat requires:
#   - data: probes × samples (expression matrix)
#   - batch: batch labels (numeric)
#   - covariates (optional): other variables to preserve (e.g., disease status)

# Create covariate matrix to preserve group effects
group_numeric = pd.Categorical(sample_combined['group']).codes
covariates_df = pd.DataFrame({
    'group': group_numeric
}, index=sample_combined.index)

# Apply ComBat
print("  Running ComBat...")

# The Combat class expects: data (samples × features), batch
# Also note: data should be samples × features, not features × samples!
# So we need to transpose

# Option 1: If Combat is a class with fit_transform
combat_model = Combat()
expr_corrected = combat_model.fit_transform(
    expr_merged.values.T,  # samples × probes (transposed!)
    batch_vector           # batch labels
).T  # Transpose back to probes × samples

# Convert back to DataFrame
expr_corrected = pd.DataFrame(
    expr_corrected,
    index=expr_merged.index,
    columns=expr_merged.columns
)


print("✓ ComBat batch correction completed")
print(f"  Output matrix: {expr_corrected.shape}")

# ============================================================================
# STEP 7: VISUALIZE BATCH EFFECT (AFTER CORRECTION)
# ============================================================================
print("\nSTEP 7: Visualizing batch effects (AFTER correction)...")

print("  Computing PCA...")
pca_after = PCA(n_components=2, random_state=42)
pca_coords_after = pca_after.fit_transform(expr_corrected.T)

pca_after_df = pd.DataFrame(
    pca_coords_after,
    columns=['PC1', 'PC2'],
    index=expr_corrected.columns
)
pca_after_df['batch'] = sample_combined['batch'].values
pca_after_df['group'] = sample_combined['group'].values

# Create comparison plot
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Top row: Before correction
batch_colors = {'GSE25628': 'red', 'GSE7305': 'blue'}
for batch in batch_colors:
    mask = pca_df['batch'] == batch
    axes[0, 0].scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], 
                      label=batch, alpha=0.7, s=100, 
                      color=batch_colors[batch])
axes[0, 0].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
axes[0, 0].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
axes[0, 0].set_title("BEFORE Correction (by dataset)")
axes[0, 0].legend()
axes[0, 0].grid(alpha=0.3)

group_colors = {'normal': 'blue', 'ectopic': 'red', 'diseased': 'orange'}
for group in group_colors:
    mask = pca_df['group'] == group
    if mask.sum() > 0:
        axes[0, 1].scatter(pca_df[mask]['PC1'], pca_df[mask]['PC2'], 
                          label=group, alpha=0.7, s=100, 
                          color=group_colors[group])
axes[0, 1].set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
axes[0, 1].set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
axes[0, 1].set_title("BEFORE Correction (by group)")
axes[0, 1].legend()
axes[0, 1].grid(alpha=0.3)

# Bottom row: After correction
for batch in batch_colors:
    mask = pca_after_df['batch'] == batch
    axes[1, 0].scatter(pca_after_df[mask]['PC1'], pca_after_df[mask]['PC2'], 
                      label=batch, alpha=0.7, s=100, 
                      color=batch_colors[batch])
axes[1, 0].set_xlabel(f"PC1 ({pca_after.explained_variance_ratio_[0]*100:.1f}%)")
axes[1, 0].set_ylabel(f"PC2 ({pca_after.explained_variance_ratio_[1]*100:.1f}%)")
axes[1, 0].set_title("AFTER Correction (by dataset)")
axes[1, 0].legend()
axes[1, 0].grid(alpha=0.3)

for group in group_colors:
    mask = pca_after_df['group'] == group
    if mask.sum() > 0:
        axes[1, 1].scatter(pca_after_df[mask]['PC1'], pca_after_df[mask]['PC2'], 
                          label=group, alpha=0.7, s=100, 
                          color=group_colors[group])
axes[1, 1].set_xlabel(f"PC1 ({pca_after.explained_variance_ratio_[0]*100:.1f}%)")
axes[1, 1].set_ylabel(f"PC2 ({pca_after.explained_variance_ratio_[1]*100:.1f}%)")
axes[1, 1].set_title("AFTER Correction (by group)")
axes[1, 1].legend()
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "02_PCA_before_after_correction.pdf"))
plt.close()

print("✓ Comparison PCA plot saved: 02_PCA_before_after_correction.pdf")

# ============================================================================
# STEP 8: SAVE HARMONIZED DATASET
# ============================================================================
print("\nSTEP 8: Saving harmonized dataset...")

# Save corrected expression matrix
harmonized_expr_file = os.path.join(DATA_PROCESSED, "harmonized_expression_matrix.csv")
expr_corrected.to_csv(harmonized_expr_file)
print(f"✓ Expression matrix saved: {harmonized_expr_file}")

# Save sample info with batch labels
harmonized_sample_file = os.path.join(DATA_PROCESSED, "harmonized_sample_info.csv")
sample_combined.to_csv(harmonized_sample_file)
print(f"✓ Sample info saved: {harmonized_sample_file}")

# Save batch-corrected RDS (for R compatibility)
print("  (Saving RDS version for R scripts...)")
try:
    import pyreadr
    pyreadr.write_rds(
        os.path.join(DATA_PROCESSED, "harmonized_expression_matrix.rds"),
        expr_corrected.T  # Transpose: samples × genes for R
    )
    print("✓ RDS file saved for R compatibility")
except ImportError:
    print("  (pyreadr not installed — skipping RDS. Use CSV instead.)")

# ============================================================================
# STEP 9: SUMMARY STATISTICS
# ============================================================================
print("\n" + "="*70)
print("PHASE 1 SUMMARY")
print("="*70)

print(f"\nDatasets harmonized:")
print(f"  GSE25628: {(sample_combined['batch'] == 'GSE25628').sum()} samples (normal, ectopic, eutopic)")
print(f"  GSE7305:  {(sample_combined['batch'] == 'GSE7305').sum()} samples (normal, diseased)")
print(f"\nCombined dataset:")
print(f"  Total samples: {expr_corrected.shape[1]}")
print(f"  Total probes: {expr_corrected.shape[0]}")
print(f"  Expression scale: log2 (after batch correction)")
print(f"\nSample group composition:")
print(sample_combined['group'].value_counts())
print(f"\nBatch effect assessment:")
print(f"  Before correction - PC1 variance: {pca.explained_variance_ratio_[0]*100:.1f}%")
print(f"  After correction  - PC1 variance: {pca_after.explained_variance_ratio_[0]*100:.1f}%")
if pca_after.explained_variance_ratio_[0] < pca.explained_variance_ratio_[0]:
    print("  ✓ Batch effect reduced (PC1 variance decreased)")
else:
    print("  ⚠ PC1 variance similar (batch effect may persist in first PC)")

print(f"\nOutput files created:")
print(f"  - harmonized_expression_matrix.csv ({expr_corrected.shape[0]} × {expr_corrected.shape[1]})")
print(f"  - harmonized_sample_info.csv")
print(f"  - PCA comparison plots in figures/harmonization/")

print(f"\n✓✓✓ PHASE 1 COMPLETE ✓✓✓")
print(f"Ready for Phase 2: DEG Analysis (use R + limma)")
