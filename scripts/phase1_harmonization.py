#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
PHASE 1: HARMONIZATION (NO BATCH CORRECTION)
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

parser = argparse.ArgumentParser(description="Phase 1: Harmonization")
parser.add_argument("--project_root", default=os.getcwd())
args = parser.parse_args()

PROJECT_DIR = args.project_root
DATA_PROCESSED = os.path.join(PROJECT_DIR, "data", "processed")
FIGURES_DIR = os.path.join(PROJECT_DIR, "figures", "harmonization")
os.makedirs(FIGURES_DIR, exist_ok=True)

print("="*70)
print("PHASE 1: HARMONIZATION (NO BATCH CORRECTION)")
print("="*70)
print(f"Project directory: {PROJECT_DIR}")

# Load data
expr_25628 = pd.read_csv(os.path.join(DATA_PROCESSED, "GSE25628_expression_matrix.csv"), index_col=0)
expr_7305 = pd.read_csv(os.path.join(DATA_PROCESSED, "GSE7305_expression_matrix.csv"), index_col=0)
sample_25628 = pd.read_csv(os.path.join(DATA_PROCESSED, "GSE25628_sample_info.csv"), index_col=0)
sample_7305 = pd.read_csv(os.path.join(DATA_PROCESSED, "GSE7305_sample_info.csv"), index_col=0)

print(f"GSE25628: {expr_25628.shape[0]} probes × {expr_25628.shape[1]} samples")
print(f"GSE7305:  {expr_7305.shape[0]} probes × {expr_7305.shape[1]} samples")

# Align probes
common_probes = expr_25628.index.intersection(expr_7305.index)
expr_25628_aligned = expr_25628.loc[common_probes].copy()
expr_7305_aligned = expr_7305.loc[common_probes].copy()

print(f"Aligned to {len(common_probes)} common probes")

# Add batch labels
sample_25628['batch'] = 'GSE25628'
sample_7305['batch'] = 'GSE7305'
sample_combined = pd.concat([sample_25628, sample_7305])

# Merge
expr_merged = pd.concat([expr_25628_aligned, expr_7305_aligned], axis=1)
print(f"Merged: {expr_merged.shape[0]} probes × {expr_merged.shape[1]} samples")

# ==== CRITICAL: Check batch x group structure ====
crosstab = pd.crosstab(sample_combined['batch'], sample_combined['group'])
print("\nBatch × Group design:")
print(crosstab.to_string())

confounded_groups = [g for g in crosstab.columns if (crosstab[g] > 0).sum() < 2]
if confounded_groups:
    print(f"\nWARNING: These groups exist in only ONE batch: {confounded_groups}")
    print("   Cross-batch comparisons involving these groups are NOT valid.")
    print("   Use per-dataset analysis only.")

# Save merged data (without correction)
expr_merged.to_csv(os.path.join(DATA_PROCESSED, "harmonized_expression_matrix.csv"))
sample_combined.to_csv(os.path.join(DATA_PROCESSED, "harmonized_sample_info.csv"))

# Save confound flags
pd.DataFrame({
    "group": crosstab.columns,
    "n_batches_present": [(crosstab[g] > 0).sum() for g in crosstab.columns],
    "confounded_with_batch": [g in confounded_groups for g in crosstab.columns],
}).to_csv(os.path.join(DATA_PROCESSED, "batch_confound_flags.csv"), index=False)

# Save PCA before (for reference)
pca = PCA(n_components=2, random_state=42)
coords = pca.fit_transform(expr_merged.T)
var = pca.explained_variance_ratio_

print(f"\nPCA before any correction: PC1: {var[0]:.1%}, PC2: {var[1]:.1%}")

print("\n" + "="*70)
print("PHASE 1 COMPLETE — NO BATCH CORRECTION APPLIED")
print("="*70)
print("Outputs saved to:", DATA_PROCESSED)
print("\nProceed to Phase 2 for per-dataset differential expression.")