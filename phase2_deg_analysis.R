#!/usr/bin/env Rscript
# ============================================================================
# PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
# PHASE 2: DIFFERENTIAL EXPRESSION ANALYSIS (DEG)
#
# Workflow:
#   1. Load harmonized expression matrix from Phase 1
#   2. Create experimental design
#   3. Fit linear models using limma (gold standard for microarray)
#   4. Extract contrasts (pairwise comparisons between groups)
#   5. Apply statistical testing + FDR correction
#   6. Generate visualizations (volcano plots, heatmaps, MA plots)
#   7. Save DEG results for Phase 3 (immune profiling + pathways)
#
# Runtime: 10-15 minutes
# Output: DEG lists + publication-quality figures
# ============================================================================

# Load required libraries
cat("Loading libraries...\n")

library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)

cat("✓ Libraries loaded\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# UPDATE THIS PATH
PROJECT_DIR <- "C:/Users/Yasna/OneDrive/Belgeler/endometriosis-transcriptomic-analysis"

# Define directories
DATA_PROCESSED <- file.path(PROJECT_DIR, "data", "processed")
FIGURES_DIR <- file.path(PROJECT_DIR, "figures", "deg_analysis")
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "deg_analysis")

# Create output directories
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

cat(strrep("=", 70), "\n")
cat("PHASE 2: DIFFERENTIAL EXPRESSION ANALYSIS (limma)\n")
cat(strrep("=", 70), "\n")
cat("Project directory:", PROJECT_DIR, "\n")
cat("Output figures:", FIGURES_DIR, "\n")
cat("Output results:", RESULTS_DIR, "\n\n")

# ============================================================================
# STEP 1: LOAD HARMONIZED DATA FROM PHASE 1
# ============================================================================
cat("STEP 1: Loading harmonized data from Phase 1...\n")

expr_file <- file.path(DATA_PROCESSED, "harmonized_expression_matrix.csv")
sample_file <- file.path(DATA_PROCESSED, "harmonized_sample_info.csv")

if (!file.exists(expr_file)) {
  stop("ERROR: Harmonized expression matrix not found at:\n", expr_file)
}

expr_data <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_matrix <- as.matrix(expr_data)

cat("✓ Expression matrix loaded:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
cat("  Data range:", round(min(expr_matrix), 2), "to", round(max(expr_matrix), 2), "\n")

sample_info <- read.csv(sample_file, row.names = 1)
cat("✓ Sample info loaded:", nrow(sample_info), "samples\n")
cat("  Groups:", paste(names(table(sample_info$group)), "=", table(sample_info$group), collapse=", "), "\n")

if (!all(colnames(expr_matrix) == rownames(sample_info))) {
  stop("ERROR: Sample order mismatch!")
}

cat("✓ Sample order verified\n\n")

# ============================================================================
# STEP 2: CREATE EXPERIMENTAL DESIGN MATRIX
# ============================================================================
cat("STEP 2: Creating experimental design...\n")

sample_info$group <- factor(sample_info$group, 
                            levels = c("normal", "ectopic", "eutopic", "diseased"))

design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- gsub("group", "", colnames(design))

cat("✓ Design matrix created\n")
cat("Samples per group:\n")
print(table(sample_info$group))
cat("\n")

# ============================================================================
# STEP 3: FIT LINEAR MODELS (limma)
# ============================================================================
cat("STEP 3: Fitting linear models with limma...\n")

fit <- lmFit(expr_matrix, design)
cat("✓ Linear models fitted for", nrow(fit), "genes\n\n")

# ============================================================================
# STEP 4: DEFINE CONTRASTS
# ============================================================================
cat("STEP 4: Defining contrasts...\n")

contrast.matrix <- makeContrasts(
  ectopic_vs_normal = ectopic - normal,
  eutopic_vs_normal = eutopic - normal,
  ectopic_vs_eutopic = ectopic - eutopic,
  diseased_vs_normal = diseased - normal,
  levels = design
)

cat("✓ Contrasts defined (4 comparisons)\n\n")

fit2 <- contrasts.fit(fit, contrast.matrix)

# ============================================================================
# STEP 5: EMPIRICAL BAYES
# ============================================================================
cat("STEP 5: Applying empirical Bayes statistics...\n")

fit2 <- eBayes(fit2)
cat("✓ Empirical Bayes applied\n\n")

# ============================================================================
# STEP 6: EXTRACT RESULTS
# ============================================================================
cat("STEP 6: Extracting DEG results...\n")

extract_deg <- function(fit2_obj, contrast_name) {
  results <- topTable(fit2_obj, coef = contrast_name, number = Inf, 
                      adjust.method = "BH", sort.by = "P")
  results$gene_id <- rownames(results)
  results <- results[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  return(results)
}

results_list <- list(
  ectopic_vs_normal = extract_deg(fit2, "ectopic_vs_normal"),
  eutopic_vs_normal = extract_deg(fit2, "eutopic_vs_normal"),
  ectopic_vs_eutopic = extract_deg(fit2, "ectopic_vs_eutopic"),
  diseased_vs_normal = extract_deg(fit2, "diseased_vs_normal")
)

cat("✓ Results extracted\n\n")

# ============================================================================
# STEP 7: IDENTIFY SIGNIFICANT GENES
# ============================================================================
cat("STEP 7: Identifying significant genes (FDR < 0.05, |logFC| > 1)...\n")

fdr_threshold <- 0.05
logfc_threshold <- 1

sig_genes <- lapply(results_list, function(res) {
  res[res$adj.P.Val < fdr_threshold & abs(res$logFC) > logfc_threshold, ]
})

for (name in names(sig_genes)) {
  n_sig <- nrow(sig_genes[[name]])
  n_up <- sum(sig_genes[[name]]$logFC > 0)
  n_down <- sum(sig_genes[[name]]$logFC < 0)
  cat(sprintf("  %s: %d (%d up, %d down)\n", name, n_sig, n_up, n_down))
}

cat("\n")

# ============================================================================
# STEP 8: SAVE RESULTS
# ============================================================================
cat("STEP 8: Saving results...\n")

for (name in names(results_list)) {
  write.csv(results_list[[name]], 
            file.path(RESULTS_DIR, paste0(name, "_full_results.csv")), 
            row.names = FALSE)
}

for (name in names(sig_genes)) {
  write.csv(sig_genes[[name]], 
            file.path(RESULTS_DIR, paste0(name, "_significant.csv")), 
            row.names = FALSE)
}

cat("✓ Results saved to results/deg_analysis/\n\n")

# ============================================================================
# STEP 9: VOLCANO PLOTS
# ============================================================================
cat("STEP 9: Creating volcano plots...\n")

plot_volcano <- function(results, title, filename) {
  plot_data <- results %>%
    mutate(
      significant = adj.P.Val < 0.05 & abs(logFC) > 1,
      direction = ifelse(logFC > 0, "Up", "Down"),
      direction = ifelse(!significant, "NS", direction)
    )
  
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = direction), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(title = title, x = "log2(FC)", y = "-log10(FDR)", color = "Regulation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  top_genes <- results %>% filter(adj.P.Val < 0.05) %>% arrange(P.Value) %>% head(5)
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene_id), size = 3)
  }
  
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
}

plot_volcano(results_list$ectopic_vs_normal, "Ectopic vs. Normal",
             file.path(FIGURES_DIR, "01_volcano_ectopic_vs_normal.pdf"))
cat("✓ Volcano plot 1\n")

plot_volcano(results_list$eutopic_vs_normal, "Eutopic vs. Normal",
             file.path(FIGURES_DIR, "02_volcano_eutopic_vs_normal.pdf"))
cat("✓ Volcano plot 2\n")

plot_volcano(results_list$ectopic_vs_eutopic, "Ectopic vs. Eutopic",
             file.path(FIGURES_DIR, "03_volcano_ectopic_vs_eutopic.pdf"))
cat("✓ Volcano plot 3\n")

plot_volcano(results_list$diseased_vs_normal, "Diseased vs. Normal",
             file.path(FIGURES_DIR, "04_volcano_diseased_vs_normal.pdf"))
cat("✓ Volcano plot 4\n\n")

# ============================================================================
# STEP 10: MA PLOTS
# ============================================================================
cat("STEP 10: Creating MA plots...\n")

pdf(file.path(FIGURES_DIR, "05_MA_plots.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2))

for (i in 1:4) {
  limma::plotMA(fit2, coef = i, main = names(results_list)[i], ylim = c(-6, 6))
  abline(h = c(-1, 1), col = "gray", lty = 2)
}

dev.off()
cat("✓ MA plots saved\n\n")

# ============================================================================
# STEP 11: HEATMAP
# ============================================================================
cat("STEP 11: Creating heatmap...\n")

top_genes_names <- sig_genes$ectopic_vs_normal %>%
  arrange(P.Value) %>%
  head(50) %>%
  pull(gene_id)

if (length(top_genes_names) > 0) {
  top_expr <- expr_matrix[top_genes_names, ]
  top_expr_scaled <- t(scale(t(top_expr)))
  
  pdf(file.path(FIGURES_DIR, "06_heatmap_top_genes.pdf"), width = 10, height = 12)
  
  sample_groups <- factor(sample_info$group, levels = c("normal", "ectopic", "eutopic", "diseased"))
  group_colors <- c(normal = "blue", ectopic = "red", eutopic = "orange", diseased = "gold")
  
  heatmap(top_expr_scaled, 
          ColSideColors = as.character(group_colors[sample_groups]),
          scale = "none",
          main = "Top 50 DEGs: Ectopic vs. Normal",
          margins = c(10, 4))
  
  legend("topright", legend = names(group_colors), col = group_colors, lty = 1, lwd = 3)
  
  dev.off()
  cat("✓ Heatmap saved\n\n")
}

# ============================================================================
# SUMMARY
# ============================================================================
cat(strrep("=", 70), "\n")
cat("PHASE 2 COMPLETE\n")
cat(strrep("=", 70), "\n")

cat("Method: limma (Linear Models for Microarray)\n")
cat("Correction: Benjamini-Hochberg FDR\n\n")

cat("Contrasts analyzed:\n")
cat("  1. Ectopic vs. Normal\n")
cat("  2. Eutopic vs. Normal\n")
cat("  3. Ectopic vs. Eutopic\n")
cat("  4. Diseased vs. Normal\n\n")

cat("Significant genes (FDR < 0.05, |logFC| > 1):\n")
for (name in names(sig_genes)) {
  cat(sprintf("  %s: %d genes\n", name, nrow(sig_genes[[name]])))
}

cat("\n✓✓✓ Ready for Phase 3: Immune profiling & Pathways ✓✓✓\n")

