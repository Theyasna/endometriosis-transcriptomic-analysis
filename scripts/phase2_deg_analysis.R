#!/usr/bin/env Rscript
# ============================================================================
# PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
# PHASE 2: DIFFERENTIAL EXPRESSION ANALYSIS (DEG) - FULLY FIXED VERSION
# ============================================================================

cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(reshape2)
})

cat("✓ Libraries loaded\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

# Use current directory as default if no argument provided
if (length(args) > 0) {
  PROJECT_DIR <- args[1]
} else {
  # Auto-detect: if running from scripts folder, go up one level
  if (file.exists("../data/processed")) {
    PROJECT_DIR <- ".."
  } else {
    PROJECT_DIR <- getwd()
  }
}

cat("Project directory:", PROJECT_DIR, "\n")

DATA_PROCESSED <- file.path(PROJECT_DIR, "data", "processed")
FIGURES_DIR <- file.path(PROJECT_DIR, "figures", "deg_analysis")
RESULTS_DIR <- file.path(PROJECT_DIR, "results", "deg_analysis")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)  # Reproducibility

cat("Output figures:", FIGURES_DIR, "\n")
cat("Output results:", RESULTS_DIR, "\n\n")

cat(strrep("=", 70), "\n")
cat("PHASE 2: DIFFERENTIAL EXPRESSION ANALYSIS (limma)\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# STEP 1: LOAD HARMONIZED DATA
# ============================================================================
cat("STEP 1: Loading harmonized data from Phase 1...\n")

expr_file <- file.path(DATA_PROCESSED, "harmonized_expression_matrix.csv")
sample_file <- file.path(DATA_PROCESSED, "harmonized_sample_info.csv")

if (!file.exists(expr_file)) {
  stop(paste("ERROR: Harmonized expression matrix not found at:", expr_file))
}
if (!file.exists(sample_file)) {
  stop(paste("ERROR: Harmonized sample info not found at:", sample_file))
}

expr_data <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_matrix <- as.matrix(expr_data)

cat("✓ Expression matrix loaded:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "samples\n")
cat("  Data range:", round(min(expr_matrix, na.rm = TRUE), 2), "to", 
    round(max(expr_matrix, na.rm = TRUE), 2), "\n")

sample_info <- read.csv(sample_file, row.names = 1)
cat("✓ Sample info loaded:", nrow(sample_info), "samples\n")

# Check group distribution
cat("  Groups:", paste(names(table(sample_info$group)), "=", 
                       table(sample_info$group), collapse=", "), "\n")

# Verify sample order
if (!all(colnames(expr_matrix) == rownames(sample_info))) {
  cat("WARNING: Sample order mismatch! Aligning...\n")
  common_samples <- intersect(colnames(expr_matrix), rownames(sample_info))
  expr_matrix <- expr_matrix[, common_samples]
  sample_info <- sample_info[common_samples, ]
  cat("  Aligned to", length(common_samples), "samples\n")
}
cat("✓ Sample order verified\n\n")

# ============================================================================
# STEP 2: DESIGN MATRIX
# ============================================================================
cat("STEP 2: Creating experimental design...\n")

# Define group levels (all possible groups in the data)
all_groups <- unique(sample_info$group)
cat("  Groups found:", paste(all_groups, collapse=", "), "\n")

# Set factor levels
sample_info$group <- factor(sample_info$group, levels = all_groups)

# Create design matrix
design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- gsub("group", "", colnames(design))

cat("✓ Design matrix created\n")
cat("  Design columns:", paste(colnames(design), collapse=", "), "\n")
cat("  Samples per group:\n")
print(table(sample_info$group))
cat("\n")

# ============================================================================
# STEP 3: FIT LINEAR MODELS
# ============================================================================
cat("STEP 3: Fitting linear models with limma...\n")

# Check for samples with complete data
if (ncol(expr_matrix) < 3) {
  stop("ERROR: Not enough samples for differential expression analysis")
}

fit <- lmFit(expr_matrix, design)
cat("✓ Linear models fitted for", nrow(fit), "probes\n\n")

# ============================================================================
# STEP 4: CONTRASTS
# ============================================================================
cat("STEP 4: Defining contrasts...\n")

# Build contrast matrix dynamically based on available groups
contrast_list <- list()

# Ectopic vs Normal (if both exist)
if ("ectopic" %in% all_groups && "normal" %in% all_groups) {
  contrast_list$ectopic_vs_normal <- "ectopic - normal"
}

# Eutopic vs Normal (if both exist)
if ("eutopic" %in% all_groups && "normal" %in% all_groups) {
  contrast_list$eutopic_vs_normal <- "eutopic - normal"
}

# Ectopic vs Eutopic (if both exist)
if ("ectopic" %in% all_groups && "eutopic" %in% all_groups) {
  contrast_list$ectopic_vs_eutopic <- "ectopic - eutopic"
}

# Diseased vs Normal (if both exist)
if ("diseased" %in% all_groups && "normal" %in% all_groups) {
  contrast_list$diseased_vs_normal <- "diseased - normal"
}

# Additional contrasts if applicable
if ("ectopic" %in% all_groups && "diseased" %in% all_groups) {
  contrast_list$ectopic_vs_diseased <- "ectopic - diseased"
}

if ("eutopic" %in% all_groups && "diseased" %in% all_groups) {
  contrast_list$eutopic_vs_diseased <- "eutopic - diseased"
}

cat("  Defining", length(contrast_list), "contrasts:\n")
for (name in names(contrast_list)) {
  cat("    -", name, ":", contrast_list[[name]], "\n")
}

if (length(contrast_list) == 0) {
  stop("ERROR: No valid contrasts can be defined with available groups")
}

contrast_matrix <- do.call(makeContrasts, c(contrast_list, list(levels = design)))
cat("✓ Contrast matrix created with", ncol(contrast_matrix), "contrasts\n\n")

# ============================================================================
# STEP 5: EMPIRICAL BAYES
# ============================================================================
cat("STEP 5: Applying empirical Bayes statistics...\n")
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
cat("✓ Empirical Bayes applied\n\n")

# ============================================================================
# STEP 6: EXTRACT RESULTS
# ============================================================================
cat("STEP 6: Extracting DEG results...\n")

extract_deg <- function(fit2_obj, contrast_name) {
  results <- topTable(fit2_obj, coef = contrast_name, number = Inf,
                      adjust.method = "BH", sort.by = "P")
  results$probe_id <- rownames(results)
  results <- results[, c("probe_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  return(results)
}

results_list <- list()
for (contrast_name in colnames(contrast_matrix)) {
  results_list[[contrast_name]] <- extract_deg(fit2, contrast_name)
}

cat("✓ Results extracted for", length(results_list), "contrasts\n\n")

# ============================================================================
# STEP 7: SIGNIFICANT GENES
# ============================================================================
cat("STEP 7: Identifying significant probes (FDR < 0.05, |logFC| > 1)...\n")

fdr_threshold <- 0.05
logfc_threshold <- 1

sig_genes <- list()
for (contrast_name in names(results_list)) {
  res <- results_list[[contrast_name]]
  sig <- res[res$adj.P.Val < fdr_threshold & abs(res$logFC) > logfc_threshold, ]
  sig_genes[[contrast_name]] <- sig
}

cat("  Significant probes per contrast:\n")
for (name in names(sig_genes)) {
  n_sig <- nrow(sig_genes[[name]])
  if (n_sig > 0) {
    n_up <- sum(sig_genes[[name]]$logFC > 0)
    n_down <- sum(sig_genes[[name]]$logFC < 0)
    cat(sprintf("    %s: %d (%d up, %d down)\n", name, n_sig, n_up, n_down))
  } else {
    cat(sprintf("    %s: 0\n", name))
  }
}
cat("\n")

# ============================================================================
# STEP 8: SAVE RESULTS
# ============================================================================
cat("STEP 8: Saving results...\n")

# Save full results
for (name in names(results_list)) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", name)
  write.csv(results_list[[name]],
            file.path(RESULTS_DIR, paste0(safe_name, "_full_results.csv")),
            row.names = FALSE)
}

# Save significant genes
for (name in names(sig_genes)) {
  if (nrow(sig_genes[[name]]) > 0) {
    safe_name <- gsub("[^A-Za-z0-9_]", "_", name)
    write.csv(sig_genes[[name]],
              file.path(RESULTS_DIR, paste0(safe_name, "_significant.csv")),
              row.names = FALSE)
  }
}

cat("✓ Results saved to:", RESULTS_DIR, "\n\n")

# ============================================================================
# STEP 9: VOLCANO PLOTS
# ============================================================================
cat("STEP 9: Creating volcano plots...\n")

plot_volcano <- function(results, title, filename) {
  if (nrow(results) == 0) {
    cat("  Warning: No results to plot for", title, "\n")
    return(NULL)
  }
  
  plot_data <- results %>%
    mutate(
      significant = adj.P.Val < 0.05 & abs(logFC) > 1,
      direction = case_when(
        !significant ~ "NS",
        logFC > 0 ~ "Up",
        logFC < 0 ~ "Down"
      )
    )
  
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = direction), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Up" = "#e74c3c", "Down" = "#2ecc71", "NS" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5, color = "darkred") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5, color = "darkblue") +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(FDR)", color = "Regulation") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  # Add top genes labels
  top_genes <- results %>% 
    filter(adj.P.Val < 0.05) %>% 
    arrange(P.Value) %>% 
    head(10)
  
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = probe_id), 
                             size = 3, max.overlaps = 10)
  }
  
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  return(p)
}

volcano_count <- 0
for (name in names(results_list)) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", name)
  filename <- file.path(FIGURES_DIR, paste0("01_volcano_", safe_name, ".pdf"))
  plot_volcano(results_list[[name]], 
               paste(gsub("_", " ", name), "- Volcano Plot"),
               filename)
  volcano_count <- volcano_count + 1
  cat("  ✓ Volcano plot", volcano_count, ":", name, "\n")
}
cat("\n")

# ============================================================================
# STEP 10: MA PLOTS
# ============================================================================
cat("STEP 10: Creating MA plots...\n")

if (ncol(fit2) > 0) {
  pdf(file.path(FIGURES_DIR, "02_MA_plots.pdf"), width = 12, height = 10)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  n_plots <- min(ncol(fit2), 4)
  for (i in 1:n_plots) {
    contrast_name <- colnames(fit2)[i]
    limma::plotMA(fit2, coef = i, 
                  main = paste(gsub("_", " ", contrast_name)), 
                  ylim = c(-6, 6))
    abline(h = c(-1, 1), col = "red", lty = 2, lwd = 2)
  }
  dev.off()
  cat("✓ MA plots saved (showing first", n_plots, "contrasts)\n\n")
} else {
  cat("  Warning: No contrasts available for MA plots\n\n")
}

# ============================================================================
# STEP 11: HEATMAP OF TOP GENES
# ============================================================================
cat("STEP 11: Creating heatmap of top DEGs...\n")

# Use first contrast with significant genes
sig_contrast <- NULL
for (name in names(sig_genes)) {
  if (nrow(sig_genes[[name]]) > 0) {
    sig_contrast <- name
    break
  }
}

if (!is.null(sig_contrast) && nrow(sig_genes[[sig_contrast]]) > 0) {
  top_genes <- sig_genes[[sig_contrast]] %>%
    arrange(P.Value) %>%
    head(50) %>%
    pull(probe_id)
  
  if (length(top_genes) > 0) {
    top_expr <- expr_matrix[top_genes, ]
    
    # Scale expression
    top_expr_scaled <- t(scale(t(top_expr)))
    
    # Define group colors
    group_colors <- c(
      normal = "#2ecc71",
      ectopic = "#e74c3c",
      eutopic = "#3498db",
      diseased = "#f1c40f"
    )
    
    sample_groups <- factor(sample_info$group, 
                            levels = c("normal", "ectopic", "eutopic", "diseased"))
    
    pdf(file.path(FIGURES_DIR, "03_heatmap_top_genes.pdf"), width = 10, height = 12)
    
    heatmap(top_expr_scaled,
            ColSideColors = as.character(group_colors[sample_groups]),
            scale = "none",
            main = paste("Top 50 DEGs:", gsub("_", " ", sig_contrast)),
            margins = c(10, 6),
            cexRow = 0.6,
            cexCol = 0.8)
    
    legend("topright", 
           legend = names(group_colors), 
           fill = group_colors, 
           title = "Group",
           cex = 0.8)
    
    dev.off()
    cat("✓ Heatmap saved\n\n")
  } else {
    cat("  No top genes found for heatmap\n\n")
  }
} else {
  cat("  No significant genes found for heatmap\n\n")
}

# ============================================================================
# STEP 12: EXPRESSION PROFILES FOR TOP GENES
# ============================================================================
cat("STEP 12: Creating expression profiles for top genes...\n")

if (!is.null(sig_contrast) && nrow(sig_genes[[sig_contrast]]) > 0) {
  # Get top 6 genes (up and down)
  top_up <- sig_genes[[sig_contrast]] %>%
    filter(logFC > 0) %>%
    arrange(P.Value) %>%
    head(3)
  
  top_down <- sig_genes[[sig_contrast]] %>%
    filter(logFC < 0) %>%
    arrange(P.Value) %>%
    head(3)
  
  top_genes_plot <- rbind(top_up, top_down)
  
  if (nrow(top_genes_plot) > 0) {
    # Prepare data for plotting
    expr_subset <- expr_matrix[top_genes_plot$probe_id, ]
    expr_melted <- melt(as.matrix(expr_subset))
    colnames(expr_melted) <- c("Probe", "Sample", "Expression")
    expr_melted$Group <- sample_info[expr_melted$Sample, "group"]
    expr_melted$Probe <- factor(expr_melted$Probe, levels = top_genes_plot$probe_id)
    
    p <- ggplot(expr_melted, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      facet_wrap(~ Probe, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c(normal = "#2ecc71", ectopic = "#e74c3c",
                                   eutopic = "#3498db", diseased = "#f1c40f")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
      ) +
      labs(title = paste("Top DEGs:", gsub("_", " ", sig_contrast)),
           x = "", y = "Expression (log2)")
    
    ggsave(file.path(FIGURES_DIR, "04_top_genes_expression_profiles.pdf"),
           p, width = 12, height = 8, dpi = 300)
    cat("✓ Expression profiles saved\n\n")
  }
} else {
  cat("  No significant genes found for expression profiles\n\n")
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
for (i in 1:length(names(results_list))) {
  cat(sprintf("  %d. %s\n", i, names(results_list)[i]))
}

cat("\nSignificant genes (FDR < 0.05, |logFC| > 1):\n")
for (name in names(sig_genes)) {
  if (nrow(sig_genes[[name]]) > 0) {
    cat(sprintf("  %s: %d genes\n", name, nrow(sig_genes[[name]])))
  } else {
    cat(sprintf("  %s: 0 genes\n", name))
  }
}

cat("\nOutput files:\n")
cat(sprintf("  - Results: %s\n", RESULTS_DIR))
cat(sprintf("  - Figures: %s\n", FIGURES_DIR))
cat("\n✓✓✓ Ready for Phase 3: Gene Mapping & Immune Profiling ✓✓✓\n")