#!/usr/bin/env Rscript
# ============================================================================
# PROJECT 1: ENDOMETRIOSIS TRANSCRIPTOMICS
# PHASE 2: DIFFERENTIAL EXPRESSION ANALYSIS (limma)
# FIXED: Per-dataset analysis with pure normal controls
# ============================================================================

cat("Loading libraries...\n")
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

cat("Libraries loaded\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  PROJECT_DIR <- args[1]
} else {
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

set.seed(42)

cat("Output figures:", FIGURES_DIR, "\n")
cat("Output results:", RESULTS_DIR, "\n\n")

cat(strrep("=", 70), "\n")
cat("PHASE 2: DIFFERENTIAL EXPRESSION ANALYSIS (limma)\n")
cat("PER-DATASET ANALYSIS WITH PURE NORMAL CONTROLS\n")
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

expr_all <- read.csv(expr_file, row.names = 1, check.names = FALSE)
expr_matrix <- as.matrix(expr_all)

sample_info <- read.csv(sample_file, row.names = 1)

cat("Total samples:", ncol(expr_matrix), "\n")
cat("Groups:", paste(names(table(sample_info$group)), "=", 
                     table(sample_info$group), collapse=", "), "\n\n")

# ============================================================================
# STEP 2: PROBE-TO-GENE SYMBOL MAPPING (once, for all datasets)
# ============================================================================
cat("STEP 2: Mapping probes to gene symbols...\n")

probe_to_symbol <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr_matrix),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  distinct(PROBEID, .keep_all = TRUE)

symbol_lookup <- setNames(probe_to_symbol$SYMBOL, probe_to_symbol$PROBEID)
cat("Mapped", length(symbol_lookup), "of", nrow(expr_matrix), "probes to gene symbols\n\n")

# ============================================================================
# STEP 3: PER-DATASET DIFFERENTIAL EXPRESSION
# ============================================================================
cat("STEP 3: Running per-dataset differential expression...\n")

# Define datasets and their valid contrasts
datasets <- list(
  GSE25628 = list(
    samples = rownames(sample_info)[sample_info$batch == "GSE25628"],
    groups = c("ectopic", "eutopic", "normal"),
    contrasts = c("ectopic_vs_normal", "eutopic_vs_normal", "ectopic_vs_eutopic")
  ),
  GSE7305 = list(
    samples = rownames(sample_info)[sample_info$batch == "GSE7305"],
    groups = c("diseased", "normal"),
    contrasts = c("diseased_vs_normal")
  )
)

all_results <- list()
all_sig <- list()
all_fits <- list()           # FIX (#11): store each dataset's fit2 + its
all_contrast_cols <- list()  #  contrast names so MA plots can be generated
                             #  per contrast below (previously the MA function
                             #  was defined but never called, and fit2 was
                             #  overwritten each loop, so no MA plot was produced).

for (ds_name in names(datasets)) {
  cat("\n", paste(rep("-", 60), collapse = ""), "\n")
  cat("Processing:", ds_name, "\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Subset to this dataset
  ds_samples <- datasets[[ds_name]]$samples
  expr_ds <- expr_matrix[, ds_samples, drop = FALSE]
  sample_ds <- sample_info[ds_samples, , drop = FALSE]
  
  cat("  Samples:", ncol(expr_ds), "\n")
  cat("  Groups:", paste(names(table(sample_ds$group)), "=",
                         table(sample_ds$group), collapse=", "), "\n")
  
  # Make sure group is a factor
  sample_ds$group <- factor(sample_ds$group)
  
  # Create design matrix
  design <- model.matrix(~ 0 + group, data = sample_ds)
  colnames(design) <- gsub("group", "", colnames(design))
  
  cat("  Design columns:", paste(colnames(design), collapse=", "), "\n")
  
  # Fit model
  fit <- lmFit(expr_ds, design)
  
  # Define contrasts for this dataset
  groups_available <- levels(sample_ds$group)
  contrast_list <- list()
  
  if (ds_name == "GSE25628") {
    if ("ectopic" %in% groups_available && "normal" %in% groups_available) {
      contrast_list$ectopic_vs_normal <- "ectopic - normal"
    }
    if ("eutopic" %in% groups_available && "normal" %in% groups_available) {
      contrast_list$eutopic_vs_normal <- "eutopic - normal"
    }
    if ("ectopic" %in% groups_available && "eutopic" %in% groups_available) {
      contrast_list$ectopic_vs_eutopic <- "ectopic - eutopic"
    }
  } else if (ds_name == "GSE7305") {
    if ("diseased" %in% groups_available && "normal" %in% groups_available) {
      contrast_list$diseased_vs_normal <- "diseased - normal"
    }
  }
  
  if (length(contrast_list) == 0) {
    cat("  No valid contrasts for", ds_name, "\n")
    next
  }
  
  cat("  Contrasts:", paste(names(contrast_list), collapse=", "), "\n")
  
  contrast_matrix <- do.call(makeContrasts, c(contrast_list, list(levels = design)))
  
  # Apply eBayes
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  # FIX (#11): retain this dataset's fit + contrast names for MA plots later
  all_fits[[ds_name]] <- fit2
  all_contrast_cols[[ds_name]] <- colnames(contrast_matrix)
  
  # Extract results for each contrast
  for (contrast_name in colnames(contrast_matrix)) {
    results <- topTable(fit2, coef = contrast_name, number = Inf,
                        adjust.method = "BH", sort.by = "P")
    results$probe_id <- rownames(results)
    results$gene_symbol <- ifelse(
      results$probe_id %in% names(symbol_lookup),
      symbol_lookup[results$probe_id],
      results$probe_id
    )
    results <- results[, c("probe_id", "gene_symbol", "logFC", "AveExpr",
                           "t", "P.Value", "adj.P.Val")]
    
    # Save full results
    safe_name <- paste0(ds_name, "_", contrast_name)
    write.csv(results,
              file.path(RESULTS_DIR, paste0(safe_name, "_full_results.csv")),
              row.names = FALSE)
    
    # Save significant genes
    sig <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
    write.csv(sig,
              file.path(RESULTS_DIR, paste0(safe_name, "_significant.csv")),
              row.names = FALSE)
    
    cat("    ✓", contrast_name, ":", nrow(sig), "significant DEGs\n")
    
    all_results[[safe_name]] <- results
    all_sig[[safe_name]] <- sig
  }
}

# ============================================================================
# STEP 4: CROSS-DATASET COMPARISON (OVERLAP)
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("CROSS-DATASET DEG OVERLAP\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# Find overlapping DEGs between datasets
if ("GSE25628_ectopic_vs_normal" %in% names(all_sig) && 
    "GSE7305_diseased_vs_normal" %in% names(all_sig)) {
  
  genes_25628 <- all_sig$GSE25628_ectopic_vs_normal$probe_id
  genes_7305 <- all_sig$GSE7305_diseased_vs_normal$probe_id
  
  overlap <- intersect(genes_25628, genes_7305)
  
  cat("\nOverlap between GSE25628 (ectopic vs normal) and GSE7305 (diseased vs normal):\n")
  cat("  GSE25628 DEGs:", length(genes_25628), "\n")
  cat("  GSE7305 DEGs: ", length(genes_7305), "\n")
  cat("  Overlap:      ", length(overlap), "\n")
  cat("  Overlap proportion:", round(length(overlap)/min(length(genes_25628), length(genes_7305))*100, 1), "%\n")
  
  # Save overlap
  overlap_df <- data.frame(probe_id = overlap)
  write.csv(overlap_df,
            file.path(RESULTS_DIR, "cross_dataset_DEG_overlap.csv"),
            row.names = FALSE)
}

# ============================================================================
# STEP 5: GENERATE FIGURES (combined)
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("GENERATING FIGURES\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

# ---- Volcano plots ----
cat("\nCreating volcano plots...\n")

plot_volcano <- function(results, title, filename) {
  if (nrow(results) == 0) {
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
  
  top_genes <- results %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(P.Value) %>%
    head(10)
  
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, aes(label = gene_symbol),
                             size = 3, max.overlaps = 10)
  }
  
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  return(p)
}

for (name in names(all_results)) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", name)
  filename <- file.path(FIGURES_DIR, paste0("01_volcano_", safe_name, ".pdf"))
  plot_volcano(all_results[[name]], 
               paste(gsub("_", " ", name), "- Volcano Plot"),
               filename)
  cat("  Volcano plot:", name, "\n")
}

# ---- MA plots ----
# FIX (#11): the previous version defined plot_ma_combined() but never
# called it, and relied on a single leftover `fit2` (only the last
# dataset's). No MA plot was ever produced despite the console printing
# "Creating MA plots...". Below, each dataset's stored fit (all_fits) is
# used directly to produce one real, saved plot per contrast.
cat("\nCreating MA plots...\n")
ma_plots_created <- 0
for (ds_name in names(all_fits)) {
  ds_fit <- all_fits[[ds_name]]
  contrast_cols <- all_contrast_cols[[ds_name]]
  n_col <- min(length(contrast_cols), 2)
  n_row <- ceiling(length(contrast_cols) / n_col)
  pdf(file.path(FIGURES_DIR, paste0("02_MA_plots_", ds_name, ".pdf")),
      width = 6 * n_col, height = 6 * n_row)
  par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 2))
  for (i in seq_along(contrast_cols)) {
    limma::plotMA(ds_fit, coef = i,
                  main = paste(ds_name, "-", gsub("_", " ", contrast_cols[i])),
                  ylim = c(-6, 6))
    abline(h = c(-1, 1), col = "red", lty = 2, lwd = 2)
    ma_plots_created <- ma_plots_created + 1
  }
  dev.off()
  cat("  MA plots saved for", ds_name, "\n")
}
cat("  Total MA plots created:", ma_plots_created, "\n")

# ---- Heatmap ----
cat("\nCreating heatmap...\n")

# Use the first contrast with significant genes
sig_contrast <- NULL
for (name in names(all_sig)) {
  if (nrow(all_sig[[name]]) > 0) {
    sig_contrast <- name
    break
  }
}

if (!is.null(sig_contrast) && nrow(all_sig[[sig_contrast]]) > 0) {
  top_genes_df <- all_sig[[sig_contrast]] %>%
    arrange(P.Value) %>%
    head(50)
  top_genes <- top_genes_df$probe_id
  
  if (length(top_genes) > 0) {
    top_expr <- expr_matrix[top_genes, ]
    top_expr_scaled <- t(scale(t(top_expr)))
    
    row_labels <- top_genes_df$gene_symbol[match(rownames(top_expr_scaled), top_genes_df$probe_id)]
    dup_labels <- duplicated(row_labels) | duplicated(row_labels, fromLast = TRUE)
    row_labels[dup_labels] <- paste0(row_labels[dup_labels], " (",
                                     rownames(top_expr_scaled)[dup_labels], ")")
    rownames(top_expr_scaled) <- row_labels
    
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
    cat("  Heatmap saved\n")
  }
}

# ---- Expression Profiles ----
cat("\nCreating expression profiles...\n")

if (!is.null(sig_contrast) && nrow(all_sig[[sig_contrast]]) > 0) {
  top_up <- all_sig[[sig_contrast]] %>%
    filter(logFC > 0) %>%
    arrange(P.Value) %>%
    head(3)
  
  top_down <- all_sig[[sig_contrast]] %>%
    filter(logFC < 0) %>%
    arrange(P.Value) %>%
    head(3)
  
  top_genes_plot <- rbind(top_up, top_down)
  
  if (nrow(top_genes_plot) > 0) {
    expr_subset <- expr_matrix[top_genes_plot$probe_id, ]
    expr_melted <- melt(as.matrix(expr_subset))
    colnames(expr_melted) <- c("Probe", "Sample", "Expression")
    expr_melted$Group <- sample_info[expr_melted$Sample, "group"]
    
    symbol_map <- setNames(top_genes_plot$gene_symbol, top_genes_plot$probe_id)
    expr_melted$Gene <- symbol_map[as.character(expr_melted$Probe)]
    expr_melted$Gene <- factor(expr_melted$Gene, levels = unique(symbol_map[top_genes_plot$probe_id]))
    
    p <- ggplot(expr_melted, aes(x = Group, y = Expression, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
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
    cat("  Expression profiles saved\n")
  }
}

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PHASE 2 COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Method: limma (Linear Models for Microarray)\n")
cat("Correction: Benjamini-Hochberg FDR\n\n")

cat("Per-dataset differential expression completed:\n")
for (name in names(all_sig)) {
  cat(sprintf("  %s: %d DEGs\n", name, nrow(all_sig[[name]])))
}

cat("\nOutput files:\n")
cat(sprintf("  - Results: %s\n", RESULTS_DIR))
cat(sprintf("  - Figures: %s\n", FIGURES_DIR))

cat("\nReady for Phase 3: Gene Mapping & Immune Profiling ")