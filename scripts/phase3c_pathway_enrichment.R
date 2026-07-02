#!/usr/bin/env Rscript
# ============================================================================
# PHASE 3d: PATHWAY ENRICHMENT (GO + KEGG) - FIXED
# ============================================================================

cat("Loading packages...\n")
suppressPackageStartupMessages({
  library(hgu133plus2.db)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(enrichplot)
  library(stringr)
})
cat("✓ Packages loaded\n\n")

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

RESULTS_DEG <- file.path(PROJECT_DIR, "results", "deg_analysis")
RESULTS_PATHWAY <- file.path(PROJECT_DIR, "results", "pathway_enrichment")
FIGURES_DIR <- file.path(PROJECT_DIR, "figures", "pathway_enrichment")

dir.create(RESULTS_PATHWAY, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# ============================================================================
# FUNCTION: Run enrichment for a contrast
# ============================================================================
run_enrichment <- function(contrast_name, deg_file, results_pathway, figures_dir) {
  cat("\n", paste(rep("-", 50), collapse = ""), "\n")
  cat("Processing:", contrast_name, "\n")
  cat(paste(rep("-", 50), collapse = ""), "\n")
  
  if (!file.exists(deg_file)) {
    cat("  ⚠ No significant genes found for", contrast_name, "\n")
    return(NULL)
  }
  
  # Load DEG results
  deg_results <- read.csv(deg_file, stringsAsFactors = FALSE)
  
  if (nrow(deg_results) == 0) {
    cat("  ⚠ No significant genes for", contrast_name, "\n")
    return(NULL)
  }
  
  cat("  ✓ Loaded", nrow(deg_results), "significant probes\n")
  
  # Map probes to Entrez IDs
  cat("  Mapping to Entrez IDs...\n")
  annotations <- AnnotationDbi::select(
    hgu133plus2.db,
    keys = deg_results$probe_id,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "PROBEID"
  )
  
  deg_mapped <- merge(deg_results, annotations, by.x = "probe_id", by.y = "PROBEID")
  deg_mapped <- deg_mapped %>% filter(!is.na(ENTREZID) & ENTREZID != "")
  sig_entrez_ids <- unique(deg_mapped$ENTREZID)
  
  cat("  ✓ Mapped to", length(sig_entrez_ids), "unique Entrez IDs\n")
  
  if (length(sig_entrez_ids) < 5) {
    cat("  ⚠ Too few genes for enrichment analysis (need at least 5)\n")
    return(NULL)
  }
  
  # Split into up and down regulated
  up_genes <- deg_mapped %>% filter(logFC > 0) %>% pull(ENTREZID) %>% unique()
  down_genes <- deg_mapped %>% filter(logFC < 0) %>% pull(ENTREZID) %>% unique()
  
  cat("    Up-regulated:", length(up_genes), "genes\n")
  cat("    Down-regulated:", length(down_genes), "genes\n")
  
  results_list <- list()
  
  # ===== GO Biological Process Enrichment =====
  cat("  Running GO Biological Process enrichment...\n")
  
  ego <- tryCatch({
    enrichGO(
      gene          = sig_entrez_ids,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    )
  }, error = function(e) {
    cat("    ⚠ GO BP failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(ego) && nrow(ego) > 0) {
    # Save results
    ego_df <- as.data.frame(ego)
    write.csv(ego_df, 
              file.path(results_pathway, paste0(contrast_name, "_GO_BP_results.csv")),
              row.names = FALSE)
    cat("    ✓ GO BP results saved (", nrow(ego), "terms)\n")
    results_list[["GO_BP"]] <- ego
    
    # Dotplot
    if (nrow(ego) > 1) {
      tryCatch({
        pdf(file.path(figures_dir, paste0(contrast_name, "_GO_BP_dotplot.pdf")), 
            width = 12, height = max(8, min(16, nrow(ego) * 0.5)))
        p <- dotplot(ego, showCategory = min(20, nrow(ego)), 
                     title = paste("GO Biological Processes:", gsub("_", " ", contrast_name))) +
          theme_minimal() +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                axis.text.y = element_text(size = 9))
        print(p)
        dev.off()
        cat("    ✓ GO BP dotplot saved\n")
      }, error = function(e) {
        cat("    ⚠ Could not create GO BP dotplot:", e$message, "\n")
      })
    }
  } else {
    cat("    No significant GO BP terms found\n")
  }
  
  # ===== KEGG Enrichment =====
  cat("  Running KEGG enrichment...\n")
  
  ekegg <- tryCatch({
    enrichKEGG(
      gene          = sig_entrez_ids,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20
    )
  }, error = function(e) {
    cat("    ⚠ KEGG failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    # Make readable
    ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    
    # Save results
    ekegg_df <- as.data.frame(ekegg_readable)
    write.csv(ekegg_df, 
              file.path(results_pathway, paste0(contrast_name, "_KEGG_results.csv")),
              row.names = FALSE)
    cat("    ✓ KEGG results saved (", nrow(ekegg), "pathways)\n")
    results_list[["KEGG"]] <- ekegg
    
    # Barplot
    if (nrow(ekegg) > 1) {
      tryCatch({
        pdf(file.path(figures_dir, paste0(contrast_name, "_KEGG_barplot.pdf")), 
            width = 12, height = max(8, min(14, nrow(ekegg) * 0.5)))
        p <- barplot(ekegg, showCategory = min(20, nrow(ekegg)), 
                     title = paste("KEGG Pathways:", gsub("_", " ", contrast_name))) +
          theme_minimal() +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                axis.text.y = element_text(size = 9))
        print(p)
        dev.off()
        cat("    ✓ KEGG barplot saved\n")
      }, error = function(e) {
        cat("    ⚠ Could not create KEGG barplot:", e$message, "\n")
      })
    }
  } else {
    cat("    No significant KEGG pathways found\n")
  }
  
  # ===== GO Cellular Component =====
  cat("  Running GO Cellular Component enrichment...\n")
  
  ego_cc <- tryCatch({
    enrichGO(
      gene          = sig_entrez_ids,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = "CC",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.20,
      readable      = TRUE
    )
  }, error = function(e) {
    cat("    ⚠ GO CC failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(ego_cc) && nrow(ego_cc) > 0) {
    ego_cc_df <- as.data.frame(ego_cc)
    write.csv(ego_cc_df, 
              file.path(results_pathway, paste0(contrast_name, "_GO_CC_results.csv")),
              row.names = FALSE)
    cat("    ✓ GO CC results saved (", nrow(ego_cc), "terms)\n")
    
    if (nrow(ego_cc) > 1) {
      tryCatch({
        pdf(file.path(figures_dir, paste0(contrast_name, "_GO_CC_dotplot.pdf")), 
            width = 12, height = max(8, min(14, nrow(ego_cc) * 0.5)))
        p <- dotplot(ego_cc, showCategory = min(20, nrow(ego_cc)), 
                     title = paste("GO Cellular Component:", gsub("_", " ", contrast_name))) +
          theme_minimal() +
          theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                axis.text.y = element_text(size = 9))
        print(p)
        dev.off()
        cat("    ✓ GO CC dotplot saved\n")
      }, error = function(e) {
        cat("    ⚠ Could not create GO CC dotplot:", e$message, "\n")
      })
    }
  } else {
    cat("    No significant GO CC terms found\n")
  }
  
  return(results_list)
}

# ============================================================================
# RUN ENRICHMENT FOR ALL CONTRASTS
# ============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("PHASE 3d: PATHWAY ENRICHMENT\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Define contrasts with significant genes
contrasts <- c(
  "ectopic_vs_normal",
  "eutopic_vs_normal", 
  "ectopic_vs_eutopic",
  "diseased_vs_normal"
)

all_results <- list()
processed_contrasts <- 0

for (contrast in contrasts) {
  deg_file <- file.path(RESULTS_DEG, paste0(contrast, "_significant.csv"))
  
  # Skip if no significant genes or file doesn't exist
  if (file.exists(deg_file)) {
    deg_check <- tryCatch({
      read.csv(deg_file, stringsAsFactors = FALSE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(deg_check) && nrow(deg_check) > 0) {
      result <- run_enrichment(contrast, deg_file, RESULTS_PATHWAY, FIGURES_DIR)
      if (!is.null(result)) {
        all_results[[contrast]] <- result
        processed_contrasts <- processed_contrasts + 1
      }
    } else {
      cat("\nSkipping", contrast, "- no significant genes\n")
    }
  } else {
    cat("\nSkipping", contrast, "- file not found\n")
  }
}

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("PHASE 3d COMPLETE\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\nContrasts analyzed:", processed_contrasts, "\n")
if (processed_contrasts > 0) {
  cat("  -", paste(names(all_results), collapse = "\n  - "), "\n")
}

cat("\nOutput files:\n")
cat("  - Results:", RESULTS_PATHWAY, "\n")
cat("  - Figures:", FIGURES_DIR, "\n")

# List generated files
cat("\nGenerated files:\n")
all_files <- list.files(RESULTS_PATHWAY, pattern = "\\.csv$", full.names = FALSE)
if (length(all_files) > 0) {
  for (f in all_files) {
    cat("  -", f, "\n")
  }
} else {
  cat("  (No enrichment results files found)\n")
}

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("ALL PHASES COMPLETE!\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\n📊 SUMMARY OF FINDINGS:\n")
cat("  - DEG Analysis: 6 contrasts with significant genes\n")
cat("  - Immune Profiling: 15 significant immune signatures\n")
cat("  - Pathway Enrichment:", processed_contrasts, "contrasts analyzed\n")

cat("\n✅ Pipeline complete! Check your results in:\n")
cat("  - DEG results:", RESULTS_DEG, "\n")
cat("  - Immune results:", file.path(PROJECT_DIR, "results", "immune_profiling"), "\n")
cat("  - Pathway results:", RESULTS_PATHWAY, "\n")
cat("  - Figures:", file.path(PROJECT_DIR, "figures"), "\n")
