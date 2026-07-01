#!/usr/bin/env Rscript
# ============================================================================
# PLACE THIS AT THE TOP OF YOUR SCRIPT WITH THE LIBRARIES
# ============================================================================
cat("Loading local annotation and enrichment packages...\n")
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)

# Network fixes for Windows SSL connection errors
options(download.file.method = "wininet")
Sys.setenv(CURL_SSL_BACKEND = "openssl")

cat("✓ All packages loaded and network protocols optimized.\n\n")
cat("Loading local annotation and enrichment packages...\n")
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)

cat("✓ All packages loaded successfully.\n\n")

# ============================================================================
# CONFIGURATION & PATHS
# ============================================================================
PROJECT_DIR <- "C:/Users/Yasna/OneDrive/Belgeler/endometriosis-transcriptomic-analysis"

RESULTS_DEG     <- file.path(PROJECT_DIR, "results", "deg_analysis")
RESULTS_PATHWAY <- file.path(PROJECT_DIR, "results", "pathway_enrichment")
FIGURES_DIR     <- file.path(PROJECT_DIR, "figures", "pathway_enrichment")

dir.create(RESULTS_PATHWAY, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# STEP 1: LOAD SIGNIFICANT DEGs
# ============================================================================
cat("STEP 1: Loading significant DEG results from Phase 2...\n")
deg_file <- file.path(RESULTS_DEG, "ectopic_vs_normal_significant.csv")

if (!file.exists(deg_file)) {
  stop("ERROR: Significant DEG file not found. Please check your Phase 2 path.")
}

deg_results <- read.csv(deg_file, stringsAsFactors = FALSE)
cat("✓ Loaded", nrow(deg_results), "statistically significant differential probes.\n")

# ============================================================================
# STEP 2: LOCAL PROBE-TO-ENTREZ MAPPING (No more placeholders!)
# ============================================================================
cat("\nSTEP 2: Mapping probe IDs to official Entrez IDs using hgu133plus2.db...\n")

# Query your local platform database directly
annotations <- AnnotationDbi::select(
  hgu133plus2.db,
  keys    = deg_results$gene_id,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)

# Merge back with fold-change data from your analysis
deg_mapped <- merge(deg_results, annotations, by.x = "gene_id", by.y = "PROBEID")

# Drop any unmapped technical probes
deg_mapped <- deg_mapped %>% filter(!is.na(ENTREZID) & ENTREZID != "")
sig_entrez_ids <- unique(deg_mapped$ENTREZID)

cat("✓ Successfully mapped significant probes to", length(sig_entrez_ids), "unique Entrez Gene IDs.\n")

# ============================================================================
# STEP 3: RUN GENE ONTOLOGY (GO) ENRICHMENT Analysis
# ============================================================================
cat("\nSTEP 3: Running authentic GO Over-Representation Analysis (Biological Process)...\n")

ego <- enrichGO(
  gene          = sig_entrez_ids,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # Focus on Biological Processes
  pAdjustMethod = "BH",            # Benjamini-Hochberg FDR correction
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE             # Automatically translates numeric Entrez IDs back to clean Gene Symbols
)

if (!is.null(ego) && nrow(ego) > 0) {
  write.csv(as.data.frame(ego), file.path(RESULTS_PATHWAY, "go_bp_enrichment_results.csv"), row.names = FALSE)
  cat("✓ GO Analysis complete! Significant terms saved.\n")
} else {
  cat("⚠ No GO terms passed the significance thresholds.\n")
}

# ============================================================================
# STEP 4: RUN KEGG PATHWAY ENRICHMENT (With Fail-Safe Network Handling)
# ============================================================================
cat("\nSTEP 4: Running authentic KEGG Pathway Enrichment Analysis...\n")

ekegg <- NULL  # Initialize as NULL

tryCatch({
  ekegg <- enrichKEGG(
    gene          = sig_entrez_ids,
    organism      = "hsa",           # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )
  
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(ekegg), file.path(RESULTS_PATHWAY, "kegg_enrichment_results.csv"), row.names = FALSE)
    cat("✓ KEGG Analysis complete! Enriched pathways saved.\n")
  } else {
    cat("⚠ No KEGG pathways passed the significance thresholds.\n")
  }
}, error = function(e) {
  cat("⚠ KEGG online connection failed due to SSL/Network constraints.\n")
  cat("Message:", e$message, "\n")
  cat("Skipping KEGG step to preserve pipeline continuity. GO results remain secure.\n")
})

# ============================================================================
# STEP 5: GENERATE HIGH-RESOLUTION VISUALIZATIONS
# ============================================================================
cat("\nSTEP 5: Generating publication-ready graphics...\n")

# 1. GO Bubble Dotplot (Always runs since GO is entirely offline and local!)
if (!is.null(ego) && nrow(ego) > 0) {
  pdf(file.path(FIGURES_DIR, "01_go_enrichment_dotplot.pdf"), width = 11, height = 8)
  p1 <- dotplot(ego, showCategory = 15, title = "Top Enriched GO Biological Processes") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10, face = "bold")
    )
  print(p1)
  dev.off()
  cat("✓ Saved: 01_go_enrichment_dotplot.pdf\n")
}

# 2. KEGG Category Barplot (Only runs if KEGG successfully connected)
if (!is.null(ekegg) && nrow(ekegg) > 0) {
  pdf(file.path(FIGURES_DIR, "02_kegg_enrichment_barplot.pdf"), width = 11, height = 8)
  p2 <- barplot(ekegg, showCategory = 15, title = "Significantly Enriched KEGG Pathways") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y = element_text(size = 10, face = "bold")
    )
  print(p2)
  dev.off()
  cat("✓ Saved: 02_kegg_enrichment_barplot.pdf\n")
}