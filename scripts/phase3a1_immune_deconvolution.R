# ============================================================================
# PHASE 3: PROBE-TO-GENE MAPPING (Using pre-installed dplyr)
# ============================================================================

library(hgu133plus2.db)
library(dplyr) # Bypassing tidyverse since dplyr is already installed!

# Set up paths
project_dir <- "C:/Users/Yasna/OneDrive/Belgeler/endometriosis-transcriptomic-analysis"
expr_path <- file.path(project_dir, "data", "processed", "harmonized_expression_matrix.csv")
output_path <- file.path(project_dir, "data", "processed", "gene_expression_matrix.csv")

# 1. Load your expression matrix
cat("Loading harmonized probe matrix...\n")
expr_matrix <- read.csv(expr_path, row.names = 1, check.names = FALSE)

# 2. Fetch probe-to-symbol mapping from the database
cat("Fetching gene annotations...\n")
annotations <- AnnotationDbi::select(hgu133plus2.db, 
                                     keys = rownames(expr_matrix), 
                                     columns = c("SYMBOL"), 
                                     keytype = "PROBEID")

# Drop any probes that don't map to a real gene symbol
annotations <- annotations %>% filter(!is.na(SYMBOL) & SYMBOL != "")

# 3. Filter and merge with expression data
expr_mapped <- expr_matrix[annotations$PROBEID, ]
expr_mapped$SYMBOL <- annotations$SYMBOL

# Calculate variance per probe to resolve duplicates elegantly
probe_variance <- apply(expr_matrix[annotations$PROBEID, ], 1, var)
expr_mapped$Variance <- probe_variance

# Collapse duplicates: Keep the probe with the highest variance per gene symbol
cat("Collapsing duplicate probes into unique Gene Symbols...\n")
expr_collapsed <- expr_mapped %>%
  arrange(SYMBOL, desc(Variance)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

# Set Gene Symbols as row names and clean up tracking columns
rownames(expr_collapsed) <- expr_collapsed$SYMBOL
expr_collapsed <- expr_collapsed %>% dplyr::select(-SYMBOL, -Variance)

# 4. Save the publication-ready gene expression matrix
write.csv(expr_collapsed, output_path)
cat("✓ Success! Saved gene-level matrix with", nrow(expr_collapsed), "unique genes to:\n", output_path, "\n")