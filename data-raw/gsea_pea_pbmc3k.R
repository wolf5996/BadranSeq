# Code to prepare `gsea_pbmc3k_cd8t` dataset
#
# This script performs differential expression analysis (DEA) on the bundled
# pbmc3k dataset, comparing CD8 T cells against all other cell types, then
# runs a pathway enrichment analysis (PEA) using GSEA on the ranked gene list.
#
# The resulting gseaGO object is saved to data/ for use as example input
# for `do_GseaPlot()`.

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the bundled pbmc3k dataset ------------------------------------------------
load("data/pbmc3k.rda")

# DEA: Compare CD8 T cells vs all others ----------------------------------------
# Create a binary grouping: CD8 T vs all other cell types
pbmc3k$group <- ifelse(
  pbmc3k$seurat_annotations == "CD8 T",
  "CD8_T",
  "Other"
)
Idents(pbmc3k) <- "group"

# Run differential expression
markers <- FindMarkers(
  pbmc3k,
  ident.1 = "CD8_T",
  ident.2 = "Other",
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0,
  verbose = TRUE
)

# Prepare ranked gene list for GSEA -------------------------------------------
# clusterProfiler expects a named numeric vector sorted in descending order
# Use average log2FC as the ranking metric
ranked_genes <- markers$avg_log2FC
names(ranked_genes) <- rownames(markers)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Remove NA values if any
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

# Keep only genes that have Entrez IDs for gseGO--------------------------------
gene_entrez <- bitr(
  names(ranked_genes),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Create named vector using Entrez IDs
ranked_entrez <- ranked_genes[gene_entrez$SYMBOL]
names(ranked_entrez) <- gene_entrez$ENTREZID
ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

# Run GSEA with GO Biological Process -------------------------------------------
set.seed(42)

gsea_pbmc3k_cd8t <- gseGO(
  geneList = ranked_entrez,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Save to package data directory ------------------------------------------------
usethis::use_data(gsea_pbmc3k_cd8t, overwrite = TRUE, compress = "xz")
