# Code to prepare `gsea_pbmc3k` dataset
#
# This script performs differential expression analysis (DEA) on the bundled
# pbmc3k dataset, comparing Cluster 2 (effectively 100%% CD14+ Monocytes)
# against all other clusters, then runs a pathway enrichment analysis (PEA)
# using GSEA on the ranked gene list.
#
# The resulting gseaGO object is saved to data/ for use as example input
# for `do_GseaPlot()`.

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load the bundled pbmc3k dataset ------------------------------------------------
load("data/pbmc3k.rda")

# DEA: Compare Cluster 2 (CD14+ Mono) vs all other clusters ----------------------
# Cluster 2 is >97%% CD14+ Monocytes, making it the cleanest mono population
Seurat::Idents(pbmc3k) <- "seurat_clusters"

markers <- Seurat::FindMarkers(
  pbmc3k,
  ident.1 = "2",
  ident.2 = NULL,
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0,
  verbose = TRUE
)

# Prepare ranked gene list for GSEA -------------------------------------------
ranked_genes <- sort(markers$avg_log2FC, decreasing = TRUE)
names(ranked_genes) <- rownames(markers)
ranked_genes <- ranked_genes[!is.na(ranked_genes)]

# Cap extreme values (prevent NA issues with fgsea)
u <- stats::quantile(abs(ranked_genes), 0.99, na.rm = TRUE)
tmp <- pmin(pmax(ranked_genes, -u), u)
names(tmp) <- names(ranked_genes)
ranked_genes <- sort(tmp, decreasing = TRUE)

# Break exact ties with negligible random noise
set.seed(42)
noise <- stats::runif(length(ranked_genes), -1e-9, 1e-9)
ranked_genes <- ranked_genes + noise
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Keep only genes that have Entrez IDs for gseGO--------------------------------
gene_entrez <- clusterProfiler::bitr(
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

gsea_pbmc3k <- clusterProfiler::gseGO(
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
usethis::use_data(gsea_pbmc3k, overwrite = TRUE, compress = "xz")
