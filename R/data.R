#' PBMC 3K Dataset
#'
#' @description
#' A subset of 3000 peripheral blood mononuclear cells (PBMCs) from a
#' healthy donor, preprocessed as a Seurat object for demonstration purposes.
#'
#' @format A Seurat object containing:
#' \describe{
#'   \item{assays}{RNA assay with count and data layers}
#'   \item{reductions}{PCA and UMAP dimensionality reductions}
#'   \item{meta.data}{Cell metadata including cluster assignments}
#' }
#'
#' @source 10x Genomics PBMC dataset, processed with standard Seurat workflow
#'
#' @examples
#' data(pbmc3k)
#' pbmc3k
"pbmc3k"

#' DEA Results: CD14+ Monocytes vs Others (PBMC 3K)
#'
#' @description
#' Wilcoxon differential expression results comparing Cluster 2
#' (>97% CD14+ Monocytes) against all other clusters in the `pbmc3k`
#' dataset. Generated via `Seurat::FindMarkers()` with `logfc.threshold = 0`
#' so all genes are retained for downstream GSEA ranking.
#'
#' @format A data.frame with 2530 rows (all detected genes) and 5 columns:
#' \describe{
#'   \item{p_val}{Raw p-value from Wilcoxon rank-sum test}
#'   \item{avg_log2FC}{Average log2 fold-change (Mono vs Others)}
#'   \item{pct.1}{Percentage of cells in Cluster 2 expressing the gene}
#'   \item{pct.2}{Percentage of cells in all other clusters expressing the gene}
#'   \item{p_val_adj}{BH-adjusted p-value}
#' }
#' Row names are HGNC gene symbols.
#'
#' @source `Seurat::FindMarkers(ident.1 = "2")` on `pbmc3k` dataset.
#'
#' @examples
#' data(dea_pbmc3k_mono)
#' head(dea_pbmc3k_mono)
"dea_pbmc3k_mono"

#' GSEA Results: CD14+ Monocytes vs Others (PBMC 3K)
#'
#' @description
#' GSEA results from comparing Cluster 2 (>97% CD14+ Monocytes) against
#' all other clusters in the `pbmc3k` dataset. Generated via
#' `clusterProfiler::gseGO()` using biological process ontology on a ranked
#' gene list from Wilcoxon differential expression analysis.
#'
#' @format A `gseaResult` object from `clusterProfiler` with slots:
#' \describe{
#'   \item{result}{data.frame with columns: ID, Description, setSize,
#'     enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge,
#'     core_enrichment}
#' }
#'
#' @source Differential expression (Cluster 2 vs others) on `pbmc3k`
#'   dataset, GSEA via `clusterProfiler::gseGO(ont = "BP")`.
#'
#' @examples
#' data(gsea_pbmc3k)
#' gsea_pbmc3k
"gsea_pbmc3k"
