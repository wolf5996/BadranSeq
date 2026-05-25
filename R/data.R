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
