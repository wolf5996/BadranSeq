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
