# R/convert.R
# Seurat to AnnData conversion

#' Convert a Seurat Object to AnnData
#'
#' @description
#' Converts a Seurat object to a Python AnnData object and optionally writes
#' it to an .h5ad file. The first element of \code{layers} becomes
#' \code{AnnData.X}; remaining layers are stored in \code{AnnData.layers}.
#' All reductions present in the Seurat object are transferred to
#' \code{AnnData.obsm} (e.g. PCA becomes \code{X_pca}, UMAP becomes
#' \code{X_umap}).
#'
#' Requires the \pkg{anndata} R package and the Python \code{anndata} package.
#'
#' @param seurat_obj Seurat object.
#' @param outfile character. Path to output .h5ad file (default: NULL = no file
#'   written, AnnData returned in memory only).
#' @param assay character. Assay to convert (default: "RNA").
#' @param layers character. Vector of layer names to include. The first
#'   element becomes AnnData.X; the rest are stored as AnnData layers
#'   (default: \code{c("counts", "data")}).
#'
#' @return The AnnData object (invisibly).
#'
#' @examples
#' \dontrun{
#' # Write to file
#' seurat_to_h5ad(obj, outfile = "data.h5ad")
#'
#' # Get in-memory AnnData (no file)
#' adata <- seurat_to_h5ad(obj, layers = "counts")
#'
#' # Custom assay and layers
#' seurat_to_h5ad(obj, outfile = "rna.h5ad",
#'                assay = "RNA", layers = c("counts", "data"))
#' }
#'
#' @export
seurat_to_h5ad <- function(
    seurat_obj,
    outfile = NULL,
    assay = "RNA",
    layers = c("counts", "data")
) {

  # --- Validate inputs ---

  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop(
      "Package 'anndata' is required for seurat_to_h5ad(). ",
      "Install with: install.packages('anndata')"
    )
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required for seurat_to_h5ad(). ",
      "Install with: install.packages('reticulate')"
    )
  }

  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop(
      "Assay '", assay, "' not found. Available: ",
      paste(SeuratObject::Assays(seurat_obj), collapse = ", ")
    )
  }

  available_layers <- SeuratObject::Layers(seurat_obj[[assay]])
  missing_layers <- setdiff(layers, available_layers)
  if (length(missing_layers) > 0) {
    stop(
      "Layer(s) not found in assay '", assay, "': ",
      paste(missing_layers, collapse = ", "),
      ". Available: ", paste(available_layers, collapse = ", ")
    )
  }

  # --- Build AnnData ---

  SeuratObject::DefaultAssay(seurat_obj) <- assay

  main <- layers[1]
  X <- t(as.matrix(SeuratObject::GetAssayData(seurat_obj, layer = main)))

  adata <- anndata::AnnData(
    X = X,
    obs = seurat_obj@meta.data,
    var = seurat_obj[[assay]][[]]
  )

  # Add remaining layers
  for (l in layers[-1]) {
    adata$layers[l] <- t(as.matrix(
      SeuratObject::GetAssayData(seurat_obj, layer = l)
    ))
  }

  # Transfer reductions
  for (red in names(seurat_obj@reductions)) {
    adata$obsm[[paste0("X_", red)]] <- SeuratObject::Embeddings(
      seurat_obj, reduction = red
    )
  }

  # --- Optionally write to file ---

  if (!is.null(outfile)) {
    adata$write_h5ad(outfile)
    message("Saved: ", outfile)
  }

  invisible(adata)
}
