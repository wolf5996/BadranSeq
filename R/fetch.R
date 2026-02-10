# R/fetch.R
# Tidy data extraction from Seurat objects

# Suppress R CMD check notes for dplyr/tidyr NSE variables
utils::globalVariables(c("feature", "cell_id", "name", "value"))

#' Build Cell Scaffold Tibble
#'
#' @description
#' Internal helper that constructs a base tibble with cell_id as a column,
#' optional embedding coordinates, and optional metadata from a Seurat object.
#'
#' @param object Seurat object.
#' @param reductions character. Reduction names to include (NULL = none).
#' @param dims integer. Dimensions to extract per reduction (default: 1:2).
#' @param metadata logical or character. TRUE = all metadata, FALSE = none,
#'   character vector = specific columns.
#'
#' @return tibble with cell_id and requested columns.
#' @keywords internal
.build_cell_scaffold <- function(
    object,
    reductions = NULL,
    dims = 1:2,
    metadata = TRUE
) {

  # Validation ----------

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  scaffold <- tibble::tibble(cell_id = colnames(object))

  # Embeddings ----------

  if (!is.null(reductions)) {
    for (red in reductions) {
      if (!red %in% names(object@reductions)) {
        stop(
          "'", red, "' not found. Available reductions: ",
          paste(names(object@reductions), collapse = ", ")
        )
      }

      embedding_mat <- object@reductions[[red]]@cell.embeddings
      max_dims <- ncol(embedding_mat)

      if (max(dims) > max_dims) {
        stop(
          "dims ", max(dims), " exceeds available dimensions (",
          max_dims, ") for reduction '", red, "'"
        )
      }

      embedding_df <- tibble::as_tibble(
        embedding_mat[, dims, drop = FALSE],
        .name_repair = "minimal"
      )

      # Column naming: PCA -> PC1, PC2; UMAP -> UMAP1, UMAP2; etc.
      prefix <- if (tolower(red) == "pca") "PC" else toupper(red)
      colnames(embedding_df) <- paste0(prefix, dims)

      scaffold <- dplyr::bind_cols(scaffold, embedding_df)
    }
  }

  # Metadata ----------

  if (is.logical(metadata) && isTRUE(metadata)) {
    meta_df <- tibble::as_tibble(object@meta.data, .name_repair = "minimal")
    scaffold <- dplyr::bind_cols(scaffold, meta_df)
  } else if (is.character(metadata)) {
    missing_cols <- setdiff(metadata, colnames(object@meta.data))
    if (length(missing_cols) > 0) {
      stop(
        "Metadata columns not found in metadata: ",
        paste(missing_cols, collapse = ", "),
        ". Available: ", paste(colnames(object@meta.data), collapse = ", ")
      )
    }
    meta_df <- tibble::as_tibble(
      object@meta.data[, metadata, drop = FALSE],
      .name_repair = "minimal"
    )
    scaffold <- dplyr::bind_cols(scaffold, meta_df)
  }
  # metadata = FALSE: do nothing, just cell_id

  scaffold
}

#' Fetch Cell-Level Data from a Seurat Object
#'
#' @description
#' Extracts a tidy tibble with one row per cell containing cell identifiers,
#' optional dimensionality reduction embeddings, and optional metadata.
#' No expression data is included — use \code{\link{fetch_feature_data}} for that.
#'
#' @param object Seurat object.
#' @param reductions character. Reduction names to include (default: NULL = none).
#'   Example: \code{c("umap", "pca")}.
#' @param dims integer. Dimensions to extract per reduction (default: 1:2).
#' @param metadata logical or character. TRUE = all metadata columns,
#'   FALSE = none, or a character vector of specific column names.
#'
#' @return A tibble with columns: cell_id, optional embedding columns, and
#'   optional metadata columns.
#'
#' @examples
#' \dontrun{
#' # All metadata, no embeddings
#' fetch_cell_data(seurat_obj)
#'
#' # With UMAP coordinates
#' fetch_cell_data(seurat_obj, reductions = "umap")
#'
#' # Specific metadata + PCA
#' fetch_cell_data(seurat_obj, reductions = "pca", dims = 1:10,
#'                 metadata = c("seurat_clusters", "condition"))
#' }
#'
#' @export
fetch_cell_data <- function(
    object,
    reductions = NULL,
    dims = 1:2,
    metadata = TRUE
) {
  .build_cell_scaffold(
    object = object,
    reductions = reductions,
    dims = dims,
    metadata = metadata
  )
}
