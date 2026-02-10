# R/fetch.R
# Tidy data extraction from Seurat objects

# Suppress R CMD check notes for dplyr/tidyr NSE variables
utils::globalVariables(c("feature"))

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

#' Fetch Feature Expression Data from a Seurat Object
#'
#' @description
#' Extracts a tidy tibble with one row per cell per feature. Expression values
#' are provided as one column per requested layer. Optionally includes
#' embedding coordinates and metadata, joined by cell_id.
#'
#' @param object Seurat object.
#' @param features character. Gene names to extract (required). Becomes a factor
#'   column preserving the input order.
#' @param layer character. Layer(s) to extract expression from.
#'   One of "counts", "data", "scale.data", or a combination
#'   (default: "data").
#' @param assay character. Assay to use (default: NULL = DefaultAssay).
#' @param reductions character. Reduction names to include (default: NULL = none).
#' @param dims integer. Dimensions to extract per reduction (default: 1:2).
#' @param metadata logical or character. TRUE = all metadata columns,
#'   FALSE = none, or a character vector of specific column names.
#'
#' @return A tibble with columns: cell_id, feature (factor), one column per layer,
#'   optional embedding columns, optional metadata columns.
#'
#' @examples
#' \dontrun{
#' # Basic usage -- normalized expression
#' fetch_feature_data(seurat_obj, features = c("CD3D", "CD8A"))
#'
#' # Multiple layers
#' fetch_feature_data(seurat_obj, features = "CD3D",
#'                    layer = c("counts", "data"))
#'
#' # With embeddings and specific metadata
#' fetch_feature_data(seurat_obj, features = c("CD3D", "CD14"),
#'                    reductions = "umap",
#'                    metadata = c("seurat_clusters", "condition"))
#'
#' # Explicit assay
#' fetch_feature_data(seurat_obj, features = "CD3D",
#'                    assay = "RNA", layer = "counts")
#' }
#'
#' @importFrom rlang %||%
#' @export
fetch_feature_data <- function(
    object,
    features,
    layer = "data",
    assay = NULL,
    reductions = NULL,
    dims = 1:2,
    metadata = TRUE
) {

  # Validation ----------

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  if (missing(features) || length(features) == 0 || !is.character(features)) {
    stop("features must be a non-empty character vector")
  }

  assay_name <- assay %||% Seurat::DefaultAssay(object)

  if (!assay_name %in% names(object@assays)) {
    stop(
      "Assay '", assay_name, "' not found. Available assays: ",
      paste(names(object@assays), collapse = ", ")
    )
  }

  available_layers <- SeuratObject::Layers(object[[assay_name]])
  missing_layers <- setdiff(layer, available_layers)
  if (length(missing_layers) > 0) {
    stop(
      "Layer(s) not found in assay '", assay_name, "': ",
      paste(missing_layers, collapse = ", "),
      ". Available: ", paste(available_layers, collapse = ", ")
    )
  }

  # Check features exist in the assay ----------

  assay_features <- rownames(object[[assay_name]])
  valid_features <- intersect(features, assay_features)
  missing_features <- setdiff(features, assay_features)

  if (length(missing_features) > 0) {
    warning(
      "Features not found in assay '", assay_name, "': ",
      paste(missing_features, collapse = ", "), ". Skipping."
    )
  }

  if (length(valid_features) == 0) {
    stop("No valid features found in assay '", assay_name, "'.")
  }

  # Build scaffold ----------

  scaffold <- .build_cell_scaffold(
    object = object,
    reductions = reductions,
    dims = dims,
    metadata = metadata
  )

  # Extract expression per layer ----------

  expression_list <- purrr::map(layer, function(lyr) {
    mat <- SeuratObject::LayerData(
      object = object,
      assay = assay_name,
      layer = lyr,
      features = valid_features
    )
    # mat is features x cells (possibly sparse) -- convert to dense then transpose
    dense_mat <- as.matrix(mat)
    df <- tibble::as_tibble(
      t(dense_mat),
      .name_repair = "minimal"
    )
    df <- dplyr::mutate(df, cell_id = colnames(object))

    # Pivot features long: one row per cell per feature
    tidyr::pivot_longer(
      df,
      cols = dplyr::all_of(valid_features),
      names_to = "feature",
      values_to = lyr
    )
  })

  # Join layers together ----------

  expr_tbl <- expression_list[[1]]

  if (length(expression_list) > 1) {
    for (i in 2:length(expression_list)) {
      expr_tbl <- dplyr::left_join(
        expr_tbl,
        expression_list[[i]],
        by = c("cell_id", "feature")
      )
    }
  }

  # Convert feature to factor preserving input order ----------

  expr_tbl <- dplyr::mutate(
    expr_tbl,
    feature = factor(feature, levels = valid_features)
  )

  # Join with scaffold ----------

  result <- dplyr::left_join(expr_tbl, scaffold, by = "cell_id")

  # Reorder: cell_id, feature, layers, embeddings, metadata ----------

  layer_cols <- intersect(layer, colnames(result))
  other_cols <- setdiff(
    colnames(result),
    c("cell_id", "feature", layer_cols)
  )
  result <- dplyr::select(
    result,
    dplyr::all_of(c("cell_id", "feature", layer_cols, other_cols))
  )

  result
}
