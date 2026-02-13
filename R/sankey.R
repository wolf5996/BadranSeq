# R/sankey.R
# Sankey (alluvial) plot for categorical metadata relationships

# Suppress R CMD check notes for NSE variables
utils::globalVariables(c("alluvium", "x", "stratum"))

#' Sankey Plot for Categorical Metadata
#'
#' @description
#' Creates a Sankey (alluvial) diagram showing the relationship between 2–3
#' categorical metadata variables from a Seurat object. Each stratum is
#' labelled and coloured; flows connect cells across variables.
#'
#' Requires the \pkg{ggalluvial} package to be installed.
#'
#' @param object Seurat object.
#' @param columns character. Vector of 2–3 metadata column names whose
#'   relationships to visualise (required).
#' @param label.size numeric. Text size for stratum labels (default: 3).
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # Relationship between cell type and cluster
#' do_SankeyPlot(seurat_obj,
#'               columns = c("seurat_annotations", "seurat_clusters"))
#'
#' # Three-way relationship
#' do_SankeyPlot(seurat_obj,
#'               columns = c("seurat_annotations", "condition", "Phase"))
#' }
#'
#' @export
do_SankeyPlot <- function(
    object,
    columns,
    label.size = 3
) {

  # --- Check ggalluvial availability ---

  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop(
      "Package 'ggalluvial' is required for do_SankeyPlot(). ",
      "Install it with: install.packages('ggalluvial')"
    )
  }

  # --- Validation ---

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  if (missing(columns) || !is.character(columns)) {
    stop("columns must be a character vector of 2-3 metadata column names")
  }

  if (length(columns) < 2 || length(columns) > 3) {
    stop(
      "columns must contain 2 or 3 metadata column names, got ",
      length(columns)
    )
  }

  available <- colnames(object@meta.data)
  missing_cols <- setdiff(columns, available)
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found in metadata: ",
      paste(missing_cols, collapse = ", "),
      ". Available: ", paste(available, collapse = ", ")
    )
  }

  # --- Fetch and prepare data ---

  df <- fetch_cell_data(object, metadata = columns)

  # Drop NAs
  complete <- stats::complete.cases(df[, columns, drop = FALSE])
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    message(n_dropped, " cell(s) with NA values excluded from plot.")
    df <- df[complete, ]
  }

  if (nrow(df) == 0) {
    stop("No cells remaining after removing NAs.")
  }

  # --- Aggregate and convert to lodes form ---

  counts <- dplyr::count(df, !!!rlang::syms(columns), name = "Freq")

  lodes <- counts |>
    ggalluvial::to_lodes_form(axes = seq_along(columns))

  ggplot2::ggplot(
    lodes,
    ggplot2::aes(
      x = .data[["x"]],
      stratum = .data[["stratum"]],
      alluvium = .data[["alluvium"]],
      y = .data[["Freq"]],
      fill = .data[["stratum"]],
      label = .data[["stratum"]]
    )
  ) +
    ggalluvial::geom_flow(alpha = 0.5) +
    ggalluvial::geom_stratum(alpha = 0.5) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      size = label.size
    ) +
    ggplot2::labs(x = "Metadata", y = "Frequency") +
    theme_badranseq() +
    ggplot2::theme(legend.position = "none")
}
