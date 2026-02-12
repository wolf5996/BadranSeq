# R/sankey.R
# Sankey (alluvial) plot for categorical metadata relationships

# Suppress R CMD check notes for NSE variables
utils::globalVariables(c("alluvium", "x", "stratum", "n", "fill_var"))

#' Sankey Plot for Categorical Metadata
#'
#' @description
#' Creates a Sankey (alluvial) diagram showing the relationship between 2–3
#' categorical metadata variables from a Seurat object. Flows are coloured by
#' the first variable in \code{columns}. Cell counts are displayed inside
#' each stratum by default.
#'
#' Requires the \pkg{ggalluvial} package to be installed.
#'
#' @param object Seurat object.
#' @param columns character. Vector of 2–3 metadata column names whose
#'   relationships to visualise (required).
#' @param colors.use character. Named vector of colours for the first column
#'   in \code{columns}. If NULL, generates colours automatically with
#'   \code{\link{generate_badranseq_colors}}.
#' @param label logical. Whether to display cell counts inside strata
#'   (default: TRUE).
#' @param label.size numeric. Text size for stratum labels (default: 3).
#' @param alpha numeric. Flow transparency (default: 0.5).
#' @param curve_type character. Curve type for flows (default: "sigmoid").
#'   Passed to \code{ggalluvial::geom_alluvium()}.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # Relationship between cell type and cluster
#' do_SankeyPlot(seurat_obj, columns = c("seurat_annotations", "seurat_clusters"))
#'
#' # Three-way relationship
#' do_SankeyPlot(seurat_obj,
#'               columns = c("seurat_annotations", "condition", "Phase"))
#'
#' # Custom colours
#' do_SankeyPlot(seurat_obj,
#'               columns = c("seurat_annotations", "seurat_clusters"),
#'               colors.use = c("B" = "blue", "T" = "red"))
#' }
#'
#' @importFrom rlang .data syms
#' @export
do_SankeyPlot <- function(
    object,
    columns,
    colors.use = NULL,
    label = TRUE,
    label.size = 3,
    alpha = 0.5,
    curve_type = "sigmoid"
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
    stop("columns must contain 2 or 3 metadata column names, got ", length(columns))
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

  # --- Fetch data ---

  df <- fetch_cell_data(object, metadata = columns)

  # --- Handle NAs ---

  complete <- stats::complete.cases(df[, columns, drop = FALSE])
  n_dropped <- sum(!complete)
  if (n_dropped > 0) {
    message(n_dropped, " cell(s) with NA values excluded from plot.")
    df <- df[complete, ]
  }

  if (nrow(df) == 0) {
    stop("No cells remaining after removing NAs.")
  }

  # --- Convert to factors and drop unused levels ---

  for (col in columns) {
    df[[col]] <- factor(df[[col]])
  }

  # --- Aggregate counts ---

  counts <- dplyr::count(df, !!!rlang::syms(columns), name = "n")

  # --- Convert to lodes (long) form ---

  lodes <- ggalluvial::to_lodes_form(
    counts,
    axes = seq_along(columns),
    id = "alluvium"
  )

  # --- Generate colours for first variable ---

  fill_col <- columns[1]
  fill_levels <- levels(df[[fill_col]])
  n_levels <- length(fill_levels)

  if (is.null(colors.use)) {
    colors.use <- generate_badranseq_colors(n_levels)
    names(colors.use) <- fill_levels
  }

  # --- Map fill variable from first axis onto all lodes rows ---

  fill_map <- counts[[fill_col]]
  names(fill_map) <- seq_len(nrow(counts))
  lodes$fill_var <- factor(
    fill_map[as.character(lodes$alluvium)],
    levels = fill_levels
  )

  # --- Build plot ---

  p <- ggplot2::ggplot(
    lodes,
    ggplot2::aes(
      x = .data[["x"]],
      stratum = .data[["stratum"]],
      alluvium = .data[["alluvium"]],
      y = .data[["n"]],
      fill = .data[["fill_var"]]
    )
  )

  p <- p +
    ggalluvial::geom_alluvium(alpha = alpha, curve_type = curve_type) +
    ggalluvial::geom_stratum(width = 1/3, color = "black", fill = "white") +
    ggplot2::scale_fill_manual(values = colors.use, name = fill_col)

  # --- Stratum labels ---

  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data[["stratum"]]),
      stat = ggalluvial::StatStratum,
      size = label.size
    )
  }

  # --- Theme and axes ---

  p <- p +
    ggplot2::scale_x_discrete(expand = c(0.1, 0.05)) +
    ggplot2::labs(x = NULL, y = "Cells") +
    theme_badranseq()

  p
}
