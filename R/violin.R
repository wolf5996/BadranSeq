# R/violin.R
# Violin plot with ggbetweenstats-style aesthetics

# Suppress R CMD check notes for dplyr/ggplot2 NSE variables
utils::globalVariables(c("centrality_label"))

# ============================================================================
# Internal: Build a single violin panel
# ============================================================================

#' Build a Single Violin Panel
#'
#' @description
#' Internal helper that constructs a single violin + boxplot + jitter + centrality
#' panel for one gene, using ggbetweenstats-style aesthetics.
#'
#' @param df data.frame. Long-format data with at least the y_col and group_col.
#' @param y_col character. Column name for the y-axis (expression values).
#' @param group_col character. Column name for the grouping variable.
#' @param title character. Plot title (typically the gene name).
#' @param colors character. Named vector of colors for groups.
#' @param pt.size numeric. Jitter point size.
#' @param pt.alpha numeric. Jitter point alpha.
#' @param violin.alpha numeric. Violin fill alpha.
#' @param violin.width numeric. Violin width.
#' @param boxplot logical. Whether to overlay a boxplot.
#' @param boxplot.width numeric. Boxplot width.
#' @param boxplot.alpha numeric. Boxplot alpha.
#' @param centrality logical. Whether to show centrality point and label.
#' @param centrality.fn character. Centrality function name ("median" or "mean").
#'
#' @return A ggplot2 object.
#' @keywords internal
.build_violin_panel <- function(
    df,
    y_col,
    group_col,
    title,
    colors,
    pt.size = 3,
    pt.alpha = 0.4,
    violin.alpha = 0.2,
    violin.width = 0.5,
    boxplot = TRUE,
    boxplot.width = 0.3,
    boxplot.alpha = 0.2,
    centrality = TRUE,
    centrality.fn = "median"
) {

  # Base plot
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[group_col]],
      y = .data[[y_col]],
      fill = .data[[group_col]],
      color = .data[[group_col]]
    )
  )

  # Layer 1: Violin

  p <- p + ggplot2::geom_violin(
    alpha = violin.alpha,
    width = violin.width,
    color = "black",
    linewidth = 0.5
  )

  # Layer 2: Boxplot overlay (optional)
  if (boxplot) {
    p <- p + ggplot2::geom_boxplot(
      width = boxplot.width,
      alpha = boxplot.alpha,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      show.legend = FALSE
    )
  }

  # Layer 3: Jitter (raw points)
  p <- p + ggplot2::geom_jitter(
    size = pt.size,
    alpha = pt.alpha,
    width = 0.15,
    show.legend = FALSE
  )

  # Layer 4: Centrality (optional)
  if (centrality) {
    fn <- switch(centrality.fn,
      "median" = stats::median,
      "mean"   = base::mean,
      stats::median
    )
    centrality_df <- dplyr::summarise(
      dplyr::group_by(df, .data[[group_col]]),
      centrality_label = round(fn(.data[[y_col]], na.rm = TRUE), 2),
      .groups = "drop"
    )

    # Diamond point for centrality
    p <- p + ggplot2::geom_point(
      data = centrality_df,
      ggplot2::aes(
        x = .data[[group_col]],
        y = .data[["centrality_label"]]
      ),
      shape = 18,
      size = 5,
      color = "black",
      inherit.aes = FALSE
    )

    # Label with numeric value
    p <- p + ggplot2::geom_label(
      data = centrality_df,
      ggplot2::aes(
        x = .data[[group_col]],
        y = .data[["centrality_label"]],
        label = .data[["centrality_label"]]
      ),
      fill = "white",
      color = "black",
      size = 3,
      fontface = "bold",
      nudge_y = diff(range(df[[y_col]], na.rm = TRUE)) * 0.05,
      inherit.aes = FALSE
    )
  }

  # Scales and theme
  p <- p +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(title = title, x = NULL, y = y_col) +
    theme_badranseq() +
    ggplot2::theme(legend.position = "none")

  p
}

# ============================================================================
# Exported: do_ViolinPlot
# ============================================================================

#' Violin Plot with ggbetweenstats-Style Aesthetics
#'
#' @description
#' Creates publication-ready violin plots for one or more features from a
#' Seurat object, with optional boxplot overlay, jittered raw data points,
#' and centrality markers (median/mean diamond with numeric label).
#'
#' For a single feature, returns a bare ggplot2 object. For multiple features,
#' returns a patchwork composition.
#'
#' @param object Seurat object.
#' @param features character. Gene names to plot (required). Must be non-empty.
#' @param group.by character. Metadata column for grouping. If NULL (default),
#'   uses active Idents.
#' @param layer character. Layer to extract expression from (default: "data").
#' @param assay character. Assay to use (default: NULL = DefaultAssay).
#' @param colors.use character. Named vector of colors for groups. If NULL,
#'   generates colors automatically with \code{\link{generate_badranseq_colors}}.
#' @param pt.size numeric. Jitter point size (default: 3).
#' @param pt.alpha numeric. Jitter point transparency (default: 0.4).
#' @param violin.alpha numeric. Violin fill alpha (default: 0.2).
#' @param violin.width numeric. Violin width (default: 0.5).
#' @param boxplot logical. Whether to overlay a boxplot (default: TRUE).
#' @param boxplot.width numeric. Boxplot width (default: 0.3).
#' @param boxplot.alpha numeric. Boxplot alpha (default: 0.2).
#' @param centrality logical. Whether to show centrality point and label
#'   (default: TRUE).
#' @param centrality.fn character. Centrality function ("median" or "mean",
#'   default: "median").
#' @param ncol numeric. Number of columns for patchwork layout when plotting
#'   multiple features (default: NULL = auto).
#'
#' @return A ggplot2 object (single feature) or patchwork object (multiple
#'   features).
#'
#' @examples
#' \dontrun{
#' # Single gene violin plot
#' do_ViolinPlot(seurat_obj, features = "CD3D")
#'
#' # Multiple genes
#' do_ViolinPlot(seurat_obj, features = c("CD3D", "CD8A", "CD14"))
#'
#' # Custom grouping and colors
#' do_ViolinPlot(seurat_obj, features = "CD3D",
#'               group.by = "seurat_clusters",
#'               colors.use = c("0" = "red", "1" = "blue"))
#'
#' # Without boxplot and centrality
#' do_ViolinPlot(seurat_obj, features = "CD3D",
#'               boxplot = FALSE, centrality = FALSE)
#' }
#'
#' @importFrom rlang .data %||%
#' @export
do_ViolinPlot <- function(
    object,
    features,
    group.by = NULL,
    layer = "data",
    assay = NULL,
    colors.use = NULL,
    pt.size = 3,
    pt.alpha = 0.4,
    violin.alpha = 0.2,
    violin.width = 0.5,
    boxplot = TRUE,
    boxplot.width = 0.3,
    boxplot.alpha = 0.2,
    centrality = TRUE,
    centrality.fn = "median",
    ncol = NULL
) {

  # --- Validation ---

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  if (missing(features) || length(features) == 0 || !is.character(features)) {
    stop("features must be a non-empty character vector")
  }

  # --- Resolve group.by ---

  if (is.null(group.by)) {
    # Use active Idents; inject as temporary metadata column
    object[["._badranseq_group"]] <- Seurat::Idents(object)
    group_col <- "._badranseq_group"
  } else {
    if (!group.by %in% colnames(object@meta.data)) {
      stop(
        "'", group.by, "' not found in metadata. Available columns: ",
        paste(colnames(object@meta.data), collapse = ", ")
      )
    }
    group_col <- group.by
  }

  # --- Fetch data ---

  df <- fetch_feature_data(
    object = object,
    features = features,
    layer = layer,
    assay = assay,
    metadata = group_col
  )

  # valid features are those that survived fetch_feature_data
  valid_features <- levels(df$feature)

  # --- Generate colors ---

  group_levels <- levels(as.factor(df[[group_col]]))
  n_groups <- length(group_levels)

  if (is.null(colors.use)) {
    colors.use <- generate_badranseq_colors(n_groups)
    names(colors.use) <- group_levels
  }

  # --- Build one panel per feature ---

  plots <- lapply(valid_features, function(feat) {
    feat_df <- df[df$feature == feat, ]

    .build_violin_panel(
      df = feat_df,
      y_col = layer,
      group_col = group_col,
      title = feat,
      colors = colors.use,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      violin.alpha = violin.alpha,
      violin.width = violin.width,
      boxplot = boxplot,
      boxplot.width = boxplot.width,
      boxplot.alpha = boxplot.alpha,
      centrality = centrality,
      centrality.fn = centrality.fn
    )
  })

  names(plots) <- valid_features

  # --- Return ---

  if (length(plots) == 1) {
    return(plots[[1]])
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required to combine multiple violin plots. ",
      "Install it with: install.packages('patchwork')"
    )
  }

  patchwork::wrap_plots(plots, ncol = ncol)
}
