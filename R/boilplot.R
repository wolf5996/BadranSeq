# R/boilplot.R
# Boil plot for zero-inflated expression, inspired by boilplot-r
# (Lara Jerman, https://github.com/laraj1/boilplot-r, MIT licence)

# Suppress R CMD check notes for dplyr/ggplot2 NSE variables
utils::globalVariables(c("n_cells", "n_expr", "pct", "pct_label",
                         "x_pos", "bar_xmin", "bar_xmax", "bar_xfill"))

# ============================================================================
# Internal: Build a single boil plot panel
# ============================================================================

#' Build a Single Boil Plot Panel
#'
#' @description
#' Internal helper that constructs a single boil plot panel for one gene:
#' a violin + boxplot + jitter + centrality layer for non-zero cells, and a
#' detection bar below zero showing the fraction of expressing cells.
#'
#' @param df data.frame. Long-format data with `.value` and the group column.
#' @param group_col character. Column name for the grouping variable (factor).
#' @param title character. Plot title (typically the gene name).
#' @param x.title character. X-axis title (the grouping variable name).
#' @param colors character. Named vector of colors for groups.
#' @param min.nonzero numeric. Minimum expressing cells for a group to be
#'   kept in the panel. Groups below it are dropped entirely.
#' @param pt.size numeric. Jitter point size.
#' @param pt.alpha numeric. Jitter point alpha.
#' @param jitter.width numeric. Horizontal jitter width.
#' @param violin.alpha numeric. Violin fill alpha.
#' @param violin.width numeric. Maximum violin width; every violin is scaled
#'   to this width (`scale = "width"`).
#' @param boxplot logical. Whether to overlay a boxplot.
#' @param boxplot.width numeric. Boxplot width.
#' @param boxplot.alpha numeric. Boxplot alpha.
#' @param centrality logical. Whether to show centrality point and label.
#' @param centrality.fn character. Centrality function name ("median" or "mean").
#'
#' @return A ggplot2 object, or NULL when no group reaches `min.nonzero`.
#' @keywords internal
.build_boil_panel <- function(
    df,
    group_col,
    title,
    x.title,
    colors,
    min.nonzero = 10,
    pt.size = 0.7,
    pt.alpha = 0.6,
    jitter.width = 0.1,
    violin.alpha = 0.8,
    violin.width = 0.9,
    boxplot = FALSE,
    boxplot.width = 0.3,
    boxplot.alpha = 0.2,
    centrality = TRUE,
    centrality.fn = "median"
) {

  group_levels <- levels(df[[group_col]])

  # --- Detection rate per group ---

  detection_df <- dplyr::summarise(
    dplyr::group_by(df, .data[[group_col]], .drop = FALSE),
    n_cells = dplyr::n(),
    n_expr = sum(.data[[".value"]] > 0, na.rm = TRUE),
    .groups = "drop"
  )
  detection_df$pct <- ifelse(
    detection_df$n_cells > 0,
    detection_df$n_expr / detection_df$n_cells,
    0
  )
  detection_df$pct_label <- paste0(round(detection_df$pct * 100), "%")

  # Drop groups with too few expressing cells from the panel entirely
  keep_groups <- as.character(
    detection_df[[group_col]][detection_df$n_expr >= min.nonzero]
  )
  if (length(keep_groups) == 0) {
    return(NULL)
  }
  group_levels <- group_levels[group_levels %in% keep_groups]
  detection_df <- detection_df[as.character(detection_df[[group_col]]) %in% group_levels, ]
  detection_df[[group_col]] <- factor(as.character(detection_df[[group_col]]),
                                      levels = group_levels)
  df <- df[as.character(df[[group_col]]) %in% group_levels, ]
  df[[group_col]] <- factor(as.character(df[[group_col]]), levels = group_levels)

  # Detection bar geometry sits below zero, scaled to the expression range
  nonzero_df <- df[df$.value > 0 & !is.na(df$.value), ]
  y_max <- if (nrow(nonzero_df) > 0) max(nonzero_df$.value) else 1

  bar_width <- 0.65
  detection_df$x_pos <- match(detection_df[[group_col]], group_levels)
  detection_df$bar_xmin <- detection_df$x_pos - bar_width / 2
  detection_df$bar_xmax <- detection_df$x_pos + bar_width / 2
  detection_df$bar_xfill <- detection_df$bar_xmin + bar_width * detection_df$pct

  bar_ymax <- -0.05 * y_max
  bar_ymin <- -0.15 * y_max
  label_y <- -0.225 * y_max

  # Base plot
  p <- ggplot2::ggplot(
    nonzero_df,
    ggplot2::aes(
      x = .data[[group_col]],
      y = .data[[".value"]],
      fill = .data[[group_col]],
      color = .data[[group_col]]
    )
  )

  # Layer 1: Violin (non-zero cells only)
  if (nrow(nonzero_df) > 0) {
    p <- p + ggplot2::geom_violin(
      alpha = violin.alpha,
      width = violin.width,
      scale = "width",
      color = "black",
      linewidth = 0.5
    )
  }

  # Layer 2: Boxplot overlay (optional)
  if (boxplot && nrow(nonzero_df) > 0) {
    p <- p + ggplot2::geom_boxplot(
      width = boxplot.width,
      alpha = boxplot.alpha,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      show.legend = FALSE
    )
  }

  # Layer 3: Jitter (non-zero raw points), filled circles with black outline
  p <- p + ggplot2::geom_jitter(
    shape = 21,
    color = "black",
    size = pt.size,
    alpha = pt.alpha,
    width = jitter.width,
    show.legend = FALSE
  )

  # Layer 4: Centrality of non-zero cells (optional)
  if (centrality && nrow(nonzero_df) > 0) {
    fn <- switch(centrality.fn,
      "median" = stats::median,
      "mean"   = base::mean,
      stats::median
    )
    centrality_df <- dplyr::summarise(
      dplyr::group_by(nonzero_df, .data[[group_col]]),
      centrality_label = round(fn(.data[[".value"]], na.rm = TRUE), 2),
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
      nudge_y = y_max * 0.05,
      inherit.aes = FALSE
    )
  }

  # Layer 5: Detection bar (grey track = all cells, fill = expressing cells).
  # Geometry and colours follow the original boilplot: no outline, flat greys.
  p <- p +
    ggplot2::geom_hline(yintercept = 0, color = "#444444", linewidth = 0.3) +
    ggplot2::geom_rect(
      data = detection_df,
      ggplot2::aes(xmin = bar_xmin, xmax = bar_xmax),
      ymin = bar_ymin,
      ymax = bar_ymax,
      fill = "#bdbdbd",
      color = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = detection_df,
      ggplot2::aes(
        xmin = bar_xmin,
        xmax = bar_xfill,
        fill = .data[[group_col]]
      ),
      ymin = bar_ymin,
      ymax = bar_ymax,
      color = NA,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = detection_df,
      ggplot2::aes(x = x_pos, label = pct_label),
      y = label_y,
      size = 3,
      color = "#333333",
      inherit.aes = FALSE
    )

  # Scales and theme
  p <- p +
    ggplot2::scale_fill_manual(values = colors, limits = group_levels) +
    ggplot2::scale_color_manual(values = colors, limits = group_levels) +
    ggplot2::scale_x_discrete(limits = group_levels, drop = FALSE) +
    ggplot2::scale_y_continuous(
      limits = c(label_y - 0.04 * y_max, NA),
      breaks = function(lims) {
        b <- pretty(c(0, lims[2]))
        b[b >= 0 & b <= lims[2]]
      }
    ) +
    ggplot2::labs(title = title, x = x.title, y = "Expression") +
    theme_badranseq() +
    ggplot2::theme(legend.position = "none")

  p
}

# ============================================================================
# Exported: do_BoilPlot
# ============================================================================

#' Boil Plot for Zero-Inflated Gene Expression
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Creates boil plots for one or more genes from a Seurat object. Single-cell
#' expression is dominated by zeros, which squash conventional violins
#' towards the baseline. A boil plot separates the two questions a violin
#' mixes together:
#'
#' * **How many cells express the gene?** A detection bar below zero shows the
#'   fraction of cells with non-zero expression (grey track = all cells,
#'   coloured fill = expressing cells, percentage label underneath).
#' * **How is the gene expressed in those cells?** The violin, jittered
#'   points, centrality marker and optional boxplot (`boxplot = TRUE`) above
#'   zero use non-zero cells only.
#'
#' Geometry and colours follow the original boil plot, drawn with BadranSeq's
#' palette and \code{\link{theme_badranseq}}. For a single gene, returns
#' a bare ggplot2 object. For multiple genes, returns a patchwork composition.
#'
#' The boil plot design comes from the boilplot package by Lara Jerman
#' (\url{https://github.com/laraj1/boilplot-r}, Python original at
#' \url{https://github.com/laraj1/boilplot-python}). This function is a
#' reimplementation in BadranSeq's visual style, with credit to the original
#' as the source of inspiration.
#'
#' @param object Seurat object.
#' @param features character. Gene names to plot (required). Metadata columns
#'   are not accepted; use \code{\link{do_ViolinPlot}} for those.
#' @param group.by character. Metadata column for grouping. If NULL (default),
#'   uses active Idents.
#' @param group.levels character. Subset of group levels to include, in the
#'   order given (default: NULL = all levels in their existing order).
#' @param layer character. Layer to extract expression from (default: "data").
#'   Use a non-negative layer such as "data" or "counts".
#' @param assay character. Assay to use (default: NULL = DefaultAssay).
#' @param colors.use character. Named vector of colors for groups. If NULL,
#'   generates colors automatically with \code{\link{generate_badranseq_colors}}.
#' @param min.nonzero numeric. Minimum number of expressing cells a group
#'   needs to appear in a gene's panel (default: 10). Groups below the
#'   threshold are removed from that panel entirely, including their
#'   detection bar and x-axis position. Genes with no qualifying group are
#'   skipped with a warning.
#' @param pt.size numeric. Jitter point size (default: 0.7).
#' @param pt.alpha numeric. Jitter point transparency (default: 0.6).
#' @param jitter.width numeric. Horizontal jitter width for points
#'   (default: 0.1).
#' @param violin.alpha numeric. Violin fill alpha (default: 0.8).
#' @param violin.width numeric. Maximum violin width (default: 0.9). All
#'   violins share this width regardless of group size, so small groups are
#'   not squashed.
#' @param boxplot logical. Whether to overlay a boxplot (default: FALSE).
#' @param boxplot.width numeric. Boxplot width (default: 0.3).
#' @param boxplot.alpha numeric. Boxplot alpha (default: 0.2).
#' @param centrality logical. Whether to show centrality point and label for
#'   non-zero cells (default: TRUE).
#' @param centrality.fn character. Centrality function ("median" or "mean",
#'   default: "median").
#' @param ncol numeric. Number of columns for patchwork layout when plotting
#'   multiple genes (default: 1, one gene per row so panels stack vertically).
#' @param nrow numeric. Number of rows for patchwork layout (default: NULL =
#'   auto). Set `ncol = NULL` together with `nrow` to control rows only.
#'
#' @return A ggplot2 object (single gene) or patchwork object (multiple genes).
#'
#' @references
#' Jerman, L. boilplot: faithful and interpretable visualisations of
#' zero-inflated distributions. \url{https://github.com/laraj1/boilplot-r}
#'
#' @seealso \code{\link{do_ViolinPlot}}, \code{\link{do_StatsViolinPlot}}
#'
#' @examples
#' \dontrun{
#' data(pbmc3k)
#'
#' # Single gene
#' do_BoilPlot(pbmc3k, features = "CD8A", group.by = "seurat_annotations")
#'
#' # Multiple genes, selected groups in a chosen order
#' do_BoilPlot(pbmc3k, features = c("GZMA", "NKG7", "PRF1"),
#'             group.by = "seurat_annotations",
#'             group.levels = c("NK", "CD8 T", "Memory CD4 T"))
#'
#' # Mean instead of median, with boxplot overlay
#' do_BoilPlot(pbmc3k, features = "MS4A1", centrality.fn = "mean",
#'             boxplot = TRUE)
#' }
#'
#' @importFrom rlang .data %||%
#' @export
do_BoilPlot <- function(
    object,
    features,
    group.by = NULL,
    group.levels = NULL,
    layer = "data",
    assay = NULL,
    colors.use = NULL,
    min.nonzero = 10,
    pt.size = 0.7,
    pt.alpha = 0.6,
    jitter.width = 0.1,
    violin.alpha = 0.8,
    violin.width = 0.9,
    boxplot = FALSE,
    boxplot.width = 0.3,
    boxplot.alpha = 0.2,
    centrality = TRUE,
    centrality.fn = "median",
    ncol = 1,
    nrow = NULL
) {

  # --- Validation ---

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  if (missing(features) || length(features) == 0 || !is.character(features)) {
    stop("features must be a non-empty character vector")
  }

  if (!is.numeric(min.nonzero) || length(min.nonzero) != 1 ||
      is.na(min.nonzero) || min.nonzero < 1) {
    stop("min.nonzero must be a single number >= 1")
  }

  features <- unique(features)

  # --- Resolve group.by ---

  if (is.null(group.by)) {
    group_col <- "._badranseq_group"
    inject_idents <- TRUE
  } else {
    if (!group.by %in% colnames(object@meta.data)) {
      stop(
        "'", group.by, "' not found in metadata. Available columns: ",
        paste(colnames(object@meta.data), collapse = ", ")
      )
    }
    group_col <- group.by
    inject_idents <- FALSE
  }

  # --- Classify features (genes only) ---

  classified <- .classify_features(object, features, assay = assay)

  if (length(classified$metadata) > 0) {
    stop(
      "do_BoilPlot() accepts genes only. Metadata columns found: ",
      paste(classified$metadata, collapse = ", "),
      ". Use do_ViolinPlot() for metadata."
    )
  }

  if (length(classified$unknown) > 0) {
    warning(
      "Features not found in assay: ",
      paste(classified$unknown, collapse = ", "), ". Skipping."
    )
  }

  valid_features <- features[features %in% classified$genes]

  if (length(valid_features) == 0) {
    stop("No valid genes found in assay.")
  }

  # --- Fetch data ---

  df <- fetch_feature_data(
    object = object,
    features = valid_features,
    layer = layer,
    assay = assay,
    metadata = if (inject_idents) FALSE else group_col
  )
  df$.value <- df[[layer]]
  df$features <- factor(as.character(df$features), levels = valid_features)

  if (any(df$.value < 0, na.rm = TRUE)) {
    warning(
      "Negative values found in layer '", layer, "'. Boil plots assume ",
      "non-negative expression (e.g. layer = 'data' or 'counts')."
    )
  }

  # Inject Idents into tibble (avoids mutating the Seurat object)
  if (inject_idents) {
    ident_lookup <- stats::setNames(
      as.character(Seurat::Idents(object)),
      colnames(object)
    )
    df[[group_col]] <- ident_lookup[df$cell_id]
  }

  # --- Filter and order group levels ---

  df <- df[!is.na(df[[group_col]]), ]

  if (!is.null(group.levels)) {
    invalid <- setdiff(group.levels, unique(as.character(df[[group_col]])))
    if (length(invalid) > 0) {
      stop(
        "group.levels not found in '", group_col, "': ",
        paste(invalid, collapse = ", ")
      )
    }
    df <- df[df[[group_col]] %in% group.levels, ]
    df[[group_col]] <- factor(as.character(df[[group_col]]),
                              levels = group.levels)
  } else {
    df[[group_col]] <- droplevels(as.factor(df[[group_col]]))
  }

  group_levels <- levels(df[[group_col]])

  # --- Generate colors ---

  if (is.null(colors.use)) {
    colors.use <- generate_badranseq_colors(length(group_levels))
    names(colors.use) <- group_levels
  }

  # --- Build one panel per gene ---

  plots <- lapply(valid_features, function(feat) {
    .build_boil_panel(
      df = df[df$features == feat, ],
      group_col = group_col,
      title = feat,
      x.title = if (inject_idents) "Identity" else group.by,
      colors = colors.use,
      min.nonzero = min.nonzero,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      jitter.width = jitter.width,
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

  # Drop genes where no group reaches min.nonzero
  empty <- vapply(plots, is.null, logical(1))
  if (any(empty)) {
    warning(
      "No group has >= ", min.nonzero, " expressing cells for: ",
      paste(valid_features[empty], collapse = ", "), ". Skipping."
    )
    plots <- plots[!empty]
  }

  if (length(plots) == 0) {
    stop("No gene has a group with >= ", min.nonzero, " expressing cells.")
  }

  # --- Return ---

  if (length(plots) == 1) {
    return(plots[[1]])
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required to combine multiple boil plots. ",
      "Install it with: install.packages('patchwork')"
    )
  }

  patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
}
