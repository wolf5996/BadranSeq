# BadranSeq - Convenience functions for single-cell RNA-seq visualization
# Native ggplot2 implementations with publication-ready aesthetics

# Suppress R CMD check notes for ggplot2 NSE variables
utils::globalVariables(c("Dim1", "Dim2", "group", "PC", "Value",
                         "dim1", "dim2", "color_var", "expression"))

# ============================================================================
# PCA Variance Helper
# ============================================================================

#' Get PCA Variance Explained Summary
#'
#' @description
#' Helper function to extract variance explained information from a Seurat object.
#'
#' @param sample Seurat object with PCA reduction.
#' @param n_pcs numeric. Number of top PCs to return. If NULL, returns all PCs.
#'
#' @return data.frame with PC number, variance explained percentage, and cumulative variance.
#'
#' @examples
#' \dontrun{
#' get_pca_variance(seurat_object, n_pcs = 10)
#' get_pca_variance(seurat_object)
#' }
#'
#' @export
get_pca_variance <- function(sample, n_pcs = NULL) {

  if (!inherits(sample, "Seurat")) {
    stop("sample must be a Seurat object")
  }

  if (!"pca" %in% names(sample@reductions)) {
    stop("PCA reduction not found in Seurat object. Please run RunPCA() first.")
  }

  pca_stdev <- sample@reductions$pca@stdev
  total_variance <- sum(pca_stdev^2)
  variance_explained <- (pca_stdev^2 / total_variance) * 100

  variance_df <- data.frame(
    PC = seq_along(variance_explained),
    variance_explained = variance_explained,
    cumulative_variance = cumsum(variance_explained)
  )

  if (!is.null(n_pcs)) {
    if (!is.numeric(n_pcs) || n_pcs < 1) {
      stop("n_pcs must be a positive integer")
    }
    if (n_pcs > nrow(variance_df)) {
      warning(paste("n_pcs exceeds available PCs. Returning all", nrow(variance_df), "PCs."))
    } else {
      variance_df <- variance_df[1:n_pcs, ]
    }
  }

  return(variance_df)
}

# ============================================================================
# Internal DimPlot Engine
# ============================================================================

#' Internal DimPlot Implementation
#'
#' @description
#' Core function that builds dimensionality reduction plots using native ggplot2.
#'
#' @param sample Seurat object.
#' @param reduction character. Name of reduction to use.
#' @param dims numeric. Vector of 2 dimensions to plot.
#' @param group.by character. Metadata column for grouping.
#' @param split.by character. Metadata column for faceting.
#' @param colors.use Named vector of colors.
#' @param label logical. Show cluster labels.
#' @param label.size numeric. Label font size.
#' @param label.color character. Label text color.
#' @param label.fill character. Label background fill.
#' @param label.box logical. Use boxed labels.
#' @param pt.size numeric. Point size.
#' @param pt.alpha numeric. Point transparency.
#' @param plot.axes logical. Show axes.
#' @param shuffle logical. Randomize cell plotting order.
#' @param raster logical. Rasterize points for large datasets.
#' @param raster.dpi numeric. DPI for rasterization.
#' @param na.value character. Color for NA values.
#'
#' @return ggplot2 object.
#' @keywords internal
.do_DimPlot_internal <- function(
    sample,
    reduction = "umap",
    dims = c(1, 2),
    group.by = NULL,
    split.by = NULL,
    colors.use = NULL,
    label = TRUE,
    label.size = 4,
    label.color = "black",
    label.fill = "white",
    label.box = TRUE,
    pt.size = 1,
    pt.alpha = 1,
    plot.axes = TRUE,
    shuffle = TRUE,
    raster = NULL,
    raster.dpi = 300,
    na.value = "grey75"
) {

  # ---- Data Extraction ----
  embeddings <- sample@reductions[[reduction]]@cell.embeddings[, dims]
  colnames(embeddings) <- c("Dim1", "Dim2")

  plot_data <- as.data.frame(embeddings)
  plot_data$cell_id <- rownames(plot_data)

  # Determine grouping variable
  if (is.null(group.by)) {
    plot_data$group <- as.factor(Seurat::Idents(sample))
    group_var <- "Identity"
  } else {
    if (!group.by %in% colnames(sample@meta.data)) {
      stop(paste0("'", group.by, "' not found in metadata. Available columns: ",
                  paste(colnames(sample@meta.data), collapse = ", ")))
    }
    plot_data$group <- as.factor(sample@meta.data[[group.by]])
    group_var <- group.by
  }

  # Add split.by variable if specified
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(sample@meta.data)) {
      stop(paste0("'", split.by, "' not found in metadata."))
    }
    plot_data$split_var <- as.factor(sample@meta.data[[split.by]])
  }

  # ---- Shuffle cells ----
  if (shuffle) {
    plot_data <- plot_data[sample(nrow(plot_data)), ]
  }

  # ---- Color setup ----
  unique_groups <- levels(plot_data$group)
  if (is.null(unique_groups)) {
    unique_groups <- unique(as.character(plot_data$group))
  }
  n_groups <- length(unique_groups)

  if (is.null(colors.use)) {
    colors.use <- generate_badranseq_colors(n_groups)
    names(colors.use) <- unique_groups
  } else {
    # Validate provided colors
    if (!all(unique_groups %in% names(colors.use))) {
      missing <- setdiff(unique_groups, names(colors.use))
      warning(paste("Missing colors for groups:", paste(missing, collapse = ", "),
                    ". Using default colors for missing groups."))
      default_colors <- generate_badranseq_colors(length(missing))
      names(default_colors) <- missing
      colors.use <- c(colors.use, default_colors)
    }
  }

  # ---- Auto rasterization for large datasets ----
  if (is.null(raster)) {
    raster <- nrow(plot_data) > 50000
  }

  # ---- Build base plot ----
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dim1, y = Dim2, color = group))

  # Add points (with optional rasterization)
  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    p <- p + ggrastr::rasterise(
      ggplot2::geom_point(size = pt.size, alpha = pt.alpha),
      dpi = raster.dpi
    )
  } else {
    p <- p + ggplot2::geom_point(size = pt.size, alpha = pt.alpha)
  }

  # ---- Add labels ----
  if (label) {
    # Calculate label positions (group centroids)
    label_data <- stats::aggregate(
      cbind(Dim1, Dim2) ~ group,
      data = plot_data,
      FUN = mean
    )

    if (label.box) {
      p <- p + ggplot2::geom_label(
        data = label_data,
        ggplot2::aes(x = Dim1, y = Dim2, label = group),
        color = label.color,
        fill = label.fill,
        size = label.size,
        fontface = "bold",
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(x = Dim1, y = Dim2, label = group),
        color = label.color,
        size = label.size,
        fontface = "bold",
        show.legend = FALSE
      )
    }
  }

  # ---- Apply colors and theme ----
  p <- p +
    ggplot2::scale_color_manual(values = colors.use, name = group_var, na.value = na.value) +
    theme_badranseq(plot.axes = plot.axes) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1)
      )
    )

  # ---- Faceting ----
  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~ split_var, scales = "fixed")
  }

  return(p)
}

# ============================================================================
# Dimensionality Reduction Plotting - Public Functions
# ============================================================================

#' Enhanced PCA Plot with Variance Explained
#'
#' @description
#' Creates a PCA dimensionality reduction plot with variance explained
#' percentages automatically added to axis labels.
#'
#' @param sample Seurat object.
#' @param dims numeric. Vector of 2 PC dimensions to plot (default: c(1, 2)).
#' @param group.by character. Metadata column for grouping (default: active identity).
#' @param split.by character. Metadata column for faceting.
#' @param colors.use Named vector of colors for groups.
#' @param label logical. Show cluster labels (default: TRUE).
#' @param label.size numeric. Label font size (default: 4).
#' @param label.color character. Label text color (default: "black").
#' @param label.fill character. Label background fill (default: "white").
#' @param pt.size numeric. Point size (default: 1).
#' @param pt.alpha numeric. Point transparency (default: 1).
#' @param plot.axes logical. Show axes (default: TRUE).
#' @param variance_digits numeric. Decimal places for variance (default: 1).
#' @param shuffle logical. Randomize cell plotting order (default: TRUE).
#' @param ... Additional arguments (currently unused).
#'
#' @return ggplot2 object with PCA plot including variance explained in axis labels.
#'
#' @examples
#' \dontrun{
#' do_PcaPlot(seurat_object, group.by = "condition")
#' do_PcaPlot(seurat_object, dims = c(2, 3), group.by = "cluster")
#' }
#'
#' @export
do_PcaPlot <- function(
    sample,
    dims = c(1, 2),
    group.by = NULL,
    split.by = NULL,
    colors.use = NULL,
    label = TRUE,
    label.size = 4,
    label.color = "black",
    label.fill = "white",
    pt.size = 1,
    pt.alpha = 1,
    plot.axes = TRUE,
    variance_digits = 1,
    shuffle = TRUE,
    ...
) {

  # Validation
  if (!inherits(sample, "Seurat")) {
    stop("sample must be a Seurat object")
  }

  if (!"pca" %in% names(sample@reductions)) {
    stop("PCA reduction not found. Please run RunPCA() first.")
  }

  if (length(dims) != 2 || !is.numeric(dims) || any(dims < 1)) {
    stop("dims must be a vector of 2 positive integers")
  }

  max_pc <- ncol(sample@reductions$pca@cell.embeddings)
  if (any(dims > max_pc)) {
    stop(paste("Requested PC dimensions exceed available PCs. Maximum PC:", max_pc))
  }

  # Calculate variance explained
  pca_stdev <- sample@reductions$pca@stdev
  total_variance <- sum(pca_stdev^2)
  variance_explained <- (pca_stdev^2 / total_variance) * 100

  pc1_var <- variance_explained[dims[1]]
  pc2_var <- variance_explained[dims[2]]

  pc1_label <- paste0("PC", dims[1], " (", round(pc1_var, variance_digits), "%)")
  pc2_label <- paste0("PC", dims[2], " (", round(pc2_var, variance_digits), "%)")

  # Build plot
  plot <- .do_DimPlot_internal(
    sample = sample,
    reduction = "pca",
    dims = dims,
    group.by = group.by,
    split.by = split.by,
    colors.use = colors.use,
    label = label,
    label.size = label.size,
    label.color = label.color,
    label.fill = label.fill,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    plot.axes = plot.axes,
    shuffle = shuffle
  )

  # Add variance labels
  plot <- plot + ggplot2::labs(x = pc1_label, y = pc2_label)

  return(plot)
}

#' Enhanced UMAP Plot
#'
#' @description
#' Creates a UMAP dimensionality reduction plot with publication-ready styling.
#'
#' @inheritParams do_PcaPlot
#'
#' @return ggplot2 object with UMAP plot.
#'
#' @examples
#' \dontrun{
#' do_UmapPlot(seurat_object, group.by = "condition")
#' do_UmapPlot(seurat_object, group.by = "condition", label = FALSE)
#' }
#'
#' @export
do_UmapPlot <- function(
    sample,
    dims = c(1, 2),
    group.by = NULL,
    split.by = NULL,
    colors.use = NULL,
    label = TRUE,
    label.size = 4,
    label.color = "black",
    label.fill = "white",
    pt.size = 1,
    pt.alpha = 1,
    plot.axes = TRUE,
    shuffle = TRUE,
    ...
) {

  # Validation
  if (!inherits(sample, "Seurat")) {
    stop("sample must be a Seurat object")
  }

  if (!"umap" %in% names(sample@reductions)) {
    stop("UMAP reduction not found. Please run RunUMAP() first.")
  }

  if (length(dims) != 2 || !is.numeric(dims) || any(dims < 1)) {
    stop("dims must be a vector of 2 positive integers")
  }

  max_dim <- ncol(sample@reductions$umap@cell.embeddings)
  if (any(dims > max_dim)) {
    stop(paste("Requested UMAP dimensions exceed available dimensions. Maximum:", max_dim))
  }

  # Build plot
  plot <- .do_DimPlot_internal(
    sample = sample,
    reduction = "umap",
    dims = dims,
    group.by = group.by,
    split.by = split.by,
    colors.use = colors.use,
    label = label,
    label.size = label.size,
    label.color = label.color,
    label.fill = label.fill,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    plot.axes = plot.axes,
    shuffle = shuffle
  )

  # Standard UMAP axis labels
  plot <- plot + ggplot2::labs(
    x = paste0("UMAP", dims[1]),
    y = paste0("UMAP", dims[2])
  )

  return(plot)
}

#' Generalized Dimensionality Reduction Plot
#'
#' @description
#' Unified wrapper for dimensionality reduction plots. Routes to specialized
#' functions for PCA and UMAP, handles other reductions directly.
#'
#' @inheritParams do_PcaPlot
#' @param reduction character. Reduction to use (default: "umap").
#'
#' @return ggplot2 object with dimensionality reduction plot.
#'
#' @examples
#' \dontrun{
#' do_DimPlot(seurat_object, group.by = "condition")
#' do_DimPlot(seurat_object, reduction = "pca", group.by = "cluster")
#' do_DimPlot(seurat_object, reduction = "harmony", group.by = "batch")
#' }
#'
#' @export
do_DimPlot <- function(
    sample,
    reduction = "umap",
    dims = c(1, 2),
    group.by = NULL,
    split.by = NULL,
    colors.use = NULL,
    label = TRUE,
    label.size = 4,
    label.color = "black",
    label.fill = "white",
    pt.size = 1,
    pt.alpha = 1,
    plot.axes = TRUE,
    variance_digits = 1,
    shuffle = TRUE,
    ...
) {

  # Validation
  if (!inherits(sample, "Seurat")) {
    stop("sample must be a Seurat object")
  }

  if (!reduction %in% names(sample@reductions)) {
    stop(paste0("'", reduction, "' not found. Available reductions: ",
                paste(names(sample@reductions), collapse = ", ")))
  }

  # Route to specialized functions
  if (reduction == "pca") {
    return(do_PcaPlot(
      sample = sample,
      dims = dims,
      group.by = group.by,
      split.by = split.by,
      colors.use = colors.use,
      label = label,
      label.size = label.size,
      label.color = label.color,
      label.fill = label.fill,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      plot.axes = plot.axes,
      variance_digits = variance_digits,
      shuffle = shuffle,
      ...
    ))
  } else if (reduction == "umap") {
    return(do_UmapPlot(
      sample = sample,
      dims = dims,
      group.by = group.by,
      split.by = split.by,
      colors.use = colors.use,
      label = label,
      label.size = label.size,
      label.color = label.color,
      label.fill = label.fill,
      pt.size = pt.size,
      pt.alpha = pt.alpha,
      plot.axes = plot.axes,
      shuffle = shuffle,
      ...
    ))
  }

  # Handle other reductions (harmony, tsne, etc.)
  plot <- .do_DimPlot_internal(
    sample = sample,
    reduction = reduction,
    dims = dims,
    group.by = group.by,
    split.by = split.by,
    colors.use = colors.use,
    label = label,
    label.size = label.size,
    label.color = label.color,
    label.fill = label.fill,
    pt.size = pt.size,
    pt.alpha = pt.alpha,
    plot.axes = plot.axes,
    shuffle = shuffle
  )

  # Standardized axis labels
  plot <- plot + ggplot2::labs(
    x = paste0(toupper(reduction), dims[1]),
    y = paste0(toupper(reduction), dims[2])
  )

  return(plot)
}

# ============================================================================
# Feature Plot
# ============================================================================

#' Enhanced Feature Plot with Viridis Support
#'
#' @description
#' Creates feature expression plots overlaid on dimensionality reductions
#' with viridis color scales and publication-ready styling.
#'
#' @inheritParams do_DimPlot
#' @param features character. Gene names or metadata columns to plot.
#' @param slot character. Data slot ("data", "counts", "scale.data"). Default: "data".
#' @param use_viridis logical. Use viridis color scale (default: TRUE).
#' @param viridis.palette character. Viridis palette "A"-"H" (default: "B").
#' @param viridis.direction numeric. Color direction -1 or 1 (default: -1).
#' @param order logical. Plot high-expressing cells on top (default: TRUE).
#' @param min.cutoff numeric or "q10". Minimum expression cutoff.
#' @param max.cutoff numeric or "q90". Maximum expression cutoff.
#' @param na.value character. Color for NA/zero values (default: "grey90").
#'
#' @return ggplot2 object or list of ggplot2 objects for multiple features.
#'
#' @examples
#' \dontrun{
#' do_FeaturePlot(seurat_object, features = "CD3D")
#' do_FeaturePlot(seurat_object, features = c("CD3D", "CD8A", "CD4"))
#' do_FeaturePlot(seurat_object, features = "CD3D", reduction = "pca")
#' }
#'
#' @export
do_FeaturePlot <- function(
    sample,
    features,
    reduction = "umap",
    dims = c(1, 2),
    slot = "data",
    use_viridis = TRUE,
    viridis.palette = "B",
    viridis.direction = -1,
    order = TRUE,
    pt.size = 1,
    pt.alpha = 1,
    plot.axes = TRUE,
    variance_digits = 1,
    min.cutoff = NA,
    max.cutoff = NA,
    na.value = "grey90",
    split.by = NULL,
    ...
) {

  # ---- Validation ----
  if (!inherits(sample, "Seurat")) {
    stop("sample must be a Seurat object")
  }

  if (!reduction %in% names(sample@reductions)) {
    stop(paste0("'", reduction, "' not found. Available reductions: ",
                paste(names(sample@reductions), collapse = ", ")))
  }

  if (missing(features) || length(features) == 0) {
    stop("features must be specified")
  }

  # ---- Get embeddings ----
  embeddings <- sample@reductions[[reduction]]@cell.embeddings[, dims]
  colnames(embeddings) <- c("Dim1", "Dim2")
  plot_data <- as.data.frame(embeddings)
  plot_data$cell_id <- rownames(plot_data)

  # ---- Add split.by if specified ----
  if (!is.null(split.by)) {
    if (!split.by %in% colnames(sample@meta.data)) {
      stop(paste0("'", split.by, "' not found in metadata."))
    }
    plot_data$split_var <- as.factor(sample@meta.data[[split.by]])
  }

  # ---- Fetch expression data ----
  expr_data <- Seurat::FetchData(sample, vars = features, slot = slot)

  # ---- Generate axis labels ----
  if (reduction == "pca") {
    pca_stdev <- sample@reductions$pca@stdev
    total_variance <- sum(pca_stdev^2)
    variance_explained <- (pca_stdev^2 / total_variance) * 100
    x_label <- paste0("PC", dims[1], " (", round(variance_explained[dims[1]], variance_digits), "%)")
    y_label <- paste0("PC", dims[2], " (", round(variance_explained[dims[2]], variance_digits), "%)")
  } else {
    x_label <- paste0(toupper(reduction), dims[1])
    y_label <- paste0(toupper(reduction), dims[2])
  }

  # ---- Create plot for each feature ----
  plot_list <- lapply(features, function(feat) {

    # Check if feature exists
    if (!feat %in% colnames(expr_data)) {
      warning(paste0("Feature '", feat, "' not found. Skipping."))
      return(NULL)
    }

    # Add expression to plot data
    feature_data <- plot_data
    feature_data$expression <- expr_data[[feat]]

    # Apply cutoffs
    if (!is.na(min.cutoff)) {
      if (is.character(min.cutoff) && min.cutoff == "q10") {
        min_val <- stats::quantile(feature_data$expression, 0.10, na.rm = TRUE)
      } else {
        min_val <- as.numeric(min.cutoff)
      }
      feature_data$expression[feature_data$expression < min_val] <- min_val
    }

    if (!is.na(max.cutoff)) {
      if (is.character(max.cutoff) && max.cutoff == "q90") {
        max_val <- stats::quantile(feature_data$expression, 0.90, na.rm = TRUE)
      } else {
        max_val <- as.numeric(max.cutoff)
      }
      feature_data$expression[feature_data$expression > max_val] <- max_val
    }

    # Order by expression (plot high values on top)
    if (order) {
      feature_data <- feature_data[order(feature_data$expression, na.last = FALSE), ]
    }

    # Build plot
    p <- ggplot2::ggplot(feature_data, ggplot2::aes(x = Dim1, y = Dim2, color = expression)) +
      ggplot2::geom_point(size = pt.size, alpha = pt.alpha) +
      ggplot2::labs(
        x = x_label,
        y = y_label,
        title = feat,
        color = "Expression"
      ) +
      theme_badranseq(plot.axes = plot.axes)

    # Color scale
    if (use_viridis) {
      p <- p + ggplot2::scale_color_viridis_c(
        option = viridis.palette,
        direction = viridis.direction,
        na.value = na.value
      )
    } else {
      p <- p + ggplot2::scale_color_gradient(
        low = "grey90",
        high = "darkblue",
        na.value = na.value
      )
    }

    # Faceting
    if (!is.null(split.by)) {
      p <- p + ggplot2::facet_wrap(~ split_var, scales = "fixed")
    }

    return(p)
  })

  # Remove NULL entries
  plot_list <- Filter(Negate(is.null), plot_list)

  if (length(plot_list) == 0) {
    stop("No valid features found.")
  }

  # Return single plot or list
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    names(plot_list) <- features[features %in% colnames(expr_data)]
    return(plot_list)
  }
}

# ============================================================================
# Elbow Plot
# ============================================================================

#' Enhanced Elbow Plot for PCA Dimensionality
#'
#' @description
#' An improved version of Seurat's ElbowPlot with enhanced aesthetics,
#' customizable titles, optional cutoff lines, and variance explained options.
#'
#' @param object Seurat object.
#' @param ndims numeric. Number of dimensions to plot (default: 50).
#' @param reduction character. Reduction to use (default: "pca").
#' @param cutoff_pc numeric. PC number to mark with vertical line.
#' @param variance logical. Plot variance explained instead of stdev (default: TRUE).
#' @param title character. Custom title.
#' @param point_size numeric. Point size (default: 3).
#' @param point_color character. Point color (default: "steelblue").
#' @param line_color character. Line color (default: "steelblue").
#' @param cutoff_color character. Cutoff line color (default: "red").
#' @param base_size numeric. Base font size (default: 14).
#'
#' @return ggplot2 object.
#'
#' @examples
#' \dontrun{
#' EnhancedElbowPlot(seurat_obj)
#' EnhancedElbowPlot(seurat_obj, cutoff_pc = 10)
#' EnhancedElbowPlot(seurat_obj, variance = FALSE)
#' }
#'
#' @export
EnhancedElbowPlot <- function(
    object,
    ndims = 50,
    reduction = "pca",
    cutoff_pc = NULL,
    variance = TRUE,
    title = NULL,
    point_size = 3,
    point_color = "steelblue",
    line_color = "steelblue",
    cutoff_color = "red",
    base_size = 14
) {

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  stdev <- Seurat::Stdev(object = object, reduction = reduction)

  if (length(stdev) == 0) {
    stop(paste("No standard deviation info stored for", reduction))
  }

  if (ndims > length(stdev)) {
    warning("The object only has information for ", length(stdev), " dimensions. ",
            "Using ndims = ", length(stdev), " instead.")
    ndims <- length(stdev)
  }

  stdev <- stdev[1:ndims]
  dims <- 1:ndims

  if (variance) {
    y_values <- (stdev^2 / sum(stdev^2)) * 100
    y_label <- "Variance Explained (%)"
    if (is.null(title)) {
      title <- "PCA Variance Explained by Each Principal Component"
    }
  } else {
    y_values <- stdev
    y_label <- "Standard Deviation"
    if (is.null(title)) {
      title <- "PCA Elbow Plot: Standard Deviation by Principal Component"
    }
  }

  plot_data <- data.frame(PC = dims, Value = y_values)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = PC, y = Value)) +
    ggplot2::geom_line(color = line_color, linewidth = 0.8) +
    ggplot2::geom_point(color = point_color, size = point_size) +
    ggplot2::labs(title = title, x = "Principal Component", y = y_label) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  if (!is.null(cutoff_pc)) {
    if (cutoff_pc > 0 && cutoff_pc <= ndims) {
      p <- p +
        ggplot2::geom_vline(xintercept = cutoff_pc, linetype = "dashed",
                            color = cutoff_color, linewidth = 1) +
        ggplot2::annotate("text", x = cutoff_pc, y = max(y_values) * 0.95,
                          label = paste0("PC ", cutoff_pc), color = cutoff_color,
                          fontface = "bold",
                          hjust = ifelse(cutoff_pc < ndims * 0.75, -0.1, 1.1))
    } else {
      warning("cutoff_pc is outside the range of plotted dimensions. Ignoring.")
    }
  }

  return(p)
}

# ============================================================================
# Sleepwalk Integration
# ============================================================================

#' SleepWalk Wrapper for Seurat Objects
#'
#' @param object Seurat object.
#' @param embedding character. Reduction for 2D embedding (default: "umap").
#' @param features character. Reduction for feature space (default: "pca").
#'
#' @return sleepwalk interactive plot.
#'
#' @examples
#' \dontrun{
#' seurat_sleepwalk(seurat_obj, embedding = "umap", features = "pca")
#' }
#'
#' @export
seurat_sleepwalk <- function(object, embedding = "umap", features = "pca") {
  sleepwalk::sleepwalk(
    embeddings = Seurat::Embeddings(object, reduction = embedding),
    featureMatrices = Seurat::Embeddings(object, reduction = features)
  )
}

# ============================================================================
# Interactive Cell Selection
# ============================================================================

#' Interactive Cell Selection for Seurat Objects
#'
#' @param seurat_obj Seurat object with dimensionality reduction computed.
#' @param reduction character. Reduction to use ("umap", "tsne", "pca").
#' @param color_by character. Metadata column for coloring (default: "Idents").
#' @param return_cells logical. Return cell names instead of Seurat subset.
#'
#' @return Subsetted Seurat object or character vector of cell names.
#'
#' @examples
#' \dontrun{
#' seu_subset <- select_cells_interactive(seurat_obj)
#' cell_names <- select_cells_interactive(seurat_obj, return_cells = TRUE)
#' }
#'
#' @export
select_cells_interactive <- function(
    seurat_obj,
    reduction = c("umap", "tsne", "pca"),
    color_by = "Idents",
    return_cells = FALSE
) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required for interactive selection")
  }

  if (!interactive()) {
    stop("This function only works in interactive R sessions")
  }

  reduction <- match.arg(reduction)

  if (is.null(seurat_obj@reductions[[reduction]])) {
    stop("Reduction '", reduction, "' not found. Please compute it first.")
  }

  if (color_by != "Idents" && !color_by %in% colnames(seurat_obj@meta.data)) {
    stop("color_by '", color_by, "' not found in metadata. Available columns: ",
         paste(colnames(seurat_obj@meta.data), collapse = ", "))
  }

  coords <- as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings[, 1:2])
  colnames(coords) <- c("dim1", "dim2")
  coords$cell_name <- rownames(coords)

  if (color_by == "Idents") {
    coords$color_var <- as.character(Seurat::Idents(seurat_obj))
    color_label <- "Identity"
  } else {
    coords$color_var <- as.character(seurat_obj@meta.data[[color_by]])
    color_label <- color_by
  }

  ui <- shiny::fluidPage(
    shiny::titlePanel("Interactive Cell Selection"),
    shiny::fluidRow(
      shiny::column(3,
        shiny::wellPanel(
          shiny::h4("Controls"),
          shiny::actionButton("select_btn", "Select Brushed Cells",
                              class = "btn-primary", width = "100%"),
          shiny::br(), shiny::br(),
          shiny::actionButton("deselect_btn", "Deselect Brushed Cells",
                              class = "btn-warning", width = "100%"),
          shiny::br(), shiny::br(),
          shiny::actionButton("clear_btn", "Clear All",
                              class = "btn-danger", width = "100%"),
          shiny::br(), shiny::br(),
          shiny::actionButton("done_btn", "Done",
                              class = "btn-success", width = "100%"),
          shiny::hr(),
          shiny::h5("Plot Settings:"),
          shiny::p(paste("Reduction:", toupper(reduction))),
          shiny::p(paste("Colored by:", color_label)),
          shiny::hr(),
          shiny::h5("Instructions:"),
          shiny::p("1. Click and drag to brush cells"),
          shiny::p("2. Click 'Select' or 'Deselect'"),
          shiny::p("3. Repeat as needed"),
          shiny::p("4. Click 'Done' when finished"),
          shiny::hr(),
          shiny::textOutput("selection_count")
        )
      ),
      shiny::column(9,
        shiny::plotOutput("main_plot", height = "600px",
                          brush = shiny::brushOpts(id = "plot_brush"))
      )
    )
  )

  server <- function(input, output, session) {
    values <- shiny::reactiveValues(selected = rep(FALSE, nrow(coords)))

    output$main_plot <- shiny::renderPlot({
      plot_data <- coords
      plot_data$selected <- values$selected

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = dim1, y = dim2)) +
        ggplot2::geom_point(ggplot2::aes(color = color_var), size = 1.5, alpha = 0.7) +
        ggplot2::labs(x = paste0(toupper(reduction), "_1"),
                      y = paste0(toupper(reduction), "_2"),
                      color = color_label) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "right",
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::ggtitle("Select Cells by Brushing")

      if (any(values$selected)) {
        selected_data <- plot_data[plot_data$selected, ]
        p <- p + ggplot2::geom_point(
          data = selected_data,
          ggplot2::aes(x = dim1, y = dim2),
          color = "red", size = 2.5, shape = 1, stroke = 1.5
        )
      }
      p
    })

    output$selection_count <- shiny::renderText({
      paste("Selected:", sum(values$selected), "cells")
    })

    shiny::observeEvent(input$select_btn, {
      if (!is.null(input$plot_brush)) {
        brushed <- shiny::brushedPoints(coords, input$plot_brush,
                                        xvar = "dim1", yvar = "dim2")
        if (nrow(brushed) > 0) {
          indices <- match(brushed$cell_name, coords$cell_name)
          values$selected[indices] <- TRUE
        }
      }
    })

    shiny::observeEvent(input$deselect_btn, {
      if (!is.null(input$plot_brush)) {
        brushed <- shiny::brushedPoints(coords, input$plot_brush,
                                        xvar = "dim1", yvar = "dim2")
        if (nrow(brushed) > 0) {
          indices <- match(brushed$cell_name, coords$cell_name)
          values$selected[indices] <- FALSE
        }
      }
    })

    shiny::observeEvent(input$clear_btn, {
      values$selected <- rep(FALSE, nrow(coords))
    })

    shiny::observeEvent(input$done_btn, {
      shiny::stopApp(values$selected)
    })
  }

  message("Starting interactive cell selector...")
  selected_logical <- shiny::runApp(shiny::shinyApp(ui, server))

  if (is.null(selected_logical) || !any(selected_logical)) {
    message("No cells selected")
    return(NULL)
  }

  selected_cells <- coords$cell_name[selected_logical]
  message("Selected ", length(selected_cells), " cells")

  if (return_cells) {
    return(selected_cells)
  } else {
    return(subset(seurat_obj, cells = selected_cells))
  }
}
