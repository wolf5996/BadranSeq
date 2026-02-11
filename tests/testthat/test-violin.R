# tests/testthat/test-violin.R

# --- do_ViolinPlot() tests ---

test_that("do_ViolinPlot returns ggplot for single gene", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_ViolinPlot(pbmc3k, features = "CD3D")
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("do_ViolinPlot returns patchwork for multiple genes", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_ViolinPlot(pbmc3k, features = c("CD3D", "CD8A"))
  expect_s3_class(p, "patchwork")
})

test_that("do_ViolinPlot uses Idents when group.by is NULL", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_ViolinPlot(pbmc3k, features = "CD3D", group.by = NULL)
  # Should work without error and produce a ggplot

  expect_s3_class(p, "gg")
  # Build the plot to access its data
  built <- ggplot2::ggplot_build(p)
  # The number of groups should match the number of Idents levels
  n_idents <- length(levels(Seurat::Idents(pbmc3k)))
  # Check that the fill scale has the right number of colors
  expect_equal(length(built$plot$scales$scales[[1]]$palette.cache), n_idents)
})

test_that("do_ViolinPlot works with single group level", {

  data(pbmc3k, package = "BadranSeq")
  # Subset to one cluster
  sub <- subset(pbmc3k, idents = "0")
  p <- do_ViolinPlot(sub, features = "CD3D")
  expect_s3_class(p, "gg")
})

test_that("do_ViolinPlot respects colors.use", {
  data(pbmc3k, package = "BadranSeq")
  custom_colors <- c("0" = "red", "1" = "blue", "2" = "green",
                     "3" = "purple", "4" = "orange", "5" = "cyan",
                     "6" = "pink", "7" = "yellow", "8" = "brown",
                     "9" = "grey")
  p <- do_ViolinPlot(pbmc3k, features = "CD3D", colors.use = custom_colors)
  built <- ggplot2::ggplot_build(p)
  # The fill scale values should include our custom colors
  scale_fills <- p$scales$scales[[1]]$palette.cache
  expect_true("red" %in% scale_fills || "red" %in% built$plot$scales$scales[[1]]$palette.cache)
})

test_that("do_ViolinPlot errors on invalid features (empty character vector)", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    do_ViolinPlot(pbmc3k, features = character(0)),
    "features must be a non-empty character vector"
  )
})

test_that("do_ViolinPlot errors on non-Seurat object", {
  expect_error(
    do_ViolinPlot("not_a_seurat", features = "CD3D"),
    "object must be a Seurat object"
  )
})

test_that("do_ViolinPlot errors on invalid group.by", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    do_ViolinPlot(pbmc3k, features = "CD3D", group.by = "nonexistent_column"),
    "not found in metadata"
  )
})

test_that("do_ViolinPlot warns on missing features (some valid, some not)", {
  data(pbmc3k, package = "BadranSeq")
  expect_warning(
    p <- do_ViolinPlot(pbmc3k, features = c("CD3D", "FAKEGENE123")),
    "not found"
  )
  # Should still return a plot for the valid feature
  expect_s3_class(p, "gg")
})

test_that("do_ViolinPlot can disable boxplot and centrality", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_ViolinPlot(pbmc3k, features = "CD3D",
                     boxplot = FALSE, centrality = FALSE)
  expect_s3_class(p, "gg")
  # Build the plot and check layers - should have fewer layers
  built <- ggplot2::ggplot_build(p)
  # With boxplot=FALSE and centrality=FALSE, we expect:
  # geom_violin + geom_jitter = 2 layers (no boxplot, no centrality diamond/label)
  # Count the layer types
  layer_types <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomBoxplot" %in% layer_types)
  # No diamond point (centrality) or label layers
  expect_false("GeomLabel" %in% layer_types)
})

test_that("do_ViolinPlot ncol controls patchwork layout", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_ViolinPlot(pbmc3k, features = c("CD3D", "CD8A"), ncol = 1)
  expect_s3_class(p, "patchwork")
  # The patchwork layout should have ncol = 1
  # patchwork stores layout info in the plot
  expect_equal(p$patches$layout$ncol, 1)
})

# --- do_StatsViolinPlot() tests ---

test_that("do_StatsViolinPlot returns ggplot for single gene", {
  skip_if_not_installed("ggstatsplot")
  data(pbmc3k, package = "BadranSeq")
  p <- do_StatsViolinPlot(pbmc3k, features = "CD3D",
                           group.by = "seurat_clusters")
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("do_StatsViolinPlot returns patchwork for multiple genes", {
  skip_if_not_installed("ggstatsplot")
  data(pbmc3k, package = "BadranSeq")
  p <- do_StatsViolinPlot(pbmc3k, features = c("CD3D", "CD8A"),
                           group.by = "seurat_clusters")
  expect_s3_class(p, "patchwork")
})

test_that("do_StatsViolinPlot uses Idents when group.by is NULL", {
  skip_if_not_installed("ggstatsplot")
  data(pbmc3k, package = "BadranSeq")
  Seurat::Idents(pbmc3k) <- "seurat_clusters"
  p <- do_StatsViolinPlot(pbmc3k, features = "CD3D", group.by = NULL)
  expect_s3_class(p, "gg")
})

test_that("do_StatsViolinPlot errors on single group level", {
  skip_if_not_installed("ggstatsplot")
  data(pbmc3k, package = "BadranSeq")
  sub <- subset(pbmc3k, idents = "0")
  expect_error(
    do_StatsViolinPlot(sub, features = "CD3D"),
    "group.by has only 1 level"
  )
})

test_that("do_StatsViolinPlot errors on non-Seurat object", {
  expect_error(
    do_StatsViolinPlot("not_a_seurat", features = "CD3D"),
    "object must be a Seurat object"
  )
})

test_that("do_StatsViolinPlot errors on invalid group.by", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    do_StatsViolinPlot(pbmc3k, features = "CD3D",
                        group.by = "nonexistent_column"),
    "not found in metadata"
  )
})

test_that("do_StatsViolinPlot errors without ggstatsplot", {
  skip_if(requireNamespace("ggstatsplot", quietly = TRUE),
          "ggstatsplot is installed, skipping missing-package test")
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    do_StatsViolinPlot(pbmc3k, features = "CD3D",
                        group.by = "seurat_clusters"),
    "ggstatsplot"
  )
})

test_that("do_StatsViolinPlot passes type argument through", {
  skip_if_not_installed("ggstatsplot")
  data(pbmc3k, package = "BadranSeq")
  p <- do_StatsViolinPlot(pbmc3k, features = "CD3D",
                           group.by = "seurat_clusters",
                           type = "parametric")
  expect_s3_class(p, "gg")
})
