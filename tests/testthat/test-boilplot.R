# tests/testthat/test-boilplot.R

# --- do_BoilPlot() tests ---

test_that("do_BoilPlot returns ggplot for single gene", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = "CD8A")
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("do_BoilPlot returns patchwork for multiple genes, stacked by default", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = c("GZMA", "NKG7"))
  expect_s3_class(p, "patchwork")
  expect_equal(p$patches$layout$ncol, 1)
})

test_that("do_BoilPlot passes ncol and nrow to patchwork", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = c("GZMA", "NKG7", "PRF1"), ncol = 3)
  expect_equal(p$patches$layout$ncol, 3)
  p <- do_BoilPlot(pbmc3k, features = c("GZMA", "NKG7", "PRF1"),
                   ncol = NULL, nrow = 1)
  expect_equal(p$patches$layout$nrow, 1)
})

test_that("do_BoilPlot detection rate matches fraction of non-zero cells", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = "CD8A", group.by = "seurat_annotations")

  expr <- SeuratObject::LayerData(pbmc3k, layer = "data")["CD8A", ]
  grp <- pbmc3k$seurat_annotations
  expected <- tapply(expr > 0, grp, mean)
  n_expr <- tapply(expr > 0, grp, sum)
  expected <- expected[!is.na(n_expr) & n_expr >= 10]  # default min.nonzero

  label_layer <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[1]]
  det <- label_layer$data
  observed <- stats::setNames(det$pct, as.character(det$seurat_annotations))

  expect_equal(observed[names(expected)], expected[names(expected)],
               ignore_attr = TRUE)
})

test_that("do_BoilPlot violin and points use non-zero cells only", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = "CD8A")
  expect_true(all(p$data$.value > 0))
})

test_that("do_BoilPlot drops groups below min.nonzero from the panel", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = "CD8A",
                   group.by = "seurat_annotations", min.nonzero = 50)

  expr <- SeuratObject::LayerData(pbmc3k, layer = "data")["CD8A", ]
  n_expr <- tapply(expr > 0, pbmc3k$seurat_annotations, sum)
  expected <- names(n_expr)[!is.na(n_expr) & n_expr >= 50]

  expect_setequal(levels(p$data$seurat_annotations), expected)
  label_layer <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[1]]
  expect_setequal(as.character(label_layer$data$seurat_annotations), expected)
})

test_that("do_BoilPlot skips genes with no qualifying group and errors when none", {
  data(pbmc3k, package = "BadranSeq")
  expect_warning(
    p <- do_BoilPlot(pbmc3k, features = c("NKG7", "CD8A"),
                     group.by = "seurat_annotations", min.nonzero = 200),
    "Skipping"
  )
  expect_s3_class(p, "ggplot")
  expect_error(
    suppressWarnings(do_BoilPlot(pbmc3k, features = "CD8A", min.nonzero = 5000)),
    "No gene has a group"
  )
})

test_that("do_BoilPlot group.levels subsets and orders groups", {
  data(pbmc3k, package = "BadranSeq")
  lv <- c("NK", "CD8 T", "B")
  p <- do_BoilPlot(pbmc3k, features = "NKG7",
                   group.by = "seurat_annotations", group.levels = lv)
  expect_equal(levels(p$data$seurat_annotations), lv)
  expect_true(all(as.character(p$data$seurat_annotations) %in% lv))
})

test_that("do_BoilPlot errors on invalid group.levels", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    do_BoilPlot(pbmc3k, features = "CD8A",
                group.by = "seurat_annotations", group.levels = "FAKE"),
    "group.levels not found"
  )
})

test_that("do_BoilPlot rejects metadata columns", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    do_BoilPlot(pbmc3k, features = "nCount_RNA"),
    "accepts genes only"
  )
})

test_that("do_BoilPlot warns on unknown genes and errors when none valid", {
  data(pbmc3k, package = "BadranSeq")
  expect_warning(
    do_BoilPlot(pbmc3k, features = c("CD8A", "FAKEGENE")),
    "not found"
  )
  expect_error(
    suppressWarnings(do_BoilPlot(pbmc3k, features = "FAKEGENE")),
    "No valid genes"
  )
})

test_that("do_BoilPlot validates inputs", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(do_BoilPlot("not_seurat", features = "CD8A"), "Seurat object")
  expect_error(do_BoilPlot(pbmc3k, features = character(0)), "non-empty")
  expect_error(do_BoilPlot(pbmc3k, features = "CD8A", group.by = "fake_col"),
               "not found in metadata")
  expect_error(do_BoilPlot(pbmc3k, features = "CD8A", min.nonzero = 0),
               "min.nonzero")
})

test_that("do_BoilPlot has no boxplot by default and adds one on request", {
  data(pbmc3k, package = "BadranSeq")
  geoms_of <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomBoxplot" %in% geoms_of(do_BoilPlot(pbmc3k, features = "CD8A")))
  expect_true("GeomBoxplot" %in% geoms_of(
    do_BoilPlot(pbmc3k, features = "CD8A", boxplot = TRUE)
  ))
})

test_that("do_BoilPlot respects centrality and boxplot toggles", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = "CD8A", centrality = FALSE, boxplot = FALSE)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomBoxplot" %in% geoms)
  expect_false("GeomLabel" %in% geoms)
})

test_that("do_BoilPlot builds without error", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_BoilPlot(pbmc3k, features = "CD8A", group.by = "seurat_annotations")
  expect_no_error(ggplot2::ggplot_build(p))
})
