# tests/testthat/test-dimplot.R

# --- show.legend tests ---

test_that("do_UmapPlot show.legend = FALSE hides legend", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_UmapPlot(pbmc3k, show.legend = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$theme$legend.position, "none")
})

test_that("do_UmapPlot show.legend = TRUE (default) keeps legend", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_UmapPlot(pbmc3k)
  built <- ggplot2::ggplot_build(p)
  # Default should NOT be "none"
  expect_true(is.null(built$plot$theme$legend.position) ||
              built$plot$theme$legend.position != "none")
})

test_that("do_PcaPlot show.legend = FALSE hides legend", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_PcaPlot(pbmc3k, show.legend = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$theme$legend.position, "none")
})

test_that("do_DimPlot show.legend = FALSE hides legend", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_DimPlot(pbmc3k, show.legend = FALSE)
  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$theme$legend.position, "none")
})

test_that("do_UmapPlot split.by + show.legend = FALSE hides legend", {
  data(pbmc3k, package = "BadranSeq")
  set.seed(42)
  pbmc3k$condition <- sample(c("A", "B"), ncol(pbmc3k), replace = TRUE)
  p <- do_UmapPlot(pbmc3k, split.by = "condition", show.legend = FALSE)
  expect_s3_class(p, "patchwork")
})
