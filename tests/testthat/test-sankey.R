# tests/testthat/test-sankey.R

data(pbmc3k, package = "BadranSeq")

# ---------------------------------------------------------------------------
# do_SankeyPlot
# ---------------------------------------------------------------------------

test_that("do_SankeyPlot returns ggplot for 2 columns", {
  skip_if_not_installed("ggalluvial")
  p <- do_SankeyPlot(
    pbmc3k,
    columns = c("seurat_annotations", "seurat_clusters")
  )
  expect_s3_class(p, "gg")
})

test_that("do_SankeyPlot returns ggplot for 3 columns", {
  skip_if_not_installed("ggalluvial")
  set.seed(1)
  pbmc3k$Phase <- sample(
    c("G1", "S", "G2M"), ncol(pbmc3k), replace = TRUE
  )
  p <- do_SankeyPlot(
    pbmc3k,
    columns = c("seurat_annotations", "seurat_clusters", "Phase")
  )
  expect_s3_class(p, "gg")
})

test_that("do_SankeyPlot errors with 1 column", {
  skip_if_not_installed("ggalluvial")
  expect_error(
    do_SankeyPlot(pbmc3k, columns = "seurat_annotations"),
    "2 or 3"
  )
})

test_that("do_SankeyPlot errors with 4 columns", {
  skip_if_not_installed("ggalluvial")
  expect_error(
    do_SankeyPlot(pbmc3k, columns = c("a", "b", "c", "d")),
    "2 or 3"
  )
})

test_that("do_SankeyPlot errors with invalid column name", {
  skip_if_not_installed("ggalluvial")
  expect_error(
    do_SankeyPlot(
      pbmc3k,
      columns = c("seurat_annotations", "nonexistent")
    ),
    "not found in metadata"
  )
})

test_that("do_SankeyPlot errors with non-Seurat object", {
  skip_if_not_installed("ggalluvial")
  expect_error(
    do_SankeyPlot(data.frame(), columns = c("a", "b")),
    "Seurat"
  )
})

test_that("do_SankeyPlot messages about dropped NAs", {
  skip_if_not_installed("ggalluvial")
  expect_message(
    do_SankeyPlot(
      pbmc3k,
      columns = c("seurat_annotations", "seurat_clusters")
    ),
    "excluded"
  )
})

test_that("do_SankeyPlot accepts custom label.size", {
  skip_if_not_installed("ggalluvial")
  p <- suppressMessages(
    do_SankeyPlot(
      pbmc3k,
      columns = c("seurat_annotations", "seurat_clusters"),
      label.size = 5
    )
  )
  expect_s3_class(p, "gg")
})
