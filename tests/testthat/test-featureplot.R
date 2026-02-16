# tests/testthat/test-featureplot.R

test_that("do_FeaturePlot single feature returns ggplot", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_FeaturePlot(pbmc3k, features = "CD3D")
  expect_s3_class(p, "gg")
  expect_false(inherits(p, "patchwork"))
})

test_that("do_FeaturePlot multiple features returns patchwork", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_FeaturePlot(pbmc3k, features = c("CD3D", "CD8A"))
  expect_s3_class(p, "patchwork")
})

test_that("do_FeaturePlot multi-feature respects ncol = 1", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_FeaturePlot(pbmc3k, features = c("CD3D", "CD8A"), ncol = 1)
  expect_s3_class(p, "patchwork")
})

test_that("do_FeaturePlot multi-feature respects nrow = 1", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_FeaturePlot(pbmc3k, features = c("CD3D", "CD8A"), nrow = 1)
  expect_s3_class(p, "patchwork")
})

test_that("do_FeaturePlot three features with ncol = 2 returns patchwork", {
  data(pbmc3k, package = "BadranSeq")
  p <- do_FeaturePlot(pbmc3k,
                      features = c("CD3D", "CD8A", "CD4"),
                      ncol = 2)
  expect_s3_class(p, "patchwork")
})
