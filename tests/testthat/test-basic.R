test_that("Package loads correctly", {
  expect_true(TRUE)
})

test_that("pbmc3k data is available", {
  data(pbmc3k, package = "BadranSeq")
  expect_true(exists("pbmc3k"))
  expect_s4_class(pbmc3k, "Seurat")
})
