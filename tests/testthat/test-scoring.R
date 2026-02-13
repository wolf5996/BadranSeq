# tests/testthat/test-scoring.R

test_that("score_signatures rejects non-Seurat input", {
  expect_error(score_signatures("x"), "must be a Seurat object")
})

test_that("score_signatures errors when msigdbr missing", {
  skip_if(requireNamespace("msigdbr", quietly = TRUE),
          "msigdbr is installed, skipping missing-package test")
  expect_error(score_signatures(BadranSeq::pbmc3k), "msigdbr")
})

test_that("score_signatures errors on empty collection", {
  skip_if_not_installed("msigdbr")
  expect_error(
    score_signatures(BadranSeq::pbmc3k, collection = "C2",
                     subcollection = "FAKE_SUB"),
    "No gene sets found"
  )
})

test_that("score_signatures detects Ensembl IDs", {
  skip_if_not_installed("msigdbr")
  mat <- matrix(1:4, nrow = 2,
                dimnames = list(c("ENSG00000141510", "ENSG00000012048"),
                                c("cell1", "cell2")))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = mat))
  expect_error(score_signatures(obj), "Ensembl")
  expect_error(score_signatures(obj), "HGNC")
})

test_that("score_signatures detects Entrez IDs", {
  skip_if_not_installed("msigdbr")
  mat <- matrix(1:4, nrow = 2,
                dimnames = list(c("7157", "675"),
                                c("cell1", "cell2")))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = mat))
  expect_error(score_signatures(obj), "Entrez")
  expect_error(score_signatures(obj), "HGNC")
})

test_that("score_signatures adds Hallmark columns", {
  skip_if_not_installed("msigdbr")

  # ctrl = 5 because pbmc3k is a small test dataset with limited gene bins
  obj <- score_signatures(BadranSeq::pbmc3k, collection = "H", ctrl = 5)

  hallmark_cols <- grep("^Hallmark_", colnames(obj@meta.data), value = TRUE)
  expect_true(length(hallmark_cols) > 0)

  for (col in hallmark_cols) {
    expect_true(is.numeric(obj@meta.data[[col]]))
  }

  internal_cols <- grep("^\\.sigscr_", colnames(obj@meta.data), value = TRUE)
  expect_length(internal_cols, 0)
})

test_that("score_signatures respects custom prefix", {
  skip_if_not_installed("msigdbr")

  obj <- score_signatures(BadranSeq::pbmc3k, collection = "H",
                          prefix = "HM_", ctrl = 5)

  hm_cols <- grep("^HM_", colnames(obj@meta.data), value = TRUE)
  expect_true(length(hm_cols) > 0)

  hallmark_cols <- grep("^Hallmark_", colnames(obj@meta.data), value = TRUE)
  expect_length(hallmark_cols, 0)
})

test_that("score_signatures report_only returns tibble", {
  skip_if_not_installed("msigdbr")

  report <- score_signatures(BadranSeq::pbmc3k, collection = "H",
                             report_only = TRUE)

  expect_s3_class(report, "tbl_df")
  expect_named(report, c("gene_set", "gene", "found"))
  expect_type(report$found, "logical")
  expect_true(any(report$found))

  expect_true(any(!report$found))
})

test_that("score_signatures works with non-default collection", {
  skip_if_not_installed("msigdbr")

  obj <- score_signatures(BadranSeq::pbmc3k, collection = "C6", ctrl = 5)

  onco_cols <- grep("^Oncogenic_", colnames(obj@meta.data), value = TRUE)
  expect_true(length(onco_cols) > 0)
})
