# tests/testthat/test-convert.R

data(pbmc3k, package = "BadranSeq")

# ---------------------------------------------------------------------------
# seurat_to_h5ad - Input validation
# ---------------------------------------------------------------------------

test_that("seurat_to_h5ad errors with non-Seurat object", {
  expect_error(
    seurat_to_h5ad(data.frame()),
    "Seurat"
  )
})

test_that("seurat_to_h5ad errors without anndata R package", {
  skip_if(requireNamespace("anndata", quietly = TRUE),
          "anndata is installed, skipping missing-package test")
  skip_if_not_installed("reticulate")
  expect_error(
    seurat_to_h5ad(pbmc3k),
    "anndata"
  )
})

test_that("seurat_to_h5ad errors with invalid assay", {
  skip_if_not_installed("anndata")
  expect_error(
    seurat_to_h5ad(pbmc3k, assay = "nonexistent"),
    "not found"
  )
})

test_that("seurat_to_h5ad errors with invalid layer", {
  skip_if_not_installed("anndata")
  expect_error(
    seurat_to_h5ad(pbmc3k, layers = c("counts", "nonexistent")),
    "not found"
  )
})

# ---------------------------------------------------------------------------
# seurat_to_h5ad - Conversion (needs anndata R + Python)
# ---------------------------------------------------------------------------

test_that("seurat_to_h5ad returns AnnData object", {
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  skip_if_not(
    reticulate::py_module_available("anndata"),
    "Python anndata not available"
  )

  adata <- seurat_to_h5ad(pbmc3k, assay = "SCT", layers = "counts")
  expect_true(inherits(adata, "anndata._core.anndata.AnnData") ||
              inherits(adata, "python.builtin.object") ||
              inherits(adata, "AnnDataR6"))
})

test_that("seurat_to_h5ad writes h5ad file when outfile given", {
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  skip_if_not(
    reticulate::py_module_available("anndata"),
    "Python anndata not available"
  )

  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)

  result <- seurat_to_h5ad(pbmc3k, outfile = tmp, assay = "SCT",
                           layers = "counts")
  expect_true(file.exists(tmp))
  expect_true(file.size(tmp) > 0)
})

test_that("seurat_to_h5ad transfers reductions to obsm", {
  skip_if_not_installed("anndata")
  skip_if_not_installed("reticulate")
  skip_if_not(
    reticulate::py_module_available("anndata"),
    "Python anndata not available"
  )

  adata <- seurat_to_h5ad(pbmc3k, assay = "SCT", layers = "counts")
  obsm_keys <- adata$obsm_keys()
  expect_true("X_pca" %in% obsm_keys)
  expect_true("X_umap" %in% obsm_keys)
})
