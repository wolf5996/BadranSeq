# tests/testthat/test-scvi.R

data(pbmc3k, package = "BadranSeq")

# ---------------------------------------------------------------------------
# run_scvi - Input validation
# ---------------------------------------------------------------------------

test_that("run_scvi errors with non-Seurat object", {
  expect_error(
    run_scvi(data.frame(), batch = "stim"),
    "Seurat"
  )
})

test_that("run_scvi errors without reticulate", {
  skip_if(requireNamespace("reticulate", quietly = TRUE),
          "reticulate is installed, skipping missing-package test")
  expect_error(
    run_scvi(pbmc3k),
    "reticulate"
  )
})

test_that("run_scvi errors with invalid assay", {
  skip_if_not_installed("reticulate")
  expect_error(
    run_scvi(pbmc3k, assay = "nonexistent"),
    "not found"
  )
})

test_that("run_scvi errors with invalid layer", {
  skip_if_not_installed("reticulate")
  expect_error(
    run_scvi(pbmc3k, layer = "nonexistent"),
    "not found in assay"
  )
})

test_that("run_scvi errors with invalid batch column", {
  skip_if_not_installed("reticulate")
  expect_error(
    run_scvi(pbmc3k, batch = "nonexistent"),
    "not found in metadata"
  )
})

test_that("run_scvi defaults layer to counts", {
  skip_if_not_installed("reticulate")
  # Should NOT error about layer when counts exists
  # Will error later at Python step, but not at layer validation
  expect_error(
    run_scvi(pbmc3k, conda_env = "fake_env_that_does_not_exist"),
    "fake_env_that_does_not_exist|conda|python|unable",
    ignore.case = TRUE
  )
})
