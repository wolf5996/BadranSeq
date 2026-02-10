# tests/testthat/test-fetch.R

test_that(".build_cell_scaffold returns tibble with cell_id", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(pbmc3k)
  expect_s3_class(result, "tbl_df")
  expect_true("cell_id" %in% colnames(result))
  expect_equal(nrow(result), ncol(pbmc3k))
  expect_equal(result$cell_id, colnames(pbmc3k))
})

test_that(".build_cell_scaffold adds embeddings when requested", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(pbmc3k, reductions = "umap", dims = 1:2)
  expect_true(all(c("UMAP1", "UMAP2") %in% colnames(result)))
})

test_that(".build_cell_scaffold adds PCA with PC prefix", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(pbmc3k, reductions = "pca", dims = 1:3)
  expect_true(all(c("PC1", "PC2", "PC3") %in% colnames(result)))
})

test_that(".build_cell_scaffold adds multiple reductions", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(
    pbmc3k,
    reductions = c("umap", "pca"),
    dims = 1:2
  )
  expect_true(all(c("UMAP1", "UMAP2", "PC1", "PC2") %in% colnames(result)))
})

test_that(".build_cell_scaffold adds all metadata with TRUE", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(pbmc3k, metadata = TRUE)
  meta_cols <- colnames(pbmc3k@meta.data)
  expect_true(all(meta_cols %in% colnames(result)))
})

test_that(".build_cell_scaffold adds selected metadata columns", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(
    pbmc3k,
    metadata = c("seurat_clusters", "nCount_RNA")
  )
  expect_true(all(c("seurat_clusters", "nCount_RNA") %in% colnames(result)))
  expect_false("nFeature_RNA" %in% colnames(result))
})

test_that(".build_cell_scaffold excludes metadata with FALSE", {
  data(pbmc3k, package = "BadranSeq")
  result <- BadranSeq:::.build_cell_scaffold(pbmc3k, metadata = FALSE)
  expect_equal(colnames(result), "cell_id")
})

test_that(".build_cell_scaffold validates inputs", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    BadranSeq:::.build_cell_scaffold("not_seurat"),
    "object must be a Seurat object"
  )
  expect_error(
    BadranSeq:::.build_cell_scaffold(pbmc3k, reductions = "nonexistent"),
    "not found"
  )
  expect_error(
    BadranSeq:::.build_cell_scaffold(pbmc3k, reductions = "umap", dims = 1:5),
    "exceeds available dimensions"
  )
  expect_error(
    BadranSeq:::.build_cell_scaffold(pbmc3k, metadata = c("fake_column")),
    "not found in metadata"
  )
})

# --- fetch_cell_data() tests ---

test_that("fetch_cell_data returns tibble with cell_id and all metadata", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_cell_data(pbmc3k)
  expect_s3_class(result, "tbl_df")
  expect_true("cell_id" %in% colnames(result))
  expect_equal(nrow(result), ncol(pbmc3k))
  # All metadata columns present
  expect_true(all(colnames(pbmc3k@meta.data) %in% colnames(result)))
  # No embedding columns (reductions = NULL)
  expect_false("UMAP1" %in% colnames(result))
})

test_that("fetch_cell_data includes embeddings when requested", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_cell_data(pbmc3k, reductions = c("umap", "pca"), dims = 1:2)
  expect_true(all(c("UMAP1", "UMAP2", "PC1", "PC2") %in% colnames(result)))
})

test_that("fetch_cell_data respects metadata selection", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_cell_data(pbmc3k, metadata = c("seurat_clusters"))
  expect_true("seurat_clusters" %in% colnames(result))
  expect_false("nCount_RNA" %in% colnames(result))
})

test_that("fetch_cell_data with no metadata and no reductions returns cell_id only", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_cell_data(pbmc3k, metadata = FALSE)
  expect_equal(colnames(result), "cell_id")
})

# --- fetch_feature_data() tests ---

test_that("fetch_feature_data returns tibble with correct structure", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(pbmc3k, features = c("CD3D", "CD8A"))
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("cell_id", "feature", "data") %in% colnames(result)))
  # One row per cell per feature
  expect_equal(nrow(result), ncol(pbmc3k) * 2)
})

test_that("fetch_feature_data feature column is a factor with correct order", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(pbmc3k, features = c("CD8A", "CD3D"))
  expect_s3_class(result$feature, "factor")
  expect_equal(levels(result$feature), c("CD8A", "CD3D"))
})

test_that("fetch_feature_data handles single layer", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(pbmc3k, features = "CD3D", layer = "counts")
  expect_true("counts" %in% colnames(result))
  expect_false("data" %in% colnames(result))
})

test_that("fetch_feature_data handles multiple layers", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(
    pbmc3k,
    features = "CD3D",
    layer = c("counts", "data")
  )
  expect_true(all(c("counts", "data") %in% colnames(result)))
})

test_that("fetch_feature_data includes embeddings when requested", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(
    pbmc3k,
    features = "CD3D",
    reductions = "umap"
  )
  expect_true(all(c("UMAP1", "UMAP2") %in% colnames(result)))
})

test_that("fetch_feature_data includes metadata when requested", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(
    pbmc3k,
    features = "CD3D",
    metadata = c("seurat_clusters")
  )
  expect_true("seurat_clusters" %in% colnames(result))
})

test_that("fetch_feature_data respects assay parameter", {
  data(pbmc3k, package = "BadranSeq")
  result <- fetch_feature_data(
    pbmc3k,
    features = "CD3D",
    assay = "RNA",
    layer = "data",
    metadata = FALSE
  )
  expect_s3_class(result, "tbl_df")
  expect_true("data" %in% colnames(result))
})

test_that("fetch_feature_data warns on missing features and drops them", {
  data(pbmc3k, package = "BadranSeq")
  expect_warning(
    result <- fetch_feature_data(
      pbmc3k,
      features = c("CD3D", "FAKEGENE"),
      metadata = FALSE
    ),
    "not found"
  )
  expect_equal(levels(result$feature), "CD3D")
})

test_that("fetch_feature_data errors when all features missing", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    fetch_feature_data(pbmc3k, features = c("FAKE1", "FAKE2")),
    "No valid features"
  )
})

test_that("fetch_feature_data errors on empty features", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    fetch_feature_data(pbmc3k, features = character(0)),
    "features must be a non-empty character vector"
  )
})

test_that("fetch_feature_data validates assay", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    fetch_feature_data(pbmc3k, features = "CD3D", assay = "nonexistent"),
    "not found"
  )
})

test_that("fetch_feature_data validates layer", {
  data(pbmc3k, package = "BadranSeq")
  expect_error(
    fetch_feature_data(pbmc3k, features = "CD3D", layer = "nonexistent"),
    "not found"
  )
})
