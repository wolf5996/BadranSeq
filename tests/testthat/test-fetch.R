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
