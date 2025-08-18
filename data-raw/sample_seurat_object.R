# Code to prepare `seu_object` dataset goes here

# Check if SeuratData is available
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  stop("SeuratData is required to build package data. Install with:\n",
       "devtools::install_github('satijalab/seurat-data')")
}

# Load the pbmc3k dataset from SeuratData
pbmc3k <- SeuratData::LoadData("pbmc3k")

# normalization
pbmc3k <- Seurat::SCTransform(pbmc3k)

# dimensionality reduction
pbmc3k <- pbmc3k %>%
  Seurat::RunPCA() %>%
  Seurat::RunUMAP(dims = 1:20)

# clustering
pbmc3k <- pbmc3k %>%
  FindNeighbors(
    reduction = "pca",
    dims = 1:20,
    verbose = T
  ) %>%
  FindClusters(
    resolution = 0.5,
    verbose = T
  )

# Save the Seurat object to the package data directory
usethis::use_data(pbmc3k, overwrite = TRUE)
