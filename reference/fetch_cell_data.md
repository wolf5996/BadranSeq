# Fetch Cell-Level Data from a Seurat Object

**\[stable\]**

Extracts a tidy tibble with one row per cell containing cell
identifiers, optional dimensionality reduction embeddings, and optional
metadata. No expression data is included — use
[`fetch_feature_data`](https://wolf5996.github.io/BadranSeq/reference/fetch_feature_data.md)
for that.

## Usage

``` r
fetch_cell_data(object, reductions = NULL, dims = 1:2, metadata = TRUE)
```

## Arguments

- object:

  Seurat object.

- reductions:

  character. Reduction names to include (default: NULL = none). Example:
  `c("umap", "pca")`.

- dims:

  integer. Dimensions to extract per reduction (default: 1:2).

- metadata:

  logical or character. TRUE = all metadata columns, FALSE = none, or a
  character vector of specific column names.

## Value

A tibble with columns: cell_id, optional embedding columns, and optional
metadata columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# All metadata, no embeddings
fetch_cell_data(seurat_obj)

# With UMAP coordinates
fetch_cell_data(seurat_obj, reductions = "umap")

# Specific metadata + PCA
fetch_cell_data(seurat_obj, reductions = "pca", dims = 1:10,
                metadata = c("seurat_clusters", "condition"))
} # }
```
