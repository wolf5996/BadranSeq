# Fetch Feature Expression Data from a Seurat Object

Extracts a tidy tibble with one row per cell per feature. Expression
values are provided as one column per requested layer. Optionally
includes embedding coordinates and metadata, joined by cell_id.

## Usage

``` r
fetch_feature_data(
  object,
  features,
  layer = "data",
  assay = NULL,
  reductions = NULL,
  dims = 1:2,
  metadata = TRUE
)
```

## Arguments

- object:

  Seurat object.

- features:

  character. Gene names to extract (required). Becomes a factor column
  preserving the input order.

- layer:

  character. Layer(s) to extract expression from. One of "counts",
  "data", "scale.data", or a combination (default: "data").

- assay:

  character. Assay to use (default: NULL = DefaultAssay).

- reductions:

  character. Reduction names to include (default: NULL = none).

- dims:

  integer. Dimensions to extract per reduction (default: 1:2).

- metadata:

  logical or character. TRUE = all metadata columns, FALSE = none, or a
  character vector of specific column names.

## Value

A tibble with columns: cell_id, features (factor), one column per layer,
optional embedding columns, optional metadata columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage -- normalized expression
fetch_feature_data(seurat_obj, features = c("CD3D", "CD8A"))

# Multiple layers
fetch_feature_data(seurat_obj, features = "CD3D",
                   layer = c("counts", "data"))

# With embeddings and specific metadata
fetch_feature_data(seurat_obj, features = c("CD3D", "CD14"),
                   reductions = "umap",
                   metadata = c("seurat_clusters", "condition"))

# Explicit assay
fetch_feature_data(seurat_obj, features = "CD3D",
                   assay = "RNA", layer = "counts")
} # }
```
