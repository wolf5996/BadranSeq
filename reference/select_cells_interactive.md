# Interactive Cell Selection for Seurat Objects

**\[stable\]**

## Usage

``` r
select_cells_interactive(
  object,
  reduction = c("umap", "tsne", "pca"),
  color_by = "Idents",
  return_cells = FALSE
)
```

## Arguments

- object:

  Seurat object with dimensionality reduction computed.

- reduction:

  character. Reduction to use ("umap", "tsne", "pca").

- color_by:

  character. Metadata column for coloring (default: "Idents").

- return_cells:

  logical. Return cell names instead of Seurat subset.

## Value

Subsetted Seurat object or character vector of cell names.

## Examples

``` r
if (FALSE) { # \dontrun{
seu_subset <- select_cells_interactive(seurat_obj)
cell_names <- select_cells_interactive(seurat_obj, return_cells = TRUE)
} # }
```
