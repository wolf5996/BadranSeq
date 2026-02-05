# Generalized Dimensionality Reduction Plot

Unified wrapper for dimensionality reduction plots. Routes to
specialized functions for PCA and UMAP, handles other reductions
directly. When split.by is used, shows each split category with a grey
silhouette of all other cells.

## Usage

``` r
do_DimPlot(
  object,
  reduction = "umap",
  dims = c(1, 2),
  group.by = NULL,
  split.by = NULL,
  colors.use = NULL,
  label = TRUE,
  label.size = 4,
  label.color = "black",
  label.fill = "white",
  pt.size = 1,
  pt.alpha = 1,
  plot.axes = TRUE,
  variance_digits = 1,
  shuffle = TRUE,
  plot_cell_borders = TRUE,
  border.size = 2,
  border.color = "black",
  legend.title = NULL,
  ncol = NULL,
  ...
)
```

## Arguments

- object:

  Seurat object.

- reduction:

  character. Reduction to use (default: "umap").

- dims:

  numeric. Vector of 2 PC dimensions to plot (default: c(1, 2)).

- group.by:

  character. Metadata column for grouping (default: active identity).

- split.by:

  character. Metadata column for splitting (creates panels with
  silhouettes).

- colors.use:

  Named vector of colors for groups.

- label:

  logical. Show cluster labels (default: TRUE).

- label.size:

  numeric. Label font size (default: 4).

- label.color:

  character. Label text color (default: "black").

- label.fill:

  character. Label background fill (default: "white").

- pt.size:

  numeric. Point size (default: 1).

- pt.alpha:

  numeric. Point transparency (default: 1).

- plot.axes:

  logical. Show axes (default: TRUE).

- variance_digits:

  numeric. Decimal places for variance (default: 1).

- shuffle:

  logical. Randomize cell plotting order (default: TRUE).

- plot_cell_borders:

  logical. Plot black borders around cells (default: TRUE).

- border.size:

  numeric. Border size multiplier (default: 2).

- border.color:

  character. Border color (default: "black").

- legend.title:

  character. Custom legend title (NULL removes title).

- ncol:

  numeric. Number of columns for split panels (default: NULL, auto).

- ...:

  Additional arguments (currently unused).

## Value

ggplot2 object with dimensionality reduction plot.

## Examples

``` r
if (FALSE) { # \dontrun{
do_DimPlot(seurat_object, group.by = "condition")
do_DimPlot(seurat_object, reduction = "pca", group.by = "cluster")
do_DimPlot(seurat_object, split.by = "sample")  # Shows silhouettes
} # }
```
