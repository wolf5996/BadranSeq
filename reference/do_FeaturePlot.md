# Enhanced Feature Plot with Viridis Support

Creates feature expression plots overlaid on dimensionality reductions
with viridis color scales and publication-ready styling.

## Usage

``` r
do_FeaturePlot(
  object,
  features,
  reduction = "umap",
  dims = c(1, 2),
  slot = "data",
  use_viridis = TRUE,
  viridis.palette = "B",
  viridis.direction = -1,
  order = TRUE,
  pt.size = 1,
  pt.alpha = 1,
  plot.axes = TRUE,
  variance_digits = 1,
  min.cutoff = NA,
  max.cutoff = NA,
  na.value = "grey90",
  split.by = NULL,
  plot_cell_borders = TRUE,
  border.size = 2,
  border.color = "black",
  ...
)
```

## Arguments

- object:

  Seurat object.

- features:

  character. Gene names or metadata columns to plot.

- reduction:

  character. Reduction to use (default: "umap").

- dims:

  numeric. Vector of 2 PC dimensions to plot (default: c(1, 2)).

- slot:

  character. Data slot ("data", "counts", "scale.data"). Default:
  "data".

- use_viridis:

  logical. Use viridis color scale (default: TRUE).

- viridis.palette:

  character. Viridis palette "A"-"H" (default: "B").

- viridis.direction:

  numeric. Color direction -1 or 1 (default: -1).

- order:

  logical. Plot high-expressing cells on top (default: TRUE).

- pt.size:

  numeric. Point size (default: 1).

- pt.alpha:

  numeric. Point transparency (default: 1).

- plot.axes:

  logical. Show axes (default: TRUE).

- variance_digits:

  numeric. Decimal places for variance (default: 1).

- min.cutoff:

  numeric or "q10". Minimum expression cutoff.

- max.cutoff:

  numeric or "q90". Maximum expression cutoff.

- na.value:

  character. Color for NA/zero values (default: "grey90").

- split.by:

  character. Metadata column for splitting (creates panels with
  silhouettes).

- plot_cell_borders:

  logical. Plot black borders around cells (default: TRUE).

- border.size:

  numeric. Border size multiplier (default: 2).

- border.color:

  character. Border color (default: "black").

- ...:

  Additional arguments (currently unused).

## Value

ggplot2 object or list of ggplot2 objects for multiple features.

## Examples

``` r
if (FALSE) { # \dontrun{
do_FeaturePlot(seurat_object, features = "CD3D")
do_FeaturePlot(seurat_object, features = c("CD3D", "CD8A", "CD4"))
do_FeaturePlot(seurat_object, features = "CD3D", reduction = "pca")
} # }
```
