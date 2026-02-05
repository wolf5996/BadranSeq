# Internal DimPlot Implementation

Core function that builds dimensionality reduction plots using native
ggplot2. When split.by is used, shows the split population in color with
a grey silhouette of all other cells in the background (SCpubr-style).

## Usage

``` r
.do_DimPlot_internal(
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
  label.box = TRUE,
  pt.size = 1,
  pt.alpha = 1,
  plot.axes = TRUE,
  shuffle = TRUE,
  raster = NULL,
  raster.dpi = 300,
  na.value = "grey75",
  plot_cell_borders = TRUE,
  border.size = 2,
  border.color = "black",
  legend.title = NULL,
  ncol = NULL
)
```

## Arguments

- object:

  Seurat object.

- reduction:

  character. Name of reduction to use.

- dims:

  numeric. Vector of 2 dimensions to plot.

- group.by:

  character. Metadata column for grouping.

- split.by:

  character. Metadata column for splitting (creates separate panels with
  silhouettes of other cells).

- colors.use:

  Named vector of colors.

- label:

  logical. Show cluster labels.

- label.size:

  numeric. Label font size.

- label.color:

  character. Label text color.

- label.fill:

  character. Label background fill.

- label.box:

  logical. Use boxed labels.

- pt.size:

  numeric. Point size.

- pt.alpha:

  numeric. Point transparency.

- plot.axes:

  logical. Show axes.

- shuffle:

  logical. Randomize cell plotting order.

- raster:

  logical. Rasterize points for large datasets.

- raster.dpi:

  numeric. DPI for rasterization.

- na.value:

  character. Color for NA/silhouette cells.

- plot_cell_borders:

  logical. Plot black borders around cells.

- border.size:

  numeric. Border size multiplier (default: 2).

- border.color:

  character. Border color (default: "black").

- legend.title:

  character. Custom legend title (NULL for default).

- ncol:

  numeric. Number of columns for split panels (default: NULL, auto).

## Value

ggplot2 object or list of ggplot2 objects when split.by is used.
