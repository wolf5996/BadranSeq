# Violin Plot with ggbetweenstats-Style Aesthetics

**\[stable\]**

Creates publication-ready violin plots for one or more features from a
Seurat object, with optional boxplot overlay, jittered raw data points,
and centrality markers (median/mean diamond with numeric label).

For a single feature, returns a bare ggplot2 object. For multiple
features, returns a patchwork composition.

## Usage

``` r
do_ViolinPlot(
  object,
  features,
  group.by = NULL,
  layer = "data",
  assay = NULL,
  colors.use = NULL,
  pt.size = 3,
  pt.alpha = 0.4,
  violin.alpha = 0.2,
  violin.width = 0.5,
  boxplot = TRUE,
  boxplot.width = 0.3,
  boxplot.alpha = 0.2,
  centrality = TRUE,
  centrality.fn = "median",
  ncol = NULL
)
```

## Arguments

- object:

  Seurat object.

- features:

  character. Gene names or numeric metadata column names to plot
  (required). Can mix both types freely — each gets its own panel.
  Examples: `c("CD3D", "nCount_RNA", "percent.mt")`.

- group.by:

  character. Metadata column for grouping. If NULL (default), uses
  active Idents.

- layer:

  character. Layer to extract expression from (default: "data"). Only
  applies to gene features.

- assay:

  character. Assay to use (default: NULL = DefaultAssay).

- colors.use:

  character. Named vector of colors for groups. If NULL, generates
  colors automatically with
  [`generate_badranseq_colors`](https://wolf5996.github.io/BadranSeq/reference/generate_badranseq_colors.md).

- pt.size:

  numeric. Jitter point size (default: 3).

- pt.alpha:

  numeric. Jitter point transparency (default: 0.4).

- violin.alpha:

  numeric. Violin fill alpha (default: 0.2).

- violin.width:

  numeric. Violin width (default: 0.5).

- boxplot:

  logical. Whether to overlay a boxplot (default: TRUE).

- boxplot.width:

  numeric. Boxplot width (default: 0.3).

- boxplot.alpha:

  numeric. Boxplot alpha (default: 0.2).

- centrality:

  logical. Whether to show centrality point and label (default: TRUE).

- centrality.fn:

  character. Centrality function ("median" or "mean", default:
  "median").

- ncol:

  numeric. Number of columns for patchwork layout when plotting multiple
  features (default: NULL = auto).

## Value

A ggplot2 object (single feature) or patchwork object (multiple
features).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single gene violin plot
do_ViolinPlot(seurat_obj, features = "CD3D")

# Multiple genes
do_ViolinPlot(seurat_obj, features = c("CD3D", "CD8A", "CD14"))

# Custom grouping and colors
do_ViolinPlot(seurat_obj, features = "CD3D",
              group.by = "seurat_clusters",
              colors.use = c("0" = "red", "1" = "blue"))

# Without boxplot and centrality
do_ViolinPlot(seurat_obj, features = "CD3D",
              boxplot = FALSE, centrality = FALSE)
} # }
```
