# Enhanced Elbow Plot for PCA Dimensionality

**\[stable\]**

An improved version of Seurat's ElbowPlot with enhanced aesthetics,
customizable titles, optional cutoff lines, and variance explained
options.

## Usage

``` r
EnhancedElbowPlot(
  object,
  ndims = 50,
  reduction = "pca",
  cutoff_pc = NULL,
  variance = TRUE,
  title = NULL,
  point_size = 3,
  point_color = "steelblue",
  line_color = "steelblue",
  cutoff_color = "red",
  base_size = 14
)
```

## Arguments

- object:

  Seurat object.

- ndims:

  numeric. Number of dimensions to plot (default: 50).

- reduction:

  character. Reduction to use (default: "pca").

- cutoff_pc:

  numeric. PC number to mark with vertical line.

- variance:

  logical. Plot variance explained instead of stdev (default: TRUE).

- title:

  character. Custom title.

- point_size:

  numeric. Point size (default: 3).

- point_color:

  character. Point color (default: "steelblue").

- line_color:

  character. Line color (default: "steelblue").

- cutoff_color:

  character. Cutoff line color (default: "red").

- base_size:

  numeric. Base font size (default: 14).

## Value

ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
EnhancedElbowPlot(seurat_obj)
EnhancedElbowPlot(seurat_obj, cutoff_pc = 10)
EnhancedElbowPlot(seurat_obj, variance = FALSE)
} # }
```
