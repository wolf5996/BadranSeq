# Sankey Plot for Categorical Metadata

Creates a Sankey (alluvial) diagram showing the relationship between 2–3
categorical metadata variables from a Seurat object. Flows are coloured
by the first variable in `columns`. Cell counts are displayed inside
each stratum by default.

Requires the ggalluvial package to be installed.

## Usage

``` r
do_SankeyPlot(
  object,
  columns,
  colors.use = NULL,
  label = TRUE,
  label.size = 3,
  alpha = 0.5,
  curve_type = "sigmoid"
)
```

## Arguments

- object:

  Seurat object.

- columns:

  character. Vector of 2–3 metadata column names whose relationships to
  visualise (required).

- colors.use:

  character. Named vector of colours for the first column in `columns`.
  If NULL, generates colours automatically with
  [`generate_badranseq_colors`](https://wolf5996.github.io/BadranSeq/reference/generate_badranseq_colors.md).

- label:

  logical. Whether to display cell counts inside strata (default: TRUE).

- label.size:

  numeric. Text size for stratum labels (default: 3).

- alpha:

  numeric. Flow transparency (default: 0.5).

- curve_type:

  character. Curve type for flows (default: "sigmoid"). Passed to
  [`ggalluvial::geom_alluvium()`](http://corybrunson.github.io/ggalluvial/reference/geom_alluvium.md).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Relationship between cell type and cluster
do_SankeyPlot(seurat_obj, columns = c("seurat_annotations", "seurat_clusters"))

# Three-way relationship
do_SankeyPlot(seurat_obj,
              columns = c("seurat_annotations", "condition", "Phase"))

# Custom colours
do_SankeyPlot(seurat_obj,
              columns = c("seurat_annotations", "seurat_clusters"),
              colors.use = c("B" = "blue", "T" = "red"))
} # }
```
