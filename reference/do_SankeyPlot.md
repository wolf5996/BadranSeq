# Sankey Plot for Categorical Metadata

**\[stable\]**

Creates a Sankey (alluvial) diagram showing the relationship between 2–3
categorical metadata variables from a Seurat object. Each stratum is
labelled and coloured; flows connect cells across variables.

Requires the ggalluvial package to be installed.

## Usage

``` r
do_SankeyPlot(object, columns, label.size = 3)
```

## Arguments

- object:

  Seurat object.

- columns:

  character. Vector of 2–3 metadata column names whose relationships to
  visualise (required).

- label.size:

  numeric. Text size for stratum labels (default: 3).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Relationship between cell type and cluster
do_SankeyPlot(seurat_obj,
              columns = c("seurat_annotations", "seurat_clusters"))

# Three-way relationship
do_SankeyPlot(seurat_obj,
              columns = c("seurat_annotations", "condition", "Phase"))
} # }
```
