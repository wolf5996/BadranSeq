# Publication-Ready GSEA Results Barplot

**\[stable\]**

Creates a horizontal barplot of top pathways from a Gene Set Enrichment
Analysis (GSEA) result. Instead of the default dotplot, this function
shows Normalized Enrichment Scores (NES) as bars with adjusted p-values
as text labels and leading-edge gene annotations.

## Usage

``` r
do_GseaPlot(
  gsea_object,
  analysis_name = "GSEA Analysis",
  simplify_go = TRUE,
  number_of_pathways = 10,
  selection = NULL,
  fill_scale = c("diverging", "viridis"),
  show_padj = TRUE,
  show_leading_edge = TRUE,
  title = NULL
)
```

## Arguments

- gsea_object:

  A `gseaResult` object from
  [`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html),
  `ReactomePA::gsePathway()`, or similar, containing results in the
  `@result` slot.

- analysis_name:

  Character. Subtitle for the plot. Defaults to `"GSEA Analysis"`.

- simplify_go:

  Logical. Whether to apply
  [`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html)
  to GO terms (default: `TRUE`). Only applies when GO term IDs are
  detected; silently skipped otherwise.

- number_of_pathways:

  Integer. Number of top upregulated and downregulated pathways to
  display (default: `10`). Ignored when `selection` is provided.

- selection:

  Character vector of pathway descriptions to plot. Overrides automatic
  top/bottom selection. Defaults to `NULL`.

- fill_scale:

  Character. Color scale: `"diverging"` (RdBu, 0-centered) or
  `"viridis"` (sequential viridis-C, original style). Defaults to
  `"diverging"` to emphasize the signed nature of NES.

- show_padj:

  Logical. Show adjusted p-value labels on bars. Defaults to `TRUE`.

- show_leading_edge:

  Logical. Show leading-edge gene count and gene symbols in y-axis
  labels. Defaults to `TRUE`.

- title:

  Character. Custom plot title. If `NULL`, auto-generated from
  `analysis_name`.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
data(gsea_pbmc3k)
do_GseaPlot(gsea_pbmc3k, analysis_name = "CD14+ Mono vs Others")

# Use viridis (original style)
do_GseaPlot(gsea_pbmc3k, fill_scale = "viridis")

# Plot specific pathways
do_GseaPlot(
  gsea_pbmc3k,
  selection = c("cell chemotaxis")
)
} # }
```
