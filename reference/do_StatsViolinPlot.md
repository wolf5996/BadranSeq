# Statistical Violin Plot via ggbetweenstats

Creates statistical violin plots for one or more features from a Seurat
object, using
[`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html)
for between-group comparisons with built-in statistical testing.

This function requires the ggstatsplot package to be installed. For
standard violin plots without statistical testing, use
[`do_ViolinPlot`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md)
instead.

For a single feature, returns a bare ggplot2 object. For multiple
features, returns a patchwork composition (requires the patchwork
package).

## Usage

``` r
do_StatsViolinPlot(
  object,
  features,
  group.by = NULL,
  group.levels = NULL,
  layer = "data",
  assay = NULL,
  type = "nonparametric",
  p.adjust.method = "holm",
  pairwise.display = "significant",
  colors.use = NULL,
  ncol = NULL,
  ...
)
```

## Arguments

- object:

  Seurat object.

- features:

  character. Gene names to plot (required). Must be non-empty.

- group.by:

  character. Metadata column for grouping. If NULL (default), uses
  active Idents. Must have 2 or more levels for statistical comparison.

- group.levels:

  character. Subset of group levels to include (default: NULL = all
  levels). Useful for comparing specific clusters or conditions, e.g.
  `group.levels = c("0", "2", "5")`.

- layer:

  character. Layer to extract expression from (default: "data").

- assay:

  character. Assay to use (default: NULL = DefaultAssay).

- type:

  character. Type of statistical test to use. One of "nonparametric"
  (default), "parametric", "robust", or "bayes". Default is
  "nonparametric" (Kruskal-Wallis + Dunn's) because single-cell
  expression data is typically non-normally distributed. Passed through
  to
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html).

- p.adjust.method:

  character. Method for p-value adjustment (default: "holm"). Passed
  through to
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html).

- pairwise.display:

  character. Which pairwise comparisons to display. One of "significant"
  (default), "non-significant", "all", or "none". Passed through to
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html).

- colors.use:

  character. Named vector of colors for groups. If NULL, generates
  colors automatically with
  [`generate_badranseq_colors`](https://wolf5996.github.io/BadranSeq/reference/generate_badranseq_colors.md).

- ncol:

  numeric. Number of columns for patchwork layout when plotting multiple
  features (default: NULL = auto).

- ...:

  Additional arguments passed to
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html).

## Value

A ggplot2 object (single feature) or patchwork object (multiple
features).

## Examples

``` r
if (FALSE) { # \dontrun{
# Single gene statistical violin plot
do_StatsViolinPlot(seurat_obj, features = "CD3D",
                   group.by = "seurat_clusters")

# Multiple genes
do_StatsViolinPlot(seurat_obj, features = c("CD3D", "CD8A"),
                   group.by = "seurat_clusters")

# Compare specific clusters only
do_StatsViolinPlot(seurat_obj, features = "CD3D",
                   group.by = "seurat_clusters",
                   group.levels = c("0", "2", "5"))

# Parametric test
do_StatsViolinPlot(seurat_obj, features = "CD3D",
                   group.by = "condition", type = "parametric")
} # }
```
