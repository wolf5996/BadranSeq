# Boil Plot for Zero-Inflated Gene Expression

**\[experimental\]**

Creates boil plots for one or more genes from a Seurat object.
Single-cell expression is dominated by zeros, which squash conventional
violins towards the baseline. A boil plot separates the two questions a
violin mixes together:

- **How many cells express the gene?** A detection bar below zero shows
  the fraction of cells with non-zero expression (grey track = all
  cells, coloured fill = expressing cells, percentage label underneath).

- **How is the gene expressed in those cells?** The violin, jittered
  points, centrality marker and optional boxplot (`boxplot = TRUE`)
  above zero use non-zero cells only.

Geometry and colours follow the original boil plot, drawn with
BadranSeq's palette and
[`theme_badranseq`](https://wolf5996.github.io/BadranSeq/reference/theme_badranseq.md).
For a single gene, returns a bare ggplot2 object. For multiple genes,
returns a patchwork composition.

The boil plot design comes from the boilplot package by Lara Jerman
(<https://github.com/laraj1/boilplot-r>, Python original at
<https://github.com/laraj1/boilplot-python>). This function is a
reimplementation in BadranSeq's visual style, with credit to the
original as the source of inspiration.

## Usage

``` r
do_BoilPlot(
  object,
  features,
  group.by = NULL,
  group.levels = NULL,
  layer = "data",
  assay = NULL,
  colors.use = NULL,
  min.nonzero = 10,
  pt.size = 0.7,
  pt.alpha = 0.6,
  jitter.width = 0.1,
  violin.alpha = 0.8,
  violin.width = 0.9,
  boxplot = FALSE,
  boxplot.width = 0.3,
  boxplot.alpha = 0.2,
  centrality = TRUE,
  centrality.fn = "median",
  ncol = 1,
  nrow = NULL
)
```

## Arguments

- object:

  Seurat object.

- features:

  character. Gene names to plot (required). Metadata columns are not
  accepted; use
  [`do_ViolinPlot`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md)
  for those.

- group.by:

  character. Metadata column for grouping. If NULL (default), uses
  active Idents.

- group.levels:

  character. Subset of group levels to include, in the order given
  (default: NULL = all levels in their existing order).

- layer:

  character. Layer to extract expression from (default: "data"). Use a
  non-negative layer such as "data" or "counts".

- assay:

  character. Assay to use (default: NULL = DefaultAssay).

- colors.use:

  character. Named vector of colors for groups. If NULL, generates
  colors automatically with
  [`generate_badranseq_colors`](https://wolf5996.github.io/BadranSeq/reference/generate_badranseq_colors.md).

- min.nonzero:

  numeric. Minimum number of expressing cells a group needs to appear in
  a gene's panel (default: 10). Groups below the threshold are removed
  from that panel entirely, including their detection bar and x-axis
  position. Genes with no qualifying group are skipped with a warning.

- pt.size:

  numeric. Jitter point size (default: 0.7).

- pt.alpha:

  numeric. Jitter point transparency (default: 0.6).

- jitter.width:

  numeric. Horizontal jitter width for points (default: 0.1).

- violin.alpha:

  numeric. Violin fill alpha (default: 0.8).

- violin.width:

  numeric. Maximum violin width (default: 0.9). All violins share this
  width regardless of group size, so small groups are not squashed.

- boxplot:

  logical. Whether to overlay a boxplot (default: FALSE).

- boxplot.width:

  numeric. Boxplot width (default: 0.3).

- boxplot.alpha:

  numeric. Boxplot alpha (default: 0.2).

- centrality:

  logical. Whether to show centrality point and label for non-zero cells
  (default: TRUE).

- centrality.fn:

  character. Centrality function ("median" or "mean", default:
  "median").

- ncol:

  numeric. Number of columns for patchwork layout when plotting multiple
  genes (default: 1, one gene per row so panels stack vertically).

- nrow:

  numeric. Number of rows for patchwork layout (default: NULL = auto).
  Set `ncol = NULL` together with `nrow` to control rows only.

## Value

A ggplot2 object (single gene) or patchwork object (multiple genes).

## References

Jerman, L. boilplot: faithful and interpretable visualisations of
zero-inflated distributions. <https://github.com/laraj1/boilplot-r>

## See also

[`do_ViolinPlot`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md),
[`do_StatsViolinPlot`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(pbmc3k)

# Single gene
do_BoilPlot(pbmc3k, features = "CD8A", group.by = "seurat_annotations")

# Multiple genes, selected groups in a chosen order
do_BoilPlot(pbmc3k, features = c("GZMA", "NKG7", "PRF1"),
            group.by = "seurat_annotations",
            group.levels = c("NK", "CD8 T", "Memory CD4 T"))

# Mean instead of median, with boxplot overlay
do_BoilPlot(pbmc3k, features = "MS4A1", centrality.fn = "mean",
            boxplot = TRUE)
} # }
```
