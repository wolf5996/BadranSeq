# Publication-Ready Single-Cell Visualisation

## Motivation

Seurat’s plotting functions are designed for exploration. Getting them
to publication quality typically requires layering `ggplot2` theme
adjustments, adding cell borders, repositioning legends, and
reformatting axis labels — the same boilerplate repeated across every
figure in a manuscript.

BadranSeq wraps this into a small set of functions with opinionated
defaults: cell borders on, cluster labels on, clean white theme, legend
at bottom. The goal is a plot that looks good in a paper without
additional customisation.

## From Seurat defaults to publication-ready

The quickest way to see the difference — the same data plotted with
[`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html)
and
[`BadranSeq::do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md):

``` r

p1 <- Seurat::DimPlot(pbmc3k, reduction = "umap") +
  ggtitle("Seurat::DimPlot()")

p2 <- do_UmapPlot(pbmc3k) +
  ggtitle("BadranSeq::do_UmapPlot()")

p1 + p2
```

![](visualisation_files/figure-html/comparison-1.png)

BadranSeq adds cell borders, boxed cluster labels, and a clean theme
automatically. No
[`theme()`](https://ggplot2.tidyverse.org/reference/theme.html) calls
needed.

## UMAP plots

[`do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md)
is the workhorse for UMAP visualisation.

### Basic usage

``` r

do_UmapPlot(pbmc3k)
```

![](visualisation_files/figure-html/umap-basic-1.png)

### Group by a metadata column

``` r

do_UmapPlot(pbmc3k, group.by = "seurat_annotations")
```

![](visualisation_files/figure-html/umap-groupby-1.png)

### Group by multiple metadata columns

Pass a vector of column names to create side-by-side panels:

``` r

do_UmapPlot(pbmc3k, group.by = c("seurat_clusters", "seurat_annotations"))
```

![](visualisation_files/figure-html/umap-multi-groupby-1.png)

### Control panel layout

When combining multiple panels, use `ncol` and `nrow` to control the
patchwork layout. For example, stack the panels vertically with
`ncol = 1`:

``` r

do_UmapPlot(pbmc3k,
            group.by = c("seurat_clusters", "seurat_annotations"),
            ncol = 1)
```

![](visualisation_files/figure-html/umap-layout-control-1.png)

### Hide the legend

``` r

do_UmapPlot(pbmc3k, show.legend = FALSE)
```

![](visualisation_files/figure-html/umap-no-legend-1.png)

### Adjust point size and remove labels

``` r

do_UmapPlot(pbmc3k, pt.size = 0.3, label = FALSE)
```

![](visualisation_files/figure-html/umap-customise-1.png)

### Remove cell borders

``` r

do_UmapPlot(pbmc3k, plot_cell_borders = FALSE)
```

![](visualisation_files/figure-html/umap-no-borders-1.png)

## PCA with automatic variance labels

[`do_PcaPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_PcaPlot.md)
does what Seurat’s PCA plot does not — it calculates variance explained
per PC and formats the axis labels as “PC1 (X.X%)”.

``` r

p1 <- Seurat::DimPlot(pbmc3k, reduction = "pca") +
  ggtitle("Seurat::DimPlot(reduction = 'pca')")

p2 <- do_PcaPlot(pbmc3k) +
  ggtitle("BadranSeq::do_PcaPlot()")

p1 + p2
```

![](visualisation_files/figure-html/pca-comparison-1.png)

The variance annotation is automatic — no manual calculation or label
formatting required.

### Different dimension pairs

Explore PC2 vs PC3:

``` r

do_PcaPlot(pbmc3k, dims = c(2, 3))
```

![](visualisation_files/figure-html/pca-dims-1.png)

## The unified entry point

[`do_DimPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_DimPlot.md)
routes to the appropriate handler based on the `reduction` argument:

| `reduction` | Routes to | Axis labels |
|----|----|----|
| `"umap"` | [`do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md) | UMAP1, UMAP2 |
| `"pca"` | [`do_PcaPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_PcaPlot.md) | PC1 (X.X%), PC2 (X.X%) |
| other | [`.do_DimPlot_internal()`](https://wolf5996.github.io/BadranSeq/reference/dot-do_DimPlot_internal.md) | REDUCTION1, REDUCTION2 |

``` r

p1 <- do_DimPlot(pbmc3k, reduction = "umap") + ggtitle("reduction = 'umap'")
p2 <- do_DimPlot(pbmc3k, reduction = "pca") + ggtitle("reduction = 'pca'")

p1 + p2
```

![](visualisation_files/figure-html/dimplot-routing-1.png)

Use
[`do_DimPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_DimPlot.md)
when writing general-purpose code that should work with any reduction.
Use the specific functions
([`do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md),
[`do_PcaPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_PcaPlot.md))
when you know the reduction type.

## Feature expression plots

[`do_FeaturePlot()`](https://wolf5996.github.io/BadranSeq/reference/do_FeaturePlot.md)
overlays gene expression on the embedding using a viridis colour scale.

``` r

do_FeaturePlot(pbmc3k, features = c("CD3D", "CD14", "MS4A1"))
```

![](visualisation_files/figure-html/featureplot-1.png)

Cells are ordered by expression by default (`order = TRUE`), so
high-expressing cells are plotted on top — avoiding the common problem
of rare populations being hidden beneath low-expressing cells.

### Control panel layout

Use `ncol` and `nrow` to arrange multi-feature panels. For example,
stack them vertically:

``` r

do_FeaturePlot(pbmc3k,
               features = c("CD3D", "CD14", "MS4A1"),
               ncol = 1)
```

![](visualisation_files/figure-html/featureplot-layout-1.png)

## Split-panel silhouettes

This is BadranSeq’s signature visualisation feature. When you pass
`split.by`, each panel shows:

1.  A **grey silhouette** of all cells (spatial context)
2.  **Black borders** on all cells
3.  **Coloured cells** only for the current split category

This preserves the spatial layout across panels, making it easy to see
where a subpopulation sits relative to the full dataset.

First, add a synthetic condition to demonstrate:

``` r

set.seed(42)
pbmc3k$condition <- sample(
  c("Control", "Treatment"),
  ncol(pbmc3k),
  replace = TRUE
)
```

Now compare Seurat’s default split with BadranSeq’s silhouette approach:

``` r

p1 <- Seurat::DimPlot(pbmc3k, reduction = "umap", split.by = "condition") +
  ggtitle("Seurat::DimPlot(split.by)")

p1
```

![](visualisation_files/figure-html/split-comparison-1.png)

``` r

do_UmapPlot(pbmc3k, split.by = "condition")
```

![](visualisation_files/figure-html/split-silhouette-1.png)

In the Seurat version, each panel only shows its own cells — you lose
the sense of where those cells sit in the overall structure. In
BadranSeq’s version, the grey background preserves this context: you can
immediately see which clusters are enriched in each condition.

The silhouette approach works with any grouping:

``` r

do_UmapPlot(pbmc3k, group.by = "seurat_annotations", split.by = "condition")
```

![](visualisation_files/figure-html/split-annotations-1.png)

## Violin plots

BadranSeq provides two violin plot functions built on
[`fetch_feature_data()`](https://wolf5996.github.io/BadranSeq/reference/fetch_feature_data.md):

- **[`do_StatsViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md)**
  — wraps
  [`ggstatsplot::ggbetweenstats()`](https://www.indrapatil.com/ggstatsplot/reference/ggbetweenstats.html)
  for between-group comparisons with built-in statistical testing.
- **[`do_ViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md)**
  — a pure ggplot2 violin with boxplot overlay, jittered points, and
  centrality markers (no statistical testing).

The primary function is
[`do_StatsViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md).
If you want the same aesthetic without statistics, use
[`do_ViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md)
or set `pairwise.display = "none"` in
[`do_StatsViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md).

### Statistical violin plots

[`do_StatsViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md)
runs a nonparametric omnibus test (Kruskal–Wallis) and pairwise
comparisons (Dunn’s test) by default, annotating the plot with p-values:

``` r

do_StatsViolinPlot(pbmc3k, features = "CD3D",
                   group.by = "seurat_annotations",
                   pairwise.display = "none")
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](visualisation_files/figure-html/stats-violin-basic-1.png)

With nine annotation groups, pairwise brackets would overwhelm the
figure. The `pairwise.display` argument controls this:

### Controlling pairwise comparisons

| Value               | Behaviour                        |
|---------------------|----------------------------------|
| `"significant"`     | Only significant pairs (default) |
| `"non-significant"` | Only non-significant pairs       |
| `"all"`             | Every pairwise comparison        |
| `"none"`            | Omnibus test only, no brackets   |

Setting `pairwise.display = "none"` keeps the subtitle showing the
overall test statistic while removing all bracket annotations — useful
for publication figures where the brackets would overlap.

### Comparing specific groups with `group.levels`

Often you do not want to compare all cell types — perhaps only a few are
biologically relevant. The `group.levels` argument filters the grouping
variable to a specified subset *before* the statistical test is run, so
the test is performed only on the groups you care about.

With fewer groups, pairwise brackets become readable again. S100A4 rises
from naive to memory T cells, so the three groups have clearly different
medians:

``` r

do_StatsViolinPlot(pbmc3k, features = "S100A4",
                   group.by = "seurat_annotations",
                   group.levels = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](visualisation_files/figure-html/stats-violin-group-levels-1.png)

This is cleaner than subsetting the Seurat object manually and avoids
the need to drop unused factor levels.

### Parametric and Bayesian alternatives

The default nonparametric test is appropriate for most scRNA-seq data,
but you can switch to parametric (ANOVA / *t*-test), robust, or Bayesian
tests via the `type` argument:

``` r

do_StatsViolinPlot(pbmc3k, features = "S100A4",
                   group.by = "seurat_annotations",
                   group.levels = c("Naive CD4 T", "Memory CD4 T", "CD8 T"),
                   type = "parametric")
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](visualisation_files/figure-html/stats-violin-parametric-1.png)

### Multiple genes

Pass a vector of gene names to get a combined patchwork layout:

``` r

do_StatsViolinPlot(pbmc3k,
                   features = c("CD3D", "CD8A"),
                   group.by = "seurat_annotations",
                   pairwise.display = "none",
                   ncol = 2)
```

    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.
    ## Scale for colour is already present.
    ## Adding another scale for colour, which will replace the existing scale.

![](visualisation_files/figure-html/stats-violin-multi-1.png)

### Plain violin plots

[`do_ViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md)
produces a violin + boxplot + jitter + centrality marker (median
diamond) without any statistical annotation. Use this when you want a
clean descriptive plot:

``` r

do_ViolinPlot(pbmc3k, features = "CD3D",
              group.by = "seurat_annotations")
```

![](visualisation_files/figure-html/plain-violin-1.png)

## Boil plots for zero-inflated genes

Most cells do not express most genes. In a conventional violin the zeros
pile up at the baseline, the cells that *do* express the gene are
squeezed into a thin tail, and the fraction of expressing cells has to
be guessed from the thickness of the base.
[`do_BoilPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_BoilPlot.md)
separates the two questions:

- **How many cells express the gene?** A detection bar below zero shows
  the percentage of cells with non-zero expression per group.
- **How is it expressed in those cells?** The violin, points and median
  diamond above zero use non-zero cells only.

The design comes from the
[boilplot](https://github.com/laraj1/boilplot-r) package by Lara Jerman;
BadranSeq redraws it with its own palette and theme and adds the
controls below.

### A single gene

``` r

do_BoilPlot(pbmc3k, features = "CCL5", group.by = "seurat_annotations")
```

![](visualisation_files/figure-html/boil-basic-1.png)

CCL5 is detected in nearly every CD8 T cell and platelet but in only a
fifth of naive CD4 T cells, and among the cells that do express it the
median rises from about one count in monocytes and B cells to well over
ten in CD8 T cells. Both facts are read directly from the panel; in a
conventional violin the zeros would hide the second one.

A median label of 0.69 is `log1p(1)`: most expressing cells in that
group carry a single corrected count. It appears often for sparsely
detected genes and is a property of the data, not of the plot.

### Filtering out groups that barely express the gene

Boil plots carry a built-in filter. A group is kept in a gene’s panel
only if at least `min.nonzero` of its cells express that gene; every
other group disappears from the panel entirely, taking its violin,
points, detection bar and x-axis slot with it. The default is 10
expressing cells. The reason is statistical: a density estimated from a
handful of points is noise, and a detection bar for a group that does
not express the gene says nothing a reader cannot infer from its
absence.

The filter is a **count of expressing cells**, not a percentage. A group
of 700 cells at 2% detection (14 cells) passes the default; a group of
30 cells at 20% (6 cells) does not. This keeps small but genuinely
expressing populations out of the plot only when they are too small to
draw.

With `min.nonzero = 1`, nothing is filtered and every group that has at
least one expressing cell appears:

``` r

do_BoilPlot(pbmc3k, features = "CCL5", group.by = "seurat_annotations",
            min.nonzero = 1)
```

    ## Warning: Groups with fewer than two datapoints have been dropped.
    ## ℹ Set `drop = FALSE` to consider such groups for position adjustment purposes.

![](visualisation_files/figure-html/boil-min-nonzero-all-1.png)

Raising the threshold strips the panel down to the groups where the gene
is actually expressed. At 100 expressing cells only the four
CCL5-positive populations remain:

``` r

do_BoilPlot(pbmc3k, features = "CCL5", group.by = "seurat_annotations",
            min.nonzero = 100)
```

![](visualisation_files/figure-html/boil-min-nonzero-1.png)

Two details matter when reading a filtered panel:

- The percentages of the remaining groups are unchanged. They are still
  computed over all cells in the group, including the non-expressing
  ones, so a bar reading 32% means 32% of all Memory CD4 T cells, not
  32% of the cells shown.
- Filtering is per gene. In a multi-gene plot each panel keeps the
  groups that pass for *that* gene, so the x-axes can differ between
  panels. A gene with no group above the threshold is skipped with a
  warning; if no gene qualifies the call stops with an error.

To force the same groups into every panel regardless of expression, use
`group.levels` instead of (or together with) `min.nonzero`; groups
listed in `group.levels` are still subject to the threshold, so set
`min.nonzero = 1` if you want them all drawn.

### Multiple genes

Several genes stack vertically by default (`ncol = 1`), one panel per
gene. Each panel keeps only the groups that pass `min.nonzero` for that
gene, so panels can show different group sets — compare the x-axes of
the three panels below:

``` r

do_BoilPlot(pbmc3k, features = c("CCL5", "HLA-DRA", "IL32"),
            group.by = "seurat_annotations")
```

![](visualisation_files/figure-html/boil-multi-1.png)

Use `ncol` and `nrow` to change the layout, as with the violin
functions.

### Selecting and ordering groups

`group.levels` works as in
[`do_StatsViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md):
it keeps only the listed groups, in the order given.

``` r

do_BoilPlot(pbmc3k, features = "NKG7", group.by = "seurat_annotations",
            group.levels = c("NK", "CD8 T", "Memory CD4 T", "B"))
```

![](visualisation_files/figure-html/boil-group-levels-1.png)

### Centrality and boxplot options

The diamond marks the median of the expressing cells; switch to the mean
with `centrality.fn = "mean"`, or remove it with `centrality = FALSE`. A
boxplot of the non-zero values can be overlaid with `boxplot = TRUE`:

``` r

do_BoilPlot(pbmc3k, features = "GZMA", group.by = "seurat_annotations",
            group.levels = c("NK", "CD8 T"),
            centrality.fn = "mean", boxplot = TRUE)
```

![](visualisation_files/figure-html/boil-options-1.png)

Point size, transparency and jitter width follow the original boil plot
(`pt.size = 0.7`, `pt.alpha = 0.6`, `jitter.width = 0.1`);
`violin.alpha` and `violin.width` control the violins. Only genes are
accepted — for metadata columns such as `nCount_RNA`, which are never
zero-inflated, use
[`do_ViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md).

## Sankey / alluvial diagrams

[`do_SankeyPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_SankeyPlot.md)
visualises relationships between categorical metadata variables using an
alluvial (Sankey) diagram. Each stratum is labelled and coloured; flows
connect cells across variables.

### Clusters to cell types

A common use case — mapping cluster IDs to biological annotations:

``` r

do_SankeyPlot(pbmc3k,
              columns = c("seurat_clusters", "seurat_annotations"))
```

    ## 62 cell(s) with NA values excluded from plot.

![](visualisation_files/figure-html/sankey-basic-1.png)

This immediately shows which clusters map cleanly to a single cell type
and which are split across multiple annotations.
