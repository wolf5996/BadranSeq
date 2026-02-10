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

| `reduction` | Routes to                                                                                             | Axis labels            |
|-------------|-------------------------------------------------------------------------------------------------------|------------------------|
| `"umap"`    | [`do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md)                      | UMAP1, UMAP2           |
| `"pca"`     | [`do_PcaPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_PcaPlot.md)                        | PC1 (X.X%), PC2 (X.X%) |
| other       | [`.do_DimPlot_internal()`](https://wolf5996.github.io/BadranSeq/reference/dot-do_DimPlot_internal.md) | REDUCTION1, REDUCTION2 |

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
p1 <- do_FeaturePlot(pbmc3k, features = "CD3D")
p2 <- do_FeaturePlot(pbmc3k, features = "CD14")
p3 <- do_FeaturePlot(pbmc3k, features = "MS4A1")

p1 + p2 + p3
```

![](visualisation_files/figure-html/featureplot-1.png)

Cells are ordered by expression by default (`order = TRUE`), so
high-expressing cells are plotted on top — avoiding the common problem
of rare populations being hidden beneath low-expressing cells.

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
