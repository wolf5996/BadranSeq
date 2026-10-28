# BadranSeq 1.5.1

## Bug fixes

* `R CMD check` on R-devel no longer fails: parallel testthat execution is
  disabled after worker results exceeded R-devel's string length limit
  (`Error in gsub("\\s+", "", x) : result string is too long`).
* `DESCRIPTION` now declares `R (>= 4.1)` because the package uses the native
  pipe, and `LICENSE` uses the DCF stub expected by `R CMD check`.
* `NES` added to `globalVariables()` to silence a `do_GseaPlot()` NOTE.

## Documentation

* README gains a `do_BoilPlot()` example figure.

# BadranSeq 1.5.0

## New features

* `do_BoilPlot()` creates boil plots for zero-inflated gene expression. A
  detection bar below zero shows the percentage of cells expressing each gene
  per group; the violin, jittered points, median/mean diamond and optional
  boxplot (`boxplot = TRUE`) above zero use non-zero cells only. Geometry and
  colours follow the original boil plot, drawn with BadranSeq's palette and
  theme, with `group.levels` to subset and order groups and `min.nonzero` to
  drop groups with too few expressing cells. Inspired by the boil plot from boilplot-r by Lara Jerman
  (<https://github.com/laraj1/boilplot-r>).

## Minor improvements

* `do_ViolinPlot()` and `do_StatsViolinPlot()` gain an `nrow` argument for
  patchwork layout control, matching `ncol`.

# BadranSeq 1.4.0

## New features

* `do_GseaPlot()` creates publication-ready horizontal barplots
  significant improvements. Creates publication-ready horizontal barplots
  of GSEA results (clusterProfiler / ReactomePA) with:
  - Diverging RdBu color scale centered at NES = 0 (emphasizes directionality)
  - `theme_badranseq()` styling for consistency with rest of package
  - Leading-edge gene count and gene symbols in pathway labels
  - Adjustable p-value labels centered on bars
  - `simplify_go` auto-detection and simplification of GO terms
  - `selection` parameter for plotting specific pathways by name
  - `fill_scale = "viridis"` option for original sequential style

* `gsea_pbmc3k` bundled dataset: CD14+ Mono (Cluster 2) vs others GSEA results
  from `clusterProfiler::gseGO(ont = "BP")` on the `pbmc3k` Seurat object.
  Core enrichment genes stored as readable HGNC symbols via `setReadable()`.
  Created via reproducible `data-raw/gsea_pea_pbmc3k.R` script.

* `dea_pbmc3k_mono` bundled dataset: Wilcoxon differential expression results
  (2,530 genes, Cluster 2 vs others) with `avg_log2FC` column for GSEA ranking.
  Row names are HGNC gene symbols.

## Dependencies

* Added `clusterProfiler` and `org.Hs.eg.db` to `Suggests` for
  GO simplification and example data generation.

# BadranSeq 1.0.0

## Visualisation

* `do_DimPlot()` provides a unified entry point for dimensionality reduction
  plots, routing to specialised handlers based on the `reduction` argument.

* `do_UmapPlot()` produces publication-ready UMAP plots with cell borders,
  boxed cluster labels, and a clean theme out of the box.

* `do_PcaPlot()` adds automatic variance-explained labels to PCA axes,
  formatted as "PC1 (X.X%)".

* `do_FeaturePlot()` overlays gene expression on embeddings using a viridis
  colour scale with cells ordered by expression.

* `EnhancedElbowPlot()` extends Seurat's elbow plot with a variance percentage
  axis and an optional visual cutoff line.

* Split-panel silhouette feature: when `split.by` is used, each panel shows
  coloured cells for the current category against a grey silhouette of all
  cells, preserving spatial context across panels.

* `do_StatsViolinPlot()` creates statistical violin plots via
  `ggstatsplot::ggbetweenstats()` with Seurat-aware data extraction, automatic
  colour generation, and nonparametric testing (Kruskal–Wallis + Dunn's) by
  default. The `group.levels` argument allows comparing a subset of groups
  without manually subsetting the Seurat object. The `pairwise.display`
  argument controls bracket visibility.

* `do_ViolinPlot()` creates descriptive violin plots with boxplot overlay,
  jittered raw data points, and centrality markers (median/mean diamond with
  numeric label). Pure ggplot2, no statistical testing.

* `do_SankeyPlot()` creates Sankey (alluvial) diagrams showing the relationship
  between 2–3 categorical metadata variables from a Seurat object. Each
  stratum is labelled and coloured; flows connect cells across variables.
  Powered by ggalluvial.

## Data extraction

* `fetch_cell_data()` extracts cell-level data from Seurat objects as a tidy
  tibble with one row per cell. Supports optional embeddings and metadata
  selection. No expression data is included.

* `fetch_feature_data()` extracts feature expression data as a tidy tibble
  with one row per cell per feature. Expression values are one column per
  requested layer (features long, layers wide). Supports explicit assay and
  layer selection, optional embeddings, and optional metadata.

* Both functions replace `Seurat::FetchData()` with a cleaner interface:
  `cell_id` as a column (not rownames), tibble output, and explicit control
  over what data is returned.

* `get_pca_variance()` returns a data.frame with variance explained and
  cumulative variance per principal component.

## Theming

* `theme_badranseq()` provides a consistent publication-ready ggplot2 theme
  with white background, no grids, bold titles, and legend at bottom.

* `generate_badranseq_colors()` generates vivid, slightly darkened colours
  using the colorspace Dark 3 qualitative palette.

## Interactive tools

* `seurat_sleepwalk()` launches an interactive sleepwalk session for exploring
  cell-to-cell distances in embedding space.

* `select_cells_interactive()` provides a Shiny app with brush selection for
  interactively selecting cells from embeddings.

## Data

* Bundled `pbmc3k` dataset: 3000 PBMCs preprocessed with SCTransform, PCA,
  UMAP, and clustered at resolution 0.5. Includes `seurat_annotations` with
  biological cell type labels.
