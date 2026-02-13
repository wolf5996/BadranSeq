# Changelog

## BadranSeq 1.0.0

### Visualisation

- [`do_DimPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_DimPlot.md)
  provides a unified entry point for dimensionality reduction plots,
  routing to specialised handlers based on the `reduction` argument.

- [`do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md)
  produces publication-ready UMAP plots with cell borders, boxed cluster
  labels, and a clean theme out of the box.

- [`do_PcaPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_PcaPlot.md)
  adds automatic variance-explained labels to PCA axes, formatted as
  “PC1 (X.X%)”.

- [`do_FeaturePlot()`](https://wolf5996.github.io/BadranSeq/reference/do_FeaturePlot.md)
  overlays gene expression on embeddings using a viridis colour scale
  with cells ordered by expression.

- [`EnhancedElbowPlot()`](https://wolf5996.github.io/BadranSeq/reference/EnhancedElbowPlot.md)
  extends Seurat’s elbow plot with a variance percentage axis and an
  optional visual cutoff line.

- Split-panel silhouette feature: when `split.by` is used, each panel
  shows coloured cells for the current category against a grey
  silhouette of all cells, preserving spatial context across panels.

- [`do_StatsViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_StatsViolinPlot.md)
  creates statistical violin plots via
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html)
  with Seurat-aware data extraction, automatic colour generation, and
  nonparametric testing (Kruskal–Wallis + Dunn’s) by default. The
  `group.levels` argument allows comparing a subset of groups without
  manually subsetting the Seurat object. The `pairwise.display` argument
  controls bracket visibility.

- [`do_ViolinPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_ViolinPlot.md)
  creates descriptive violin plots with boxplot overlay, jittered raw
  data points, and centrality markers (median/mean diamond with numeric
  label). Pure ggplot2, no statistical testing.

- [`do_SankeyPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_SankeyPlot.md)
  creates Sankey (alluvial) diagrams showing the relationship between
  2–3 categorical metadata variables from a Seurat object. Each stratum
  is labelled and coloured; flows connect cells across variables.
  Powered by ggalluvial.

### Data extraction

- [`fetch_cell_data()`](https://wolf5996.github.io/BadranSeq/reference/fetch_cell_data.md)
  extracts cell-level data from Seurat objects as a tidy tibble with one
  row per cell. Supports optional embeddings and metadata selection. No
  expression data is included.

- [`fetch_feature_data()`](https://wolf5996.github.io/BadranSeq/reference/fetch_feature_data.md)
  extracts feature expression data as a tidy tibble with one row per
  cell per feature. Expression values are one column per requested layer
  (features long, layers wide). Supports explicit assay and layer
  selection, optional embeddings, and optional metadata.

- Both functions replace
  [`Seurat::FetchData()`](https://satijalab.github.io/seurat-object/reference/FetchData.html)
  with a cleaner interface: `cell_id` as a column (not rownames), tibble
  output, and explicit control over what data is returned.

- [`get_pca_variance()`](https://wolf5996.github.io/BadranSeq/reference/get_pca_variance.md)
  returns a data.frame with variance explained and cumulative variance
  per principal component.

### Theming

- [`theme_badranseq()`](https://wolf5996.github.io/BadranSeq/reference/theme_badranseq.md)
  provides a consistent publication-ready ggplot2 theme with white
  background, no grids, bold titles, and legend at bottom.

- [`generate_badranseq_colors()`](https://wolf5996.github.io/BadranSeq/reference/generate_badranseq_colors.md)
  generates vivid, slightly darkened colours using the colorspace Dark 3
  qualitative palette.

### Interactive tools

- [`seurat_sleepwalk()`](https://wolf5996.github.io/BadranSeq/reference/seurat_sleepwalk.md)
  launches an interactive sleepwalk session for exploring cell-to-cell
  distances in embedding space.

- [`select_cells_interactive()`](https://wolf5996.github.io/BadranSeq/reference/select_cells_interactive.md)
  provides a Shiny app with brush selection for interactively selecting
  cells from embeddings.

### Data

- Bundled `pbmc3k` dataset: 3000 PBMCs preprocessed with SCTransform,
  PCA, UMAP, and clustered at resolution 0.5. Includes
  `seurat_annotations` with biological cell type labels.
