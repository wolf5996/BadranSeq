# BadranSeq 1.1.0

## New features

* `do_SankeyPlot()` creates Sankey (alluvial) diagrams showing the relationship
  between 2–3 categorical metadata variables from a Seurat object. Each
  stratum is labelled and coloured; flows connect cells across variables.
  Powered by ggalluvial.

* `seurat_to_h5ad()` converts a Seurat object to a Python AnnData object
  (and optionally writes an .h5ad file). The first requested layer becomes
  AnnData.X; additional layers and all reductions are transferred. Replaces
  sceasy, which does not support Seurat v5. Requires `anndata` and
  `reticulate` (both in Suggests).

---

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
