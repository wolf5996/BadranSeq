# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Build & Development Commands

``` bash
# Generate documentation (roxygen2)
Rscript -e "devtools::document()"

# Load package for interactive development
Rscript -e "devtools::load_all()"

# Run tests (45 tests across 2 test files)
Rscript -e "devtools::test()"

# Check package (comprehensive) — expect 0 errors, 0 warnings, 0 notes
Rscript -e "devtools::check()"

# Build and install locally
R CMD build .
R CMD INSTALL BadranSeq_*.tar.gz

# Build pkgdown site (renders articles, requires package installed)
Rscript -e "devtools::install(quick = TRUE); pkgdown::build_site()"

# Build pkgdown articles only
Rscript -e "pkgdown::build_articles()"
```

## Package Purpose

BadranSeq (v1.0.0) addresses common pain points in single-cell RNA-seq
analysis workflows:

1.  **Repetitive boilerplate** — Standard Seurat plots require repeated
    customization for publication quality
2.  **Missing information** — PCA plots typically lack variance
    explained percentages on axes
3.  **Split comparisons** — Silhouette feature shows split populations
    with grey context of all cells
4.  **No tidy extraction** — Seurat’s `FetchData()` returns data.frames
    with rownames, no layer control, no tidy format
5.  **Manual cell selection** — No easy interactive way to select cell
    subsets from embeddings

## Architecture

### Source Files

- `R/BadranSeq.R` — Core visualization functions (DimPlot, PCA, UMAP,
  FeaturePlot, ElbowPlot)
- `R/fetch.R` — Tidy data extraction functions (fetch_cell_data,
  fetch_feature_data, .build_cell_scaffold)
- `R/theme-helpers.R` — Shared theming infrastructure (colors, theme)
- `R/utils-pipe.R` — Magrittr pipe export
- `R/data.R` — Data documentation for bundled pbmc3k

### Test Files

- `tests/testthat/test-basic.R` — 3 basic package tests
- `tests/testthat/test-fetch.R` — 24 tests for scaffold +
  fetch_cell_data + fetch_feature_data (includes 1 expected warning)

### pkgdown Articles

Located in `vignettes/articles/` as `.Rmd` files (NOT formal vignettes —
pkgdown-only):

- `data-extraction.Rmd` — Covers fetch_cell_data, fetch_feature_data,
  get_pca_variance, EnhancedElbowPlot with FetchData comparisons
- `visualisation.Rmd` — Covers do_UmapPlot, do_PcaPlot, do_DimPlot,
  do_FeaturePlot, split.by silhouettes with Seurat comparisons

These are excluded from `R CMD build` via `.Rbuildignore` and rendered
only by
[`pkgdown::build_articles()`](https://pkgdown.r-lib.org/reference/build_articles.html).
They use `#|` hashpipe chunk options with `fig.width`/`fig.height`
(knitr dot syntax, not Quarto hyphen syntax).

### JOSS Paper

Located in `paper/` (excluded from package build):

- `paper/paper.qmd` — JOSS paper draft
- `paper/badranseq_joss_figures/` — Figure generation scripts and output

### Parameter Naming Convention

All functions use consistent parameter names: - `object` — Seurat object
(not `sample` or `seurat_obj`) - `reduction` — Dimensionality reduction
name - `group.by` — Metadata column for grouping - `split.by` — Metadata
column for split panels with silhouettes

### Theming Infrastructure (theme-helpers.R)

All plots use a consistent publication-ready theme matching SCpubr
aesthetics:

- `generate_badranseq_colors(n)` — Generates colors using
  `colorspace::qualitative_hcl(n, palette = "Dark 3")` with HSV
  adjustments (saturation +0.2, value -0.1) for more vivid, slightly
  darker colors
- `theme_badranseq(base_size, plot.axes)` — Consistent theme: white
  background, no grids, bold titles, legend at bottom

### Single-Cell Visualization (BadranSeq.R)

**All plotting is native ggplot2** — matches SCpubr aesthetics without
the dependency.

**Key defaults:** - `label = TRUE` — Show cluster labels by default -
`plot.axes = TRUE` — Show axes by default - `plot_cell_borders = TRUE` —
Black borders around cells - `border.size = 2` — Border size multiplier

**Internal engine**: -
[`.do_DimPlot_internal()`](https://wolf5996.github.io/BadranSeq/reference/dot-do_DimPlot_internal.md)
— Workhorse function for standard plots -
[`.do_DimPlot_split_silhouette()`](https://wolf5996.github.io/BadranSeq/reference/dot-do_DimPlot_split_silhouette.md)
— Handles split.by with silhouette approach

**Split.by silhouette feature:** When `split.by` is specified, each
panel shows: 1. Grey silhouette of ALL cells (background context) 2.
Black borders for ALL cells 3. Colored cells only for the current split
category This provides spatial context when comparing populations.

**Function hierarchy**:

    do_DimPlot(reduction="umap"|"pca"|other)
        ├── reduction="pca"  → do_PcaPlot() → native ggplot2 + variance labels
        ├── reduction="umap" → do_UmapPlot() → native ggplot2 + standard labels
        └── other            → .do_DimPlot_internal() + uppercase labels

    When split.by is used:
        └── .do_DimPlot_split_silhouette() → patchwork combined panels

**Visualization functions:** -
[`do_DimPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_DimPlot.md)
— Unified entry point; routes to specialized handlers -
[`do_PcaPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_PcaPlot.md)
— Calculates variance %, formats as “PC1 (X.X%)” -
[`do_UmapPlot()`](https://wolf5996.github.io/BadranSeq/reference/do_UmapPlot.md)
— UMAP with labels and axes by default -
[`do_FeaturePlot()`](https://wolf5996.github.io/BadranSeq/reference/do_FeaturePlot.md)
— Gene expression overlay with viridis palette -
[`get_pca_variance()`](https://wolf5996.github.io/BadranSeq/reference/get_pca_variance.md)
— Returns variance explained data.frame -
[`EnhancedElbowPlot()`](https://wolf5996.github.io/BadranSeq/reference/EnhancedElbowPlot.md)
— Variance explained with optional cutoff

### Tidy Data Extraction (fetch.R)

**Replaces
[`Seurat::FetchData()`](https://satijalab.github.io/seurat-object/reference/FetchData.html)**
with explicit, tidy alternatives.

- `.build_cell_scaffold(object, reductions, dims, metadata)` — Internal
  helper building base tibble with cell_id, optional embeddings (prefix:
  PCA→PC, UMAP→UMAP), optional metadata
- `fetch_cell_data(object, reductions, dims, metadata)` — Exported thin
  wrapper. One row per cell, no expression data
- `fetch_feature_data(object, features, layer, assay, reductions, dims, metadata)`
  — One row per cell per feature. Features long, layers wide. Uses
  [`SeuratObject::Layers()`](https://satijalab.github.io/seurat-object/reference/Layers.html)
  and
  [`SeuratObject::LayerData()`](https://satijalab.github.io/seurat-object/reference/Layers.html)
  (not `Seurat::`)

**Key implementation details:** -
[`SeuratObject::LayerData()`](https://satijalab.github.io/seurat-object/reference/Layers.html)
returns `dgCMatrix` (sparse) — must call
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) before
transposition - `rlang::%||%` used for null-coalescing (imported via
`@importFrom rlang %||%`) - `utils::globalVariables(c("feature"))`
suppresses R CMD check NSE notes

### Interactive Tools

- [`seurat_sleepwalk()`](https://wolf5996.github.io/BadranSeq/reference/seurat_sleepwalk.md)
  — Interactive distance exploration via sleepwalk
- [`select_cells_interactive()`](https://wolf5996.github.io/BadranSeq/reference/select_cells_interactive.md)
  — Shiny app with brush selection; supports additive selection,
  deselection, and clear; returns cell names or subsetted Seurat object

## Bundled Data

Located in `data/`, created by script in `data-raw/`:

- `pbmc3k` — PBMC 3k dataset (~80MB). SCTransform normalized, PCA (50
  dims), UMAP (2 dims), clustered at resolution 0.5. 10 clusters.
  Metadata: `orig.ident`, `nCount_RNA`, `nFeature_RNA`,
  `seurat_annotations`, `nCount_SCT`, `nFeature_SCT`, `SCT_snn_res.0.5`,
  `seurat_clusters`. Assays: RNA (layers: counts, data), SCT.

## Key Dependencies

**Imports:** - **Seurat** (\>= 5.0.0) — scRNA-seq data structures and
methods - **SeuratObject** — `Layers()`, `LayerData()` for fetch
functions - **ggplot2** — all plotting (native implementation, no
SCpubr) - **dplyr**, **tibble**, **tidyr**, **purrr** — tidy data
extraction pipeline - **rlang** — `%||%` operator - **colorspace** —
color palette generation (Dark 3 qualitative palette) - **sleepwalk** —
interactive embedding exploration - **shiny** — interactive cell
selection UI - **magrittr** — pipe operator - **stats** — aggregate for
label positioning

**Suggests:** - **patchwork** — combining split.by panels - **ggrastr**
— rasterization for large datasets - **testthat** (\>= 3.0.0) — testing
framework

## pkgdown Configuration

`_pkgdown.yml` defines the site at
<https://wolf5996.github.io/BadranSeq/>. Articles in
`vignettes/articles/` require backtick-escaping in the YAML `contents`
field because pkgdown evaluates names as R expressions (hyphens become
minus operators):

``` yaml
articles:
  - title: Getting Started
    contents:
      - "`articles/data-extraction`"
      - "`articles/visualisation`"
```

## Backup Files

- `BadranSeq.R.backup` — Original SCpubr-dependent implementation (kept
  for reference)
