# BadranSeq (development version)

## New features

* Added `fetch_cell_data()` for tidy extraction of cell-level data from Seurat
  objects. Returns a tibble with one row per cell containing cell identifiers,
  optional dimensionality reduction embeddings, and optional metadata. No
  expression data is included.

* Added `fetch_feature_data()` for tidy extraction of feature expression data
  from Seurat objects. Returns a tibble with one row per cell per feature, with
  expression values as one column per requested layer (features long, layers
  wide). Supports explicit assay and layer selection, optional embeddings, and
  optional metadata.

* Both functions replace `Seurat::FetchData()` with a cleaner interface:
  `cell_id` as a column (not rownames), tibble output, and explicit control
  over what data is returned.

## Internal

* Added `.build_cell_scaffold()` internal helper shared by both public
  functions.

## Dependencies

* Added `tibble`, `dplyr`, `tidyr`, `rlang`, `purrr`, and `SeuratObject` to
  Imports.
