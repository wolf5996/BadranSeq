# Build Cell Scaffold Tibble

Internal helper that constructs a base tibble with cell_id as a column,
optional embedding coordinates, and optional metadata from a Seurat
object.

## Usage

``` r
.build_cell_scaffold(object, reductions = NULL, dims = 1:2, metadata = TRUE)
```

## Arguments

- object:

  Seurat object.

- reductions:

  character. Reduction names to include (NULL = none).

- dims:

  integer. Dimensions to extract per reduction (default: 1:2).

- metadata:

  logical or character. TRUE = all metadata, FALSE = none, character
  vector = specific columns.

## Value

tibble with cell_id and requested columns.
