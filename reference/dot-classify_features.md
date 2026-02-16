# Classify Features as Genes or Metadata Columns

Internal helper that splits a character vector of feature names into
gene features (found in the assay rownames) and metadata columns (found
in object metadata). Genes take priority if a name exists in both.

## Usage

``` r
.classify_features(object, features, assay = NULL)
```

## Arguments

- object:

  Seurat object.

- features:

  character. Feature names to classify.

- assay:

  character. Assay to check (default: NULL = DefaultAssay).

## Value

A list with components: `genes`, `metadata`, `unknown`.
