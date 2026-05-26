# Check if GSEA object contains GO terms

Detects whether the result table contains GO identifiers (GO: followed
by 7 digits, e.g. GO:0000000).

## Usage

``` r
.is_go_gsea(gsea_object)
```

## Arguments

- gsea_object:

  A GSEA result object with `@result` slot.

## Value

Logical.
