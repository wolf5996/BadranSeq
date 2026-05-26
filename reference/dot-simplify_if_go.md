# Simplify GO terms if applicable

Conditionally applies
[`clusterProfiler::simplify()`](https://rdrr.io/pkg/clusterProfiler/man/simplify-methods.html)
when GO terms are detected.

## Usage

``` r
.simplify_if_go(gsea_object, cutoff = 0.7, by = "p.adjust", select_fun = min)
```

## Arguments

- gsea_object:

  A GSEA result object.

- cutoff:

  Numeric. Similarity cutoff for simplify().

- by:

  Character. Column to select representative terms.

- select_fun:

  Function to select terms.

## Value

GSEA object (possibly simplified).
