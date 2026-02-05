# Get PCA Variance Explained Summary

Helper function to extract variance explained information from a Seurat
object.

## Usage

``` r
get_pca_variance(object, n_pcs = NULL)
```

## Arguments

- object:

  Seurat object with PCA reduction.

- n_pcs:

  numeric. Number of top PCs to return. If NULL, returns all PCs.

## Value

data.frame with PC number, variance explained percentage, and cumulative
variance.

## Examples

``` r
if (FALSE) { # \dontrun{
get_pca_variance(seurat_object, n_pcs = 10)
get_pca_variance(seurat_object)
} # }
```
