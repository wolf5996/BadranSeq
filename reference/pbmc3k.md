# PBMC 3K Dataset

A subset of 3000 peripheral blood mononuclear cells (PBMCs) from a
healthy donor, preprocessed as a Seurat object for demonstration
purposes.

## Usage

``` r
pbmc3k
```

## Format

A Seurat object containing:

- assays:

  RNA assay with count and data layers

- reductions:

  PCA and UMAP dimensionality reductions

- meta.data:

  Cell metadata including cluster assignments

## Source

10x Genomics PBMC dataset, processed with standard Seurat workflow

## Examples

``` r
data(pbmc3k)
pbmc3k
#> An object of class Seurat 
#> 26286 features across 2700 samples within 2 assays 
#> Active assay: SCT (12572 features, 3000 variable features)
#>  3 layers present: counts, data, scale.data
#>  1 other assay present: RNA
#>  2 dimensional reductions calculated: pca, umap
```
