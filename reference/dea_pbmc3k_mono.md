# DEA Results: CD14+ Monocytes vs Others (PBMC 3K)

Wilcoxon differential expression results comparing Cluster 2 (\>97%
CD14+ Monocytes) against all other clusters in the `pbmc3k` dataset.
Generated via
[`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
with `logfc.threshold = 0` so all genes are retained for downstream GSEA
ranking.

## Usage

``` r
dea_pbmc3k_mono
```

## Format

A data.frame with 2530 rows (all detected genes) and 5 columns:

- p_val:

  Raw p-value from Wilcoxon rank-sum test

- avg_log2FC:

  Average log2 fold-change (Mono vs Others)

- pct.1:

  Percentage of cells in Cluster 2 expressing the gene

- pct.2:

  Percentage of cells in all other clusters expressing the gene

- p_val_adj:

  BH-adjusted p-value

Row names are HGNC gene symbols.

## Source

`Seurat::FindMarkers(ident.1 = "2")` on `pbmc3k` dataset.

## Examples

``` r
data(dea_pbmc3k_mono)
head(dea_pbmc3k_mono)
#>                p_val avg_log2FC pct.1 pct.2     p_val_adj
#> S100A8  0.000000e+00  107.06441 0.970 0.120  0.000000e+00
#> LGALS2  0.000000e+00   44.27975 0.905 0.058  0.000000e+00
#> S100A9  0.000000e+00  186.82737 0.996 0.213  0.000000e+00
#> FCN1   1.196034e-319   31.00664 0.955 0.146 1.640241e-315
#> CD14   2.584069e-299   16.66226 0.665 0.027 3.543793e-295
#> MS4A6A 1.714710e-276   42.55068 0.682 0.041 2.351553e-272
```
