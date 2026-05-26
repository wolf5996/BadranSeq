# GSEA Results: CD14+ Monocytes vs Others (PBMC 3K)

GSEA results from comparing Cluster 2 (\>97% CD14+ Monocytes) against
all other clusters in the `pbmc3k` dataset. Generated via
[`clusterProfiler::gseGO()`](https://rdrr.io/pkg/clusterProfiler/man/gseGO.html)
using biological process ontology on a ranked gene list from Wilcoxon
differential expression analysis.

## Usage

``` r
gsea_pbmc3k
```

## Format

A `gseaResult` object from `clusterProfiler` with slots:

- result:

  data.frame with columns: ID, Description, setSize, enrichmentScore,
  NES, pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment

## Source

Differential expression (Cluster 2 vs others) on `pbmc3k` dataset, GSEA
via `clusterProfiler::gseGO(ont = "BP")`.

## Examples

``` r
data(gsea_pbmc3k)
gsea_pbmc3k
#> #
#> # Gene Set Enrichment Analysis
#> #
#> #...@organism     Homo sapiens 
#> #...@setType      BP 
#> #...@keytype      ENTREZID 
#> #...@geneList     Named num [1:2351] 103.4 103.4 103.4 103.4 94.5 ...
#>  - attr(*, "names")= chr [1:2351] "3957" "6279" "2219" "6280" ...
#> #...nPerm     
#> #...pvalues adjusted by 'BH' with cutoff < 0.05
#> #...101 enriched terms found
#> 'data.frame':    101 obs. of  11 variables:
#>  $ ID             : chr  "GO:0060326" "GO:0030316" "GO:0030595" "GO:0006959" ...
#>  $ Description    : chr  "cell chemotaxis" "osteoclast differentiation" "leukocyte chemotaxis" "humoral immune response" ...
#>  $ setSize        : int  66 22 55 45 19 26 54 51 51 50 ...
#>  $ enrichmentScore: num  0.575 0.865 0.589 0.642 0.813 ...
#>  $ NES            : num  2.24 2.56 2.26 2.46 2.3 ...
#>  $ pvalue         : num  6.17e-08 1.26e-07 7.59e-07 1.49e-06 1.75e-05 ...
#>  $ p.adjust       : num  0.000141 0.000144 0.000578 0.000853 0.002984 ...
#>  $ qvalue         : num  0.000123 0.000125 0.000504 0.000745 0.002605 ...
#>  $ rank           : num  272 204 226 194 45 212 182 226 226 226 ...
#>  $ leading_edge   : chr  "tags=33%, list=12%, signal=30%" "tags=36%, list=9%, signal=34%" "tags=33%, list=10%, signal=30%" "tags=38%, list=8%, signal=35%" ...
#>  $ core_enrichment: chr  "S100A8/S100A9/AIF1/GSTP1/FCER1G/LGALS3/CSF3R/SPI1/RAC1/S100A12/BST1/CXCR4/CCL3/IL1B/RPS19/RPL13A/C5AR1/CD74/LGA"| __truncated__ "TYROBP/FCER1G/SPI1/CEBPB/ANXA2/CCL3/MAFB/OSCAR" "S100A8/S100A9/AIF1/FCER1G/LGALS3/CSF3R/SPI1/RAC1/S100A12/BST1/CXCR4/CCL3/IL1B/RPS19/RPL13A/C5AR1/CD74/LGALS9" "FCN1/S100A9/LYZ/CFD/LGALS3/CFP/S100A12/BST1/GAPDH/HLA-DRB1/RPL30/CCL3/IL1B/SLC11A1/RPS19/RNASE6/B2M" ...
#> #...Citation
#> S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu. Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320 
#> 
```
