# SleepWalk Wrapper for Seurat Objects

**\[stable\]**

## Usage

``` r
seurat_sleepwalk(object, embedding = "umap", features = "pca")
```

## Arguments

- object:

  Seurat object.

- embedding:

  character. Reduction for 2D embedding (default: "umap").

- features:

  character. Reduction for feature space (default: "pca").

## Value

sleepwalk interactive plot.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_sleepwalk(seurat_obj, embedding = "umap", features = "pca")
} # }
```
