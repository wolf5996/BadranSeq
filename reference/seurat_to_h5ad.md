# Convert a Seurat Object to AnnData

Converts a Seurat object to a Python AnnData object and optionally
writes it to an .h5ad file. The first element of `layers` becomes
`AnnData.X`; remaining layers are stored in `AnnData.layers`. All
reductions present in the Seurat object are transferred to
`AnnData.obsm` (e.g. PCA becomes `X_pca`, UMAP becomes `X_umap`).

Requires the anndata R package and the Python `anndata` package.

## Usage

``` r
seurat_to_h5ad(
  seurat_obj,
  outfile = NULL,
  assay = "RNA",
  layers = c("counts", "data")
)
```

## Arguments

- seurat_obj:

  Seurat object.

- outfile:

  character. Path to output .h5ad file (default: NULL = no file written,
  AnnData returned in memory only).

- assay:

  character. Assay to convert (default: "RNA").

- layers:

  character. Vector of layer names to include. The first element becomes
  AnnData.X; the rest are stored as AnnData layers (default:
  `c("counts", "data")`).

## Value

The AnnData object (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
# Write to file
seurat_to_h5ad(obj, outfile = "data.h5ad")

# Get in-memory AnnData (no file)
adata <- seurat_to_h5ad(obj, layers = "counts")

# Custom assay and layers
seurat_to_h5ad(obj, outfile = "rna.h5ad",
               assay = "RNA", layers = c("counts", "data"))
} # }
```
