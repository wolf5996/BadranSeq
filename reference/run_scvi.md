# Run scVI Batch Correction

Wraps the scVI variational autoencoder workflow into a single function
call. Takes a Seurat object, converts to AnnData via
[`seurat_to_h5ad`](https://wolf5996.github.io/BadranSeq/reference/seurat_to_h5ad.md),
trains an scVI model, and stores the latent representation back in the
Seurat object as a dimensional reduction.

Requires Python packages `scvi-tools` and `anndata`, plus R packages
`reticulate` and `anndata`.

## Usage

``` r
run_scvi(
  object,
  batch = NULL,
  assay = NULL,
  layer = "counts",
  conda_env = NULL,
  n_hvg = 2000L,
  max_epochs = 200L,
  n_latent = 30L,
  reduction.name = "scvi",
  reduction.key = "scvi_",
  seed = 42L
)
```

## Arguments

- object:

  Seurat object.

- batch:

  character. Metadata column name for batch correction (default: NULL
  for no batch correction).

- assay:

  character. Assay to use (default: NULL = DefaultAssay).

- layer:

  character. Layer containing raw counts (default: "counts"). scVI
  requires raw (unnormalized) counts.

- conda_env:

  character. Name of conda environment containing scvi-tools (default:
  NULL = use current Python).

- n_hvg:

  integer. Number of highly variable features to select (default: 2000).
  The object is normalised and variable features are identified
  internally; only the top `n_hvg` genes are passed to scVI.

- max_epochs:

  integer. Maximum training epochs (default: 200).

- n_latent:

  integer. Number of latent dimensions (default: 30).

- reduction.name:

  character. Name for the DimReduc object (default: "scvi").

- reduction.key:

  character. Key prefix for dimension names (default: "scvi\_").

- seed:

  integer. Random seed for reproducibility (default: 42).

## Value

The Seurat object with a new dimensional reduction stored under
`reduction.name`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic batch correction
obj <- run_scvi(obj, batch = "stim", conda_env = "scvi-env")

# Custom settings
obj <- run_scvi(obj, batch = "stim", n_hvg = 3000L,
                n_latent = 50L, max_epochs = 400L,
                conda_env = "scvi-env")

# Use for downstream analysis
obj <- Seurat::FindNeighbors(obj, reduction = "scvi", dims = 1:30)
obj <- Seurat::RunUMAP(obj, reduction = "scvi", dims = 1:30)
} # }
```
