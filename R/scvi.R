# R/scvi.R
# scVI batch correction wrapper via reticulate

#' Run scVI Batch Correction
#'
#' @description
#' Wraps the scVI variational autoencoder workflow into a single function call.
#' Takes a Seurat object, converts to AnnData via \code{\link{seurat_to_h5ad}},
#' trains an scVI model, and stores the latent representation back in the Seurat
#' object as a dimensional reduction.
#'
#' Requires Python packages \code{scvi-tools} and \code{anndata}, plus R
#' packages \code{reticulate} and \code{anndata}.
#'
#' @param object Seurat object.
#' @param batch character. Metadata column name for batch correction
#'   (default: NULL for no batch correction).
#' @param assay character. Assay to use (default: NULL = DefaultAssay).
#' @param layer character. Layer containing raw counts (default: "counts").
#'   scVI requires raw (unnormalized) counts.
#' @param conda_env character. Name of conda environment containing scvi-tools
#'   (default: NULL = use current Python).
#' @param n_hvg integer. Number of highly variable features to select
#'   (default: 2000). The object is normalised and variable features are
#'   identified internally; only the top \code{n_hvg} genes are passed to scVI.
#' @param max_epochs integer. Maximum training epochs (default: 200).
#' @param n_latent integer. Number of latent dimensions (default: 30).
#' @param reduction.name character. Name for the DimReduc object
#'   (default: "scvi").
#' @param reduction.key character. Key prefix for dimension names
#'   (default: "scvi_").
#' @param seed integer. Random seed for reproducibility (default: 42).
#'
#' @return The Seurat object with a new dimensional reduction stored under
#'   \code{reduction.name}.
#'
#' @examples
#' \dontrun{
#' # Basic batch correction
#' obj <- run_scvi(obj, batch = "stim", conda_env = "scvi-env")
#'
#' # Custom settings
#' obj <- run_scvi(obj, batch = "stim", n_hvg = 3000L,
#'                 n_latent = 50L, max_epochs = 400L,
#'                 conda_env = "scvi-env")
#'
#' # Use for downstream analysis
#' obj <- Seurat::FindNeighbors(obj, reduction = "scvi", dims = 1:30)
#' obj <- Seurat::RunUMAP(obj, reduction = "scvi", dims = 1:30)
#' }
#'
#' @export
run_scvi <- function(
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
) {

  # --- Validate Seurat object ---

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  # --- Check R dependencies ---

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required for run_scvi(). ",
      "Install it with: install.packages('reticulate')"
    )
  }

  # --- Validate assay and layer ---

  assay <- assay %||% SeuratObject::DefaultAssay(object)

  if (!assay %in% SeuratObject::Assays(object)) {
    stop(
      "Assay '", assay, "' not found. Available: ",
      paste(SeuratObject::Assays(object), collapse = ", ")
    )
  }

  available_layers <- SeuratObject::Layers(object[[assay]])
  if (!layer %in% available_layers) {
    stop(
      "Layer '", layer, "' not found in assay '", assay,
      "'. Available: ", paste(available_layers, collapse = ", ")
    )
  }

  # --- Validate batch column ---

  if (!is.null(batch)) {
    if (!batch %in% colnames(object@meta.data)) {
      stop(
        "Batch column '", batch, "' not found in metadata. Available: ",
        paste(colnames(object@meta.data), collapse = ", ")
      )
    }
  }

  # --- Set up Python environment ---

  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  }

  # --- Select highly variable features ---

  message("Selecting top ", n_hvg, " variable features...")

  work_obj <- Seurat::NormalizeData(object, verbose = FALSE)
  work_obj <- Seurat::FindVariableFeatures(
    work_obj, nfeatures = n_hvg, verbose = FALSE
  )
  hvg <- SeuratObject::VariableFeatures(work_obj)
  work_obj <- work_obj[hvg, ]

  message(length(hvg), " variable features selected.")

  # --- Convert Seurat to AnnData ---

  adata <- seurat_to_h5ad(work_obj, assay = assay, layers = layer)

  # --- Set up and train scVI ---

  scvi <- reticulate::import("scvi")
  scvi$settings$seed <- as.integer(seed)

  message("Setting up scVI model...")
  scvi$model$SCVI$setup_anndata(adata, batch_key = batch)

  model <- scvi$model$SCVI(adata, n_latent = as.integer(n_latent))

  message("Training scVI model (max_epochs: ", max_epochs, ")...")
  model$train(max_epochs = as.integer(max_epochs))

  # --- Extract latent representation ---

  message("Extracting latent representation...")
  latent <- model$get_latent_representation()

  latent_mat <- as.matrix(latent)
  rownames(latent_mat) <- colnames(object)
  colnames(latent_mat) <- paste0(reduction.key, seq_len(ncol(latent_mat)))

  # --- Store in Seurat ---

  object[[reduction.name]] <- SeuratObject::CreateDimReducObject(
    embeddings = latent_mat,
    key = reduction.key,
    assay = assay
  )

  message(
    "Done. Latent representation stored in object[['",
    reduction.name, "']]"
  )

  object
}
