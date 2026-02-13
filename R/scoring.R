# R/scoring.R
# Signature scoring using MSigDB gene sets

#' Score Gene Set Signatures
#'
#' @description
#' Scores cells against MSigDB gene sets using
#' \code{\link[Seurat]{AddModuleScore}}. Fetches gene sets via
#' \code{\link[msigdbr]{msigdbr}}, computes per-cell module scores, and stores
#' them as metadata columns with clean pathway names.
#'
#' @param object A Seurat object. Gene names (rownames) must be HGNC symbols.
#'   Ensembl IDs and Entrez IDs are detected and rejected with an informative
#'   error.
#' @param species character. Species name for msigdbr (default: \code{"Homo
#'   sapiens"}). Available species:
#'   \itemize{
#'     \item \code{"Homo sapiens"} -- human
#'     \item \code{"Mus musculus"} -- mouse
#'     \item \code{"Rattus norvegicus"} -- rat
#'     \item \code{"Danio rerio"} -- zebrafish
#'     \item \code{"Drosophila melanogaster"} -- fruit fly
#'     \item \code{"Saccharomyces cerevisiae"} -- yeast
#'     \item \code{"Caenorhabditis elegans"}
#'     \item \code{"Sus scrofa"} -- pig
#'     \item \code{"Bos taurus"} -- cow
#'     \item \code{"Canis lupus familiaris"} -- dog
#'     \item \code{"Felis catus"} -- cat
#'     \item \code{"Gallus gallus"} -- chicken
#'     \item \code{"Equus caballus"} -- horse
#'     \item \code{"Pan troglodytes"} -- chimpanzee
#'     \item \code{"Macaca mulatta"} -- rhesus macaque
#'     \item \code{"Anolis carolinensis"} -- green anole
#'     \item \code{"Monodelphis domestica"} -- opossum
#'     \item \code{"Ornithorhynchus anatinus"} -- platypus
#'     \item \code{"Xenopus tropicalis"} -- tropical clawed frog
#'     \item \code{"Schizosaccharomyces pombe 972h-"}
#'   }
#' @param collection character. MSigDB collection (default: \code{"H"}).
#'   Available collections:
#'   \itemize{
#'     \item \code{"H"} -- Hallmark gene sets (50 sets). Well-defined biological
#'       states and processes with reduced redundancy.
#'     \item \code{"C1"} -- Positional gene sets (302 sets). Genes organised by
#'       chromosomal position.
#'     \item \code{"C2"} -- Curated gene sets. From publications,
#'       pathway databases, and expert knowledge. Subcollections:
#'       \itemize{
#'         \item \code{"CGP"} -- Chemical and Genetic Perturbations
#'         \item \code{"CP"} -- Canonical Pathways
#'         \item \code{"CP:BIOCARTA"} -- BioCarta Pathways
#'         \item \code{"CP:KEGG_LEGACY"} -- KEGG Legacy Pathways
#'         \item \code{"CP:KEGG_MEDICUS"} -- KEGG Medicus Pathways
#'         \item \code{"CP:PID"} -- PID Pathways
#'         \item \code{"CP:REACTOME"} -- Reactome Pathways
#'         \item \code{"CP:WIKIPATHWAYS"} -- WikiPathways
#'       }
#'     \item \code{"C3"} -- Regulatory target gene sets. Genes sharing a
#'       regulatory motif. Subcollections:
#'       \itemize{
#'         \item \code{"MIR:MIRDB"} -- miRDB microRNA targets
#'         \item \code{"MIR:MIR_LEGACY"} -- Legacy microRNA targets
#'         \item \code{"TFT:GTRD"} -- GTRD transcription factor targets
#'         \item \code{"TFT:TFT_LEGACY"} -- Legacy transcription factor targets
#'       }
#'     \item \code{"C4"} -- Computational gene sets. Defined by mining large
#'       collections of cancer-oriented microarray data. Subcollections:
#'       \itemize{
#'         \item \code{"3CA"} -- Curated Cancer Cell Atlas
#'         \item \code{"CGN"} -- Cancer Gene Neighborhoods
#'         \item \code{"CM"} -- Cancer Modules
#'       }
#'     \item \code{"C5"} -- Ontology gene sets. Derived from controlled
#'       vocabularies. Subcollections:
#'       \itemize{
#'         \item \code{"GO:BP"} -- GO Biological Process
#'         \item \code{"GO:CC"} -- GO Cellular Component
#'         \item \code{"GO:MF"} -- GO Molecular Function
#'         \item \code{"HPO"} -- Human Phenotype Ontology
#'       }
#'     \item \code{"C6"} -- Oncogenic signature gene sets (189 sets). Defined
#'       from microarray data of cancer gene perturbations.
#'     \item \code{"C7"} -- Immunologic signature gene sets. From immunology
#'       studies. Subcollections:
#'       \itemize{
#'         \item \code{"IMMUNESIGDB"} -- ImmuneSigDB
#'         \item \code{"VAX"} -- HIPC Vaccine Response
#'       }
#'     \item \code{"C8"} -- Cell type signature gene sets (866 sets). Curated
#'       from single-cell sequencing studies.
#'   }
#' @param subcollection character. MSigDB sub-collection (default: \code{NULL}).
#'   Only applicable for collections with sub-collections (C2, C3, C4, C5, C7).
#'   See \code{collection} for all valid subcollection values.
#' @param prefix character. Prefix for metadata column names (default:
#'   \code{NULL} = automatically derived from \code{collection}, e.g.
#'   \code{"Hallmark_"} for collection \code{"H"}).
#' @param report_only logical. If \code{TRUE}, return a tibble reporting which
#'   genes from each gene set are found in the object instead of scoring. The
#'   tibble has columns \code{gene_set}, \code{gene}, and \code{found}. All
#'   original genes from the msigdbr fetch are included (not just those present
#'   in the object). Default: \code{FALSE}.
#' @param ... Additional arguments passed to
#'   \code{\link[Seurat]{AddModuleScore}}.
#'
#' @return If \code{report_only = FALSE} (default), the Seurat object with new
#'   metadata columns (one per gene set). Column names are cleaned pathway names
#'   with the specified prefix. If \code{report_only = TRUE}, a tibble with
#'   columns \code{gene_set} (character), \code{gene} (character), and
#'   \code{found} (logical).
#'
#' @examples
#' \dontrun{
#' # Score Hallmark gene sets (default)
#' obj <- score_signatures(obj)
#'
#' # Score KEGG pathways
#' obj <- score_signatures(obj, collection = "C2",
#'                         subcollection = "CP:KEGG_MEDICUS")
#'
#' # Score mouse cell type signatures
#' obj <- score_signatures(obj, species = "Mus musculus", collection = "C8")
#'
#' # Custom prefix
#' obj <- score_signatures(obj, collection = "C6", prefix = "Onco_")
#'
#' # Report-only mode: see which genes are found
#' report <- score_signatures(obj, collection = "H", report_only = TRUE)
#' }
#'
#' @export
score_signatures <- function(
    object,
    species = "Homo sapiens",
    collection = "H",
    subcollection = NULL,
    prefix = NULL,
    report_only = FALSE,
    ...
) {

  # --- 1. Validate inputs ---

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "Package 'msigdbr' is required for score_signatures(). ",
      "Install with: install.packages('msigdbr')"
    )
  }

  # --- 2. Fetch gene sets ---

  msigdbr_args <- list(species = species, collection = collection)
  if (!is.null(subcollection)) {
    msigdbr_args$subcollection <- subcollection
  }

  gene_sets_df <- tryCatch(
    do.call(msigdbr::msigdbr, msigdbr_args),
    error = function(e) {
      stop(
        "No gene sets found for species = '", species,
        "', collection = '", collection, "'",
        if (!is.null(subcollection)) paste0(", subcollection = '", subcollection, "'"),
        ". Check msigdbr::msigdbr_species() and msigdbr::msigdbr_collections(). ",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  if (nrow(gene_sets_df) == 0) {
    stop(
      "No gene sets found for species = '", species,
      "', collection = '", collection, "'",
      if (!is.null(subcollection)) paste0(", subcollection = '", subcollection, "'"),
      ". Check msigdbr::msigdbr_species() and msigdbr::msigdbr_collections()."
    )
  }

  # --- 3. Gene name identifier check ---

  sample_genes <- utils::head(rownames(object), 50)

  ensembl_frac <- mean(grepl("^ENS[A-Z]+[0-9]", sample_genes))
  if (ensembl_frac > 0.5) {
    stop(
      "Gene names appear to be Ensembl IDs (e.g. '", sample_genes[1], "'). ",
      "msigdbr returns HGNC symbols. Please convert your Seurat object ",
      "rownames from Ensembl IDs to HGNC symbols before scoring."
    )
  }

  entrez_frac <- mean(grepl("^[0-9]+$", sample_genes))
  if (entrez_frac > 0.5) {
    stop(
      "Gene names appear to be Entrez IDs (e.g. '", sample_genes[1], "'). ",
      "msigdbr returns HGNC symbols. Please convert your Seurat object ",
      "rownames from Entrez IDs to HGNC symbols before scoring."
    )
  }

  # --- 4. Split into gene set list ---

  gene_sets_all <- split(gene_sets_df$gene_symbol, gene_sets_df$gs_name)

  # --- 5. Filter to present genes ---

  object_genes <- rownames(object)
  gene_sets_filtered <- lapply(gene_sets_all, function(gs) {
    intersect(gs, object_genes)
  })

  keep <- vapply(gene_sets_filtered, length, integer(1)) >= 2L
  n_dropped <- sum(!keep)
  gene_sets_filtered <- gene_sets_filtered[keep]

  if (n_dropped > 0) {
    message(n_dropped, " gene set(s) dropped (< 2 genes found in object).")
  }

  if (length(gene_sets_filtered) == 0) {
    stop("No gene sets had >= 2 genes present in the object.")
  }

  # --- 6. Report-only mode ---

  if (report_only) {
    object_genes_set <- rownames(object)
    rows <- lapply(names(gene_sets_all), function(gs_name) {
      genes <- gene_sets_all[[gs_name]]
      tibble::tibble(
        gene_set = gs_name,
        gene     = genes,
        found    = genes %in% object_genes_set
      )
    })
    report <- do.call(rbind, rows)
    return(report)
  }

  # --- 7. Determine prefix ---

  if (is.null(prefix)) {
    prefix <- switch(collection,
      "H"  = "Hallmark_",
      "C1" = "Positional_",
      "C2" = "Curated_",
      "C3" = "Regulatory_",
      "C4" = "Computational_",
      "C5" = "Ontology_",
      "C6" = "Oncogenic_",
      "C7" = "Immunologic_",
      "C8" = "CellType_",
      paste0(collection, "_")
    )
  }

  message(
    "Scoring ", length(gene_sets_filtered), " gene sets ",
    "(", collection,
    if (!is.null(subcollection)) paste0(" / ", subcollection),
    ")..."
  )

  # --- 8. AddModuleScore ---

  internal_prefix <- paste0(
    ".sigscr_", format(Sys.time(), "%Y%m%d%H%M%S"), "_"
  )

  object <- Seurat::AddModuleScore(
    object,
    features = gene_sets_filtered,
    name = internal_prefix,
    ...
  )

  # --- 9. Rename columns ---

  score_cols <- paste0(internal_prefix, seq_along(gene_sets_filtered))

  pathway_names <- names(gene_sets_filtered)
  clean_names <- sub("^HALLMARK_", "", pathway_names)
  clean_names <- sub("^GO_", "", clean_names)
  clean_names <- sub("^REACTOME_", "", clean_names)
  clean_names <- sub("^KEGG_MEDICUS_", "", clean_names)
  clean_names <- sub("^KEGG_", "", clean_names)
  clean_names <- sub("^BIOCARTA_", "", clean_names)
  clean_names <- sub("^PID_", "", clean_names)
  clean_names <- sub("^WP_", "", clean_names)
  clean_names <- sub("^HP_", "", clean_names)

  final_names <- paste0(prefix, clean_names)

  idx <- match(score_cols, colnames(object@meta.data))
  colnames(object@meta.data)[idx] <- final_names

  message("Added ", length(final_names), " signature scores to metadata.")

  # --- 10. Return ---

  object
}
