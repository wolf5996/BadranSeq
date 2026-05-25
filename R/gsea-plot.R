# BadranSeq - GSEA Results Plotting
# Publication-ready barplots for Gene Set Enrichment Analysis results
# Adapted from create_enhanced_plot() in the original .BadranSeq prototype

# ============================================================================
# Internal Helpers
# ============================================================================

#' Check if GSEA object contains GO terms
#'
  #' Detects whether the result table contains GO identifiers
  #' (GO: followed by 7 digits, e.g. GO:0000000).
#'
#' @param gsea_object A GSEA result object with `@result` slot.
#' @return Logical.
#' @keywords internal
.is_go_gsea <- function(gsea_object) {
  if (is.null(gsea_object) || nrow(gsea_object@result) == 0) {
    return(FALSE)
  }
  go_ids <- gsea_object@result$ID
  any(grepl("^GO:\\d{7}$", go_ids))
}

#' Simplify GO terms if applicable
#'
#' Conditionally applies `clusterProfiler::simplify()` when GO terms are detected.
#'
#' @param gsea_object A GSEA result object.
#' @param cutoff Numeric. Similarity cutoff for simplify().
#' @param by Character. Column to select representative terms.
#' @param select_fun Function to select terms.
#' @return GSEA object (possibly simplified).
#' @keywords internal
.simplify_if_go <- function(gsea_object, cutoff = 0.7, by = "p.adjust",
                               select_fun = min) {
  if (.is_go_gsea(gsea_object)) {
    message("GO terms detected. Applying clusterProfiler::simplify()...")
    return(clusterProfiler::simplify(
      gsea_object,
      cutoff = cutoff,
      by = by,
      select_fun = select_fun
    ))
  } else {
    message("Non-GO terms detected. Skipping simplification.")
    return(gsea_object)
  }
}

# ============================================================================
# do_GseaPlot
# ============================================================================

#' Publication-Ready GSEA Results Barplot
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Creates a horizontal barplot of top pathways from a Gene Set Enrichment
#' Analysis (GSEA) result. Instead of the default dotplot, this function
#' shows Normalized Enrichment Scores (NES) as bars with adjusted p-values
#' as text labels and leading-edge gene annotations.
#'
#' @param gsea_object A `gseaResult` object from `clusterProfiler::gseGO()`,
#'   `ReactomePA::gsePathway()`, or similar, containing results in the
#'   `@result` slot.
#' @param analysis_name Character. Subtitle for the plot. Defaults to
#'   ``"GSEA Analysis"``.
#' @param simplify_go Logical. Whether to apply `clusterProfiler::simplify()`
#'   to GO terms (default: `TRUE`). Only applies when GO term IDs are
#'   detected; silently skipped otherwise.
#' @param number_of_pathways Integer. Number of top upregulated and
#'   downregulated pathways to display (default: `10`). Ignored when
#'   `selection` is provided.
#' @param selection Character vector of pathway descriptions to plot.
#'   Overrides automatic top/bottom selection. Defaults to `NULL`.
#' @param fill_scale Character. Color scale: `"diverging"` (RdBu, 0-centered)
#'   or `"viridis"` (sequential viridis-C, original style). Defaults to
#'   `"diverging"` to emphasize the signed nature of NES.
#' @param show_padj Logical. Show adjusted p-value labels on bars.
#'   Defaults to `TRUE`.
#' @param show_leading_edge Logical. Show leading-edge gene count and
#'   gene symbols in y-axis labels. Defaults to `TRUE`.
#' @param title Character. Custom plot title. If `NULL`, auto-generated
#'   from `analysis_name`.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' data(gsea_pbmc3k_cd8t)
#' do_GseaPlot(gsea_pbmc3k_cd8t, analysis_name = "CD8 T vs Others")
#'
#' # Use viridis (original style)
#' do_GseaPlot(gsea_pbmc3k_cd8t, fill_scale = "viridis")
#'
#' # Plot specific pathways
#' do_GseaPlot(
#'   gsea_pbmc3k_cd8t,
#'   selection = c("defense response to Gram-negative bacterium")
#' )
#' }
#'
#' @export
do_GseaPlot <- function(
    gsea_object,
    analysis_name = "GSEA Analysis",
    simplify_go = TRUE,
    number_of_pathways = 10,
    selection = NULL,
    fill_scale = c("diverging", "viridis"),
    show_padj = TRUE,
    show_leading_edge = TRUE,
    title = NULL
) {

  fill_scale <- match.arg(fill_scale)

  # Simplify GO terms if requested and applicable
  if (simplify_go) {
    gsea_object <- .simplify_if_go(gsea_object)
  }

  # Extract results data frame
  gsea_df <- gsea_object@result

  if (nrow(gsea_df) == 0) {
    message("No significant pathways found for ", analysis_name)
    return(NULL)
  }

  # Prepare data for visualization
  gsea_df <- gsea_df |>
    dplyr::mutate(
      pathway = Description,
      padj = p.adjust,
      # Count leading-edge genes
      leading_edge_count = purrr::map_int(
        core_enrichment,
        function(x) length(base::strsplit(x, "/")[[1]])
      ),
      leading_edge_clean = purrr::map_chr(core_enrichment, function(x) {
        genes <- base::strsplit(x, "/")[[1]]
        if (length(genes) > 4) {
          paste(paste(genes[1:3], collapse = ", "), "...", sep = ", ")
        } else {
          paste(genes, collapse = ", ")
        }
      })
    ) |>
    dplyr::select(
      pathway, NES, padj, leading_edge_clean, leading_edge_count
    )

  # Pathway selection logic
  if (!is.null(selection)) {
    both <- gsea_df |>
      dplyr::filter(pathway %in% selection)

    if (nrow(both) == 0) {
      message("None of the selected pathways found in results for ", analysis_name)
      return(NULL)
    }

    missing_pathways <- setdiff(selection, gsea_df$pathway)
    if (length(missing_pathways) > 0) {
      warning(
        "The following selected pathways were not found:\n",
        paste(missing_pathways, collapse = "\n"),
        call. = FALSE
      )
    }
  } else {
    top <- gsea_df |>
      dplyr::filter(NES > 0) |>
      dplyr::slice_max(order_by = NES, n = number_of_pathways)

    bottom <- gsea_df |>
      dplyr::filter(NES < 0) |>
      dplyr::slice_min(order_by = NES, n = number_of_pathways)

    both <- rbind(top, bottom)
  }

  # Build y-axis label
  if (show_leading_edge) {
    both$pathway_with_genes <- paste0(
      both$pathway,
      "\n(n=", both$leading_edge_count, ": ", both$leading_edge_clean, ")"
    )
  } else {
    both$pathway_with_genes <- both$pathway
  }

  # Compute label position and color
  both <- both |>
    dplyr::mutate(
      pos = NES / 2,
      padj_formatted = formatC(padj, 2, format = "e"),
      pathway_with_genes = factor(
        pathway_with_genes,
        levels = pathway_with_genes[order(NES)]
      ),
      colortext = dplyr::if_else(pos > 0, "black", "white")
    )

  # Determine title
  if (is.null(title)) {
    title <- paste("GSEA Results -", analysis_name)
  }

  # Build plot
  p <- ggplot2::ggplot(
    both,
    ggplot2::aes(x = pathway_with_genes, y = NES, label = padj_formatted)
  ) +
    ggplot2::geom_col(ggplot2::aes(fill = NES)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = ifelse(show_leading_edge,
                   "Pathway (n=gene count: Leading Edge Genes)",
                   "Pathway"),
      y = "Normalized Enrichment Score",
      title = title
    ) +
    theme_badranseq()

  # Color scale
  if (fill_scale == "diverging") {
    p <- p +
      ggplot2::scale_fill_distiller(
        palette = "RdBu",
        direction = -1,
        name = "NES"
      )
  } else {
    p <- p +
      ggplot2::scale_fill_viridis_c(option = "C", name = "NES")
  }

  # Add padj labels
  if (show_padj) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(y = pos, color = colortext),
        size = 4, fontface = "bold"
      ) +
      ggplot2::scale_color_identity()
  }

  return(p)
}
