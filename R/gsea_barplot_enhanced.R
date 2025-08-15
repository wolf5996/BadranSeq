#' GSEA Barplot Enhanced
#'
#' @param gsea_object A `GSEA` object containing the results of a Gene Set Enrichment Analysis. Usually created using the `clusterProfiler`  or `ReactomePA` packages.
#' @param analysis_name A character string representing the name of the analysis for labeling purposes.
#'
#' @returns A ggplot object visualizing GSEA results.
#' Instead of the default dotplots, which do not clearly convey normalized
#' enrichment scores (NES), this function generates a barplot of the top 10
#' upregulated and downregulated pathways. The plot displays NES values along
#' with adjusted p-values (padj), and includes leading-edge genes in the labels
#' for easier interpretation and to facilitate identification of duplicate terms.
#'
#' @export
#'
#' @examples
create_enhanced_plot <- function(gsea_object, analysis_name) {

  # extract results dataframe
  gsea_df <- gsea_object@result

  if (nrow(gsea_df) == 0) {
    cat("No significant pathways found for", analysis_name, "\n")
    return(NULL)
  }

  # prepare data for visualization
  gsea_df <- gsea_df %>%
    dplyr::mutate(
      pathway = Description,
      padj = p.adjust,
      # extract leading edge genes and count them
      leading_edge_count = purrr::map_int(core_enrichment, function(x) {
        length(stringr::str_split(x, "/")[[1]])
      }),
      leading_edge_clean = purrr::map_chr(core_enrichment, function(x) {
        genes <- stringr::str_split(x, "/")[[1]]
        # take first 3-4 genes to avoid overcrowding
        if (length(genes) > 4) {
          paste(paste(genes[1:3], collapse = ", "), "...", sep = ", ")
        } else {
          paste(genes, collapse = ", ")
        }
      })
    ) %>%
    dplyr::select(pathway, NES, padj, leading_edge_clean, leading_edge_count)

  top <- gsea_df %>%
    dplyr::filter(NES > 0) %>%
    dplyr::slice_max(order_by = NES, n = 10)

  bottom <- gsea_df %>%
    dplyr::filter(NES < 0) %>%
    dplyr::slice_min(order_by = NES, n = 10)

  both <- rbind(top, bottom)

  p1 <- both %>%
    dplyr::mutate(
      pos = NES/2,
      padj_formatted = formatC(padj, 2, format = "e"),
      # create pathway labels with gene count and leading edge genes
      pathway_with_genes = paste0(pathway, "\n(n=", leading_edge_count, ": ", leading_edge_clean, ")")
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(
      pathway_with_genes = factor(
        pathway_with_genes,
        levels = pathway_with_genes
      ),
      colortext = dplyr::case_when(
        pos > 0 ~ "black",
        pos < 0 ~ "white"
      )
    ) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = pathway_with_genes, y = NES, label = padj_formatted)
    ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = NES)
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "C"
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Pathway (n=gene count: Leading Edge Genes)",
      y = "Normalized Enrichment Score",
      title = paste("GSEA Results -", analysis_name)
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      lty = 2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = pos, color = colortext),
      size = 4,
      fontface = "bold"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::theme(
      text = ggplot2::element_text(
        size = 24,
        family = "sans"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10,
        lineheight = 0.8
      ),
      plot.margin = ggplot2::margin(10, 50, 10, 10)
    )

  return(p1)
}
