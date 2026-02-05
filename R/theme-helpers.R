# BadranSeq Theme Helpers
# Internal functions for consistent publication-ready styling

# ============================================================================
# Color Palette
# ============================================================================

#' Default BadranSeq Color Palette
#'
#' @description
#' Predefined colors inspired by colorspace "Dark 3" palette.
#' Extended to 20 distinct colors for larger datasets.
#'
#' @keywords internal
.badranseq_palette <- c(

"#E16A86", "#909800", "#00AD9A", "#9183E6",
  "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#666666", "#1B9E77",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"
)

#' Generate Categorical Color Palette
#'
#' @description
#' Internal function to generate publication-ready categorical colors
#' similar to SCpubr's default palette.
#'
#' @param n numeric. Number of colors needed.
#' @param custom.colors character. Optional custom color vector.
#'
#' @return Character vector of hex colors.
#' @keywords internal
generate_badranseq_colors <- function(n, custom.colors = NULL) {

  if (!is.null(custom.colors)) {
    if (length(custom.colors) < n) {
      warning("Not enough custom colors provided. Recycling colors.")
      custom.colors <- rep_len(custom.colors, n)
    }
    return(custom.colors[1:n])
  }

  if (n <= length(.badranseq_palette)) {
    return(.badranseq_palette[1:n])
  } else {
    # Generate more colors using scales
    return(scales::hue_pal()(n))
  }
}

# ============================================================================
# Theme
# ============================================================================

#' BadranSeq Publication-Ready Theme
#'
#' @description
#' Internal theme function providing consistent SCpubr-like aesthetics
#' for publication-ready plots.
#'
#' @param base_size numeric. Base font size (default: 14).
#' @param plot.axes logical. Whether to show axes (default: TRUE).
#'
#' @return ggplot2 theme object.
#' @keywords internal
theme_badranseq <- function(base_size = 14, plot.axes = TRUE) {

  # Base theme
theme_base <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # Background - white throughout
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),

      # Remove panel grid
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),

      # Title styling
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(face = "plain", hjust = 0),
      plot.caption = ggplot2::element_text(face = "italic", hjust = 1),
      plot.title.position = "plot",

      # Legend positioning (bottom, title on top)
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.text = ggplot2::element_text(face = "plain"),

      # Minimal margins
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    )

  # Conditional axis styling
  if (plot.axes) {
    theme_base <- theme_base +
      ggplot2::theme(
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(face = "plain", color = "black"),
        axis.ticks = ggplot2::element_line(color = "black"),
        axis.line = ggplot2::element_line(color = "black", linewidth = 0.5)
      )
  } else {
    theme_base <- theme_base +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank()
      )
  }

  return(theme_base)
}
