# BadranSeq Theme Helpers
# Internal functions for consistent publication-ready styling

# ============================================================================
# Color Palette
# ============================================================================

#' Generate Categorical Color Palette
#'
#' @description
#' Internal function to generate publication-ready categorical colors
#' matching SCpubr's default palette. Uses colorspace "Dark 3" palette
#' with adjusted saturation and value for more vivid, slightly darker colors.
#'
#' @param n numeric. Number of colors needed.
#' @param custom.colors character. Optional custom color vector.
#'
#' @return Character vector of hex colors.
#' @keywords internal
#' @export
generate_badranseq_colors <- function(n, custom.colors = NULL) {

  if (!is.null(custom.colors)) {
    if (length(custom.colors) < n) {
      warning("Not enough custom colors provided. Recycling colors.")
      custom.colors <- rep_len(custom.colors, n)
    }
    return(custom.colors[1:n])
  }

  # Generate colors using colorspace Dark 3 palette (same as SCpubr)
  colors <- colorspace::qualitative_hcl(n, palette = "Dark 3")

  # Convert to RGB then HSV for adjustment
  colors_rgb <- grDevices::col2rgb(colors)
  colors_hsv <- grDevices::rgb2hsv(colors_rgb)

  # Adjust saturation and value (SCpubr approach)
  # Decrease value by 0.1 (slightly darker)
  colors_hsv["v", ] <- colors_hsv["v", ] - 0.1
  # Increase saturation by 0.2 (more vivid), cap at 1
 colors_hsv["s", ] <- colors_hsv["s", ] + 0.2
  colors_hsv["s", ][colors_hsv["s", ] > 1] <- 1

  # Convert back to hex
  colors <- grDevices::hsv(
    h = colors_hsv["h", ],
    s = colors_hsv["s", ],
    v = colors_hsv["v", ],
    alpha = 1
  )

  return(colors)
}

# ============================================================================
# Theme
# ============================================================================

#' BadranSeq Publication-Ready Theme
#'
#' @description
#' Theme function providing consistent SCpubr-like aesthetics
#' for publication-ready plots. Can be applied to any ggplot2 object.
#'
#' @param base_size numeric. Base font size (default: 14).
#' @param plot.axes logical. Whether to show axes (default: TRUE).
#'
#' @return ggplot2 theme object.
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_badranseq()
#' }
#'
#' @export
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
