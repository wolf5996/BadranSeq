#' Interactive Cell Selection for Seurat Objects
#'
#' @param seurat_obj A Seurat object with dimensionality reduction computed
#' @param reduction Which reduction to use for plotting. Default is "umap". Options: "umap", "tsne", "pca"
#' @param color_by Column name from metadata to color cells by. Default is `Idents()`.
#' @param return_cells Logical, if TRUE returns selected cell names, if FALSE returns subsetted Seurat object
#'
#' @returns Either a subsetted Seurat object or character vector of selected cell names
#' @export
#'
#' @examples
#' \dontrun{
#' # Interactive selection with default settings (UMAP, colored by `Idents()`)
#' seu_subset <- select_cells_interactive(seurat_obj)
#'
#' # Use specific reduction and color by cell identities
#' seu_subset <- select_cells_interactive(seurat_obj, reduction = "tsne", color_by = "Idents")
#'
#' # Color by a custom metadata column
#' seu_subset <- select_cells_interactive(seurat_obj, color_by = "treatment")
#'
#' # Get just the cell names
#' selected_cell_names <- select_cells_interactive(seurat_obj, return_cells = TRUE)
#' }
select_cells_interactive <- function(seurat_obj,
                                     reduction = c("umap", "tsne", "pca"),
                                     color_by = "Idents",
                                     return_cells = FALSE) {

  # Check requirements
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("shiny package is required for interactive selection")
  }

  if (!interactive()) {
    stop("This function only works in interactive R sessions")
  }

  # Validate reduction argument
  reduction <- match.arg(reduction)

  if (is.null(seurat_obj@reductions[[reduction]])) {
    stop("Reduction '", reduction, "' not found. Please compute it first.")
  }

  # Validate color_by argument
  if (color_by != "Idents" && !color_by %in% base::colnames(seurat_obj@meta.data)) {
    stop("color_by '", color_by, "' not found in metadata. Available columns: ",
         base::paste(base::colnames(seurat_obj@meta.data), collapse = ", "))
  }

  # Get reduction coordinates and metadata
  coords <- base::as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings[, 1:2])
  base::colnames(coords) <- c("dim1", "dim2")
  coords$cell_name <- base::rownames(coords)

  # Add coloring variable
  if (color_by == "Idents") {
    coords$color_var <- base::as.character(Seurat::Idents(seurat_obj))
    color_label <- "Identity"
  } else {
    coords$color_var <- base::as.character(seurat_obj@meta.data[[color_by]])
    color_label <- color_by
  }

  # Shiny UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("Interactive Cell Selection"),

    shiny::fluidRow(
      shiny::column(3,
                    shiny::wellPanel(
                      shiny::h4("Controls"),
                      shiny::actionButton("select_btn", "Select Brushed Cells",
                                          class = "btn-primary", width = "100%"),
                      shiny::br(), shiny::br(),
                      shiny::actionButton("deselect_btn", "Deselect Brushed Cells",
                                          class = "btn-warning", width = "100%"),
                      shiny::br(), shiny::br(),
                      shiny::actionButton("clear_btn", "Clear All",
                                          class = "btn-danger", width = "100%"),
                      shiny::br(), shiny::br(),
                      shiny::actionButton("done_btn", "Done",
                                          class = "btn-success", width = "100%"),
                      shiny::hr(),
                      shiny::h5("Plot Settings:"),
                      shiny::p(base::paste("Reduction:", base::toupper(reduction))),
                      shiny::p(base::paste("Colored by:", color_label)),
                      shiny::hr(),
                      shiny::h5("Instructions:"),
                      shiny::p("1. Click and drag to brush cells"),
                      shiny::p("2. Click 'Select' or 'Deselect'"),
                      shiny::p("3. Repeat as needed"),
                      shiny::p("4. Click 'Done' when finished"),
                      shiny::hr(),
                      shiny::textOutput("selection_count")
                    )
      ),

      shiny::column(9,
                    shiny::plotOutput("main_plot",
                                      height = "600px",
                                      brush = shiny::brushOpts(id = "plot_brush"))
      )
    )
  )

  # Shiny Server
  server <- function(input, output, session) {

    # Reactive values to track selection
    values <- shiny::reactiveValues(
      selected = base::rep(FALSE, base::nrow(coords))
    )

    # Main plot
    output$main_plot <- shiny::renderPlot({

      plot_data <- coords
      plot_data$selected <- values$selected

      # Base plot with color variable
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = dim1, y = dim2)) +
        ggplot2::geom_point(ggplot2::aes(color = color_var),
                            size = 1.5, alpha = 0.7) +
        ggplot2::labs(
          x = base::paste0(base::toupper(reduction), "_1"),
          y = base::paste0(base::toupper(reduction), "_2"),
          color = color_label
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "right",
          plot.title = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::ggtitle("Select Cells by Brushing")

      # Add selected cells overlay
      if (base::any(values$selected)) {
        selected_data <- plot_data[plot_data$selected, ]
        p <- p + ggplot2::geom_point(
          data = selected_data,
          ggplot2::aes(x = dim1, y = dim2),
          color = "red", size = 2.5, shape = 1, stroke = 1.5
        )
      }

      p
    })

    # Selection count display
    output$selection_count <- shiny::renderText({
      n_selected <- base::sum(values$selected)
      base::paste("Selected:", n_selected, "cells")
    })

    # Select brushed cells
    shiny::observeEvent(input$select_btn, {
      if (!is.null(input$plot_brush)) {
        brushed <- shiny::brushedPoints(coords, input$plot_brush,
                                        xvar = "dim1", yvar = "dim2")
        if (base::nrow(brushed) > 0) {
          # Mark brushed cells as selected
          indices <- base::match(brushed$cell_name, coords$cell_name)
          values$selected[indices] <- TRUE
        }
      }
    })

    # Deselect brushed cells
    shiny::observeEvent(input$deselect_btn, {
      if (!is.null(input$plot_brush)) {
        brushed <- shiny::brushedPoints(coords, input$plot_brush,
                                        xvar = "dim1", yvar = "dim2")
        if (base::nrow(brushed) > 0) {
          # Mark brushed cells as not selected
          indices <- base::match(brushed$cell_name, coords$cell_name)
          values$selected[indices] <- FALSE
        }
      }
    })

    # Clear all selections
    shiny::observeEvent(input$clear_btn, {
      values$selected <- base::rep(FALSE, base::nrow(coords))
    })

    # Done - return results
    shiny::observeEvent(input$done_btn, {
      shiny::stopApp(values$selected)
    })
  }

  # Run the app and get selection
  base::message("Starting interactive cell selector...")
  selected_logical <- shiny::runApp(shiny::shinyApp(ui, server))

  # Process results
  if (is.null(selected_logical) || !any(selected_logical)) {
    message("No cells selected")
    return(NULL)
  }

  selected_cells <- coords$cell_name[selected_logical]
  message("Selected ", length(selected_cells), " cells")

  if (return_cells) {
    return(selected_cells)
  } else {
    return(subset(seurat_obj, cells = selected_cells))
  }
}
