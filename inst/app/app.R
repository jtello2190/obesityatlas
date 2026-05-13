library(shiny)
library(obesityatlas)
library(dplyr)
library(tidyr)
library(ggplot2)

# Read the data path injected by run_app()
data_path <- getOption("obesityatlas.data_path")
if (is.null(data_path)) {
  stop(
    "This app must be launched via obesityatlas::run_app(data_path = ...). ",
    "Do not source app.R directly."
  )
}

atlas <- obesityatlas::load_atlas(data_path)

ui <- fluidPage(
  titlePanel("Mouse multi-tissue obesity atlas"),

  tabsetPanel(
    tabPanel(
      "UMAP",
      sidebarLayout(
        sidebarPanel(
          selectInput("tissue_umap", "Tissue:",
                      choices = c("All", sort(unique(atlas$obs$tissue)))),
          selectInput("color_by", "Color by:",
                      choices = c("cell_type", "leiden", "diet",
                                  "phase", "S_score", "G2M_score",
                                  "total_counts"),
                      selected = "cell_type"),
          sliderInput("point_size", "Point size:", 0.1, 2, 0.5, 0.1),
          width = 3
        ),
        mainPanel(
          plotOutput("umap_plot", height = "600px"),
          width = 9
        )
      )
    ),

    tabPanel(
      "Markers",
      sidebarLayout(
        sidebarPanel(
          selectInput("tissue_mk", "Tissue:",
                      choices = sort(unique(atlas$markers$tissue))),
          uiOutput("cluster_picker"),
          numericInput("top_n", "Top N markers:", 20, min = 5, max = 100),
          width = 3
        ),
        mainPanel(
          DT::dataTableOutput("marker_table"),
          width = 9
        )
      )
    ),

    tabPanel(
      "Composition",
      plotOutput("comp_plot", height = "500px")
    )
  )
)

server <- function(input, output, session) {

  filtered_obs <- reactive({
    df <- atlas$obs
    if (input$tissue_umap != "All") {
      df <- dplyr::filter(df, tissue == input$tissue_umap)
    }
    df
  })

  output$umap_plot <- renderPlot({
    df <- filtered_obs()
    col_var <- input$color_by

    p <- ggplot(df, aes(umap_1, umap_2, color = .data[[col_var]])) +
      geom_point(size = input$point_size, alpha = 0.6) +
      theme_minimal() +
      labs(x = "UMAP 1", y = "UMAP 2",
           title = if (input$tissue_umap == "All") "All tissues"
           else input$tissue_umap)

    if (is.numeric(df[[col_var]])) {
      p + scale_color_viridis_c()
    } else {
      p + guides(color = guide_legend(override.aes = list(size = 3)))
    }
  })

  output$cluster_picker <- renderUI({
    clusters <- atlas$markers |>
      dplyr::filter(tissue == input$tissue_mk) |>
      dplyr::pull(leiden) |>
      unique() |>
      sort()
    selectInput("cluster_mk", "Leiden cluster:", choices = clusters)
  })

  output$marker_table <- DT::renderDataTable({
    req(input$cluster_mk)
    obesityatlas::top_markers(
      atlas$markers,
      input$tissue_mk,
      input$cluster_mk,
      top_n = input$top_n
    ) |>
      dplyr::mutate(
        dplyr::across(c(avg_log2FC, p_val_adj, pct.1, pct.2),
                      ~ signif(.x, 3))
      )
  })

  output$comp_plot <- renderPlot({
    comp <- obesityatlas::cell_composition(
      atlas$obs,
      by = c("tissue", "diet")
    )

    ggplot(comp, aes(tissue, prop, fill = cell_type)) +
      geom_col() +
      facet_wrap(~ diet) +
      theme_minimal() +
      labs(y = "Proportion", x = "Tissue") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
}

shinyApp(ui, server)
