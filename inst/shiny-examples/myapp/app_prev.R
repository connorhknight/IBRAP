library(shiny)
library(shinydashboard)
library(plotly)
library(cowplot)
library(SingleCellExperiment)
# library(IBRAP)

ui <- dashboardPage(
  dashboardHeader(title = 'IBRAP Visual Centre'),
  dashboardSidebar(
    fluidPage(
      h4('Upload previously analysed data from IBRAP:'),
      fileInput(inputId = 'rds_file', label = 'RDS file upload', multiple = FALSE,
                accept = c('.rds')),
      actionButton(inputId = 'generate_metadata', 'Activate'),
      hr(),
      menuItem("Clustering", tabName = "Clustering", icon = icon('glyphicon')),
      menuItem("Features", tabName = "Features", icon = icon('glyphicon'))
    )
  ),
  dashboardBody(
    body <- dashboardBody(
      tabItems(
        tabItem(tabName = "Clustering",
                h2("Cluster_exploration"),
                box(      h4('Select between methods here:'),
                          uiOutput(outputId = 'cluster_selector'),
                          hr(),
                          h4('Generate a dimensionality reduced plot:'),
                          uiOutput(outputId = 'cluster.column'),
                          uiOutput(outputId = 'reduction_selector'),
                          selectInput(inputId = 'dimensions', label = '2D or 3D?', choices = c('2D', '3D'), multiple = FALSE),
                          numericInput(inputId = 'pt_size', label = 'Point size', value = 5, min = 0.1, max = 10),
                          br(),
                          actionButton(inputId = 'plot_DR', label = 'Plot'),
                          hr(),
                          h4('Please save the file in png format:'),
                          textInput(inputId = 'filename', label = "File name"),
                          br(),
                          downloadButton(outputId = 'downloadPlot')),
                box(column(width = 12, plotlyOutput(outputId = 'int_DR_plot', width = 500, height = 500)), height = "520px"),
                box(column(width = 12, plotOutput(outputId = 'benchmark', width = 500, height = 500)), height = "520px")
        ),
        
        tabItem(tabName = "Features",
                h2("Feature_exploration")
        )
      )
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)

  forout_reactive <- reactiveValues()

  observeEvent(input$generate_metadata, {
    withProgress(message = 'Loading RDS file', {
      print('start')
      req(input$rds_file)
      direc <- input$rds_file
      print('loading .rds file')
      print(input$rds_file)
      x <- readRDS(file = as.character(direc$datapath))
      print('.rds file loaded')
      print('object attached')
      metadata(x)[['clustering']][['metadata']] <- colData(x)
      alternative.cluster <- names(metadata(x)[['clustering']])
      cluster.info <- metadata(x)[['clustering']]
      forout_reactive$alternative.cluster <- alternative.cluster
      print(paste0('detected ', alternative.cluster, ' clustering techniques'))
      clust.bench <- list()
      clust <- list()
      reduction <- list()
      cluster.names <- list()
      assays.list <- list()
      clust.bench <- metadata(x)[['benchmarking_clustering']]
      for(n in names(clust.bench)) {
        r <- clust.bench[[n]]
        clust.bench[[n]] <- r[ , colSums(is.na(r)) == 0]
      }
      rednames <- reducedDimNames(x)
      for(t in rednames) {
        reduction[[as.character(t)]] <- reducedDim(x, as.character(t))
      }
      assnames <- assayNames(x)
      for(t in assnames) {
        assays.list[[as.character(t)]] <- assay(x, as.character(t))
      }
      forout_reactive$clustering_benchmarking <- clust.bench
      forout_reactive$clustering_assignments <- cluster.info
      forout_reactive$reduction <- reduction
      forout_reactive$rednames <- rednames
      forout_reactive$assays <- assays.list
      forout_reactive$object <- x
      rm(assays.list, clust, clust.bench, cluster.names, 
         r, reduction, x, alternative.cluster, assnames, 
         n, rednames, t, cluster.info)
    })
  })

  observeEvent(forout_reactive$reduction, {
    showNotification("Project now active", closeButton = TRUE)
  })

  output$cluster_selector <- renderUI({
    choices <- unlist(forout_reactive$alternative.cluster)
    selectInput(inputId = 'cluster_technique',
                label = 'Select cluster technique',
                choices = choices, multiple = FALSE)
  })

  output$reduction_selector <- renderUI({
    choices <- unlist(forout_reactive$rednames)
    selectInput(inputId = 'reduction_technique',
                label = 'Select reduction technique',
                choices = choices, multiple = FALSE)
  })

  output$cluster.column <- renderUI({
    print('success')
    print(input$cluster_technique)
    choices <- names(forout_reactive$clustering_assignments[[as.character(input$cluster_technique)]])
    print(choices)
    selectInput(inputId = 'cluster_column',
                label = 'Select cluster column',
                choices = choices, multiple = FALSE)
  })

  observeEvent(input$plot_DR, {
    g <- forout_reactive$object
    print(input$reduction_technique)
    print(input$cluster_technique)
    print(input$cluster_column)
    p <- plot.reduced.dim(object = forout_reactive$object, 
                          reduction = input$reduction_technique, 
                          pt.size = input$pt_size, metadata.access = 'clustering', 
                          sub.access = input$cluster_technique, 
                          group.by = input$cluster_column, 
                          dimensions = input$dimensions)
    forout_reactive$DR_plot

  })

  output$int_DR_plot <- renderPlotly({
    forout_reactive$DR_plot
  })

  output$benchmark <- renderPlot({
    req(forout_reactive$clustering_benchmarking)
    g <- forout_reactive$object
    temp <- plot_benchmarking(object = g, clust.method = input$cluster_technique, ARI = TRUE)
  })

  # output$downloadPlot <- downloadHandler(
  #   filename = function(){paste0(as.character(input$filename),'.png',sep='')},
  #   content = function(file){
  #     ggsave(file, plot = input$DR_plot, device = "png")
  #   }
  # )
}

shinyApp(ui, server)





