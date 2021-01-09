library(shiny)
library(shinydashboard)
library(plotly)
library(cowplot)
library(SingleCellExperiment)
library(IBRAP)

ui <- dashboardPage(
  dashboardHeader(title = 'IBRAP'),
  dashboardSidebar(
    fluidPage(
      h4('Upload previously analysed data from IBRAP:'),
      fileInput(inputId = 'rds_file', label = 'RDS file upload', multiple = FALSE,
                accept = c('.rds')),
      actionButton(inputId = 'generate_metadata', 'Activate'),
      #hr(),
      #h4('Select between methods here:'),
      #uiOutput(outputId = 'norm_selector'),
      #uiOutput(outputId = 'cluster_selector'),
      hr(),
      h4('Generate a dimensionality reduced plot:'),
      textInput(inputId = 'plot.title', label = 'Plot title:', value = NULL),
      uiOutput(outputId = 'cluster.column'),
      uiOutput(outputId = 'reduction_selector'),
      br(),
      actionButton(inputId = 'plot_DR', label = 'Plot'),
      hr(),
      h4('Please save the file in png format:'),
      textInput(inputId = 'filename', label = "File name"),
      numericInput(inputId = 'width', label = 'Width', value = 600, min = 200),
      numericInput(inputId = 'height', label = 'Height', value = 600, min = 200),
      br(),
      downloadButton(outputId = 'downloadPlot')
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(
      '.myClass {
        font-size: 20px;
        line-height: 50px;
        text-align: left;
        font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
        padding: 0 15px;
        overflow: hidden;
        color: white;
      }
    '))),
    tags$script(HTML('
      $(document).ready(function() {
        $("header").find("nav").append(\'<span class="myClass"> Cluster Exploration Panel </span>\');
      })
     ')),
    box(column(width = 12, plotlyOutput(outputId = 'int_DR_plot', width = 500, height = 500)), height = "520px"),
    box(column(width = 12, plotOutput(outputId = 'benchmark', width = 500, height = 500)), height = "520px")
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
      forout_reactive$object <- x
      print('object attached')
      x <- forout_reactive$object
      #alternative.experiments <- altExpNames(x)
      #forout_reactive$alternative.experiments <- alternative.experiments
      alternative.cluster <- names(colData(x))
      forout_reactive$alternative.cluster <- alternative.cluster
      #print(paste0('detected ', alternative.experiments, ' nornalisation techniques'))
      print(paste0('detected ', alternative.cluster, ' clustering techniques'))
      #norm.bench <- list()
      clust.bench <- list()
      clust <- list()
      reduction <- list()
      cluster.names <- list()
      clust.bench <- metadata(x)[['benchmarking_clustering']]
      clust <- colData(x)
      #norm.bench[[as.character(r)]] <- metadata(altExp(x, r))[['benchmarking_normalisation']]
      rednames <- reducedDimNames(x)
        for(t in rednames) {
          reduction[[as.character(t)]] <- reducedDim(x, as.character(t))
      }
      forout_reactive$clusters <- clust
      forout_reactive$clustering_benchmarking <- clust.bench
      #forout_reactive$normalisation_benchmarking <- norm.bench
      forout_reactive$reduction <- reduction
      #forour_reactive$rednames <- rednames
    })
  })

  observeEvent(forout_reactive$reduction, {
    showNotification("Project now active", closeButton = TRUE)
  })

  # output$norm_selector <- renderUI({
  #   choices <- as.list(forout_reactive$alternative.experiments)
  #   selectInput(inputId = 'normalisation_technique',
  #               label = 'Select normalisation technique',
  #               choices = choices, multiple = FALSE)
  # })
  #
  # output$cluster_selector <- renderUI({
  #   choices <- as.list(forout_reactive$alternative.cluster)
  #   selectInput(inputId = 'cluster_technique',
  #               label = 'Select cluster technique',
  #               choices = choices, multiple = FALSE)
  # })

  output$reduction_selector <- renderUI({
    choices <- as.list(names(forout_reactive$reduction))
    selectInput(inputId = 'reduction_technique',
                label = 'Select reduction technique',
                choices = choices, multiple = FALSE)
  })

  output$cluster.column <- renderUI({
    #req(input$cluster_technique)
    choices <- names(colData(forout_reactive$object))
    selectInput(inputId = 'cluster_column',
                label = 'Select cluster column',
                choices = choices, multiple = FALSE)
  })

  observeEvent(input$plot_DR, {
    g <- forout_reactive$object
    print(input$reduction_technique)
    #print(input$x_axis_dimension)
    #print(input$y_axis_dimension)
    #print(input$normalisation_technique)
    #print(input$cluster_technique)
    print(input$cluster_column)
    print(input$plot.title)
    p <- plot_cluster_dr(object = g, reduction = as.character(input$reduction_technique),
                         pt.size = 0.1, group.by = as.character(input$cluster_column))
    ni <- p + labs(title = input$plot.title)
    pl <- ggplotly(ni, width=input$width, height=input$height)

    forout_reactive$DR_plot <- p
    forout_reactive$int_DR_plot <- pl
  })

  output$int_DR_plot <- renderPlotly({
    forout_reactive$DR_plot
  })

  output$benchmark <- renderPlot({
    req(forout_reactive$clustering_benchmarking)
    g <- forout_reactive$object
    temp <- plot_benchmarking(object = g, clust.method = as.character('seurat'), ARI = FALSE)
  })



  output$downloadPlot <- downloadHandler(
    filename = function(){paste0(as.character(input$filename),'.png',sep='')},
    content = function(file){
      ggsave(file, plot = input$DR_plot, device = "png")
    }
  )

}

shinyApp(ui, server)





