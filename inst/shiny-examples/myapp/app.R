library(shiny)
library(shinydashboard)
library(plotly)
library(cowplot)
library(SingleCellExperiment)
library(egg)

shinyApp(
  ui <- dashboardPage(
    dashboardHeader(title = "IBRAP Exploration Centre"),
    dashboardSidebar(
      sidebarMenu(
        h4('Upload IBRAP item:'),
        fileInput(inputId = 'rds_file', label = 'RDS file upload', multiple = FALSE,
                  accept = c('.rds')),
        actionButton(inputId = 'generate_metadata', 'Activate'),
        hr(),
        h4('Select tab:'),
        menuItem("Clustering", tabName = "Clustering"),
        menuItem("Features", tabName = "Features")
      )
    ),
    dashboardBody(
      fluidRow(
      tabItems(
        tabItem(tabName = "Clustering",
                column( 12,
                        box(height = 900, width = 900, solidHeader = TRUE, status = "primary", title = 'Clustering plots',
                            splitLayout(cellWidths = c("25%", "75%"),
                                        box(height = 820, width = 300, title = 'Select metadata',
                                            h4('Select between methods here:'),
                                            uiOutput(outputId = 'cluster_selector'),
                                            uiOutput(outputId = 'group_by'),
                                            hr(),
                                            h4('Generate a dimensionality reduced plot:'),
                                            uiOutput(outputId = 'cluster.column'),
                                            uiOutput(outputId = 'reduction_selector'),
                                            selectInput(inputId = 'cluster_dimensions', label = '2D or 3D?', choices = c('2D', '3D'), multiple = FALSE),
                                            numericInput(inputId = 'pt_size', label = 'Point size', value = 5, min = 0.1, max = 10),
                                            br(),
                                            actionButton(inputId = 'plot_DR', label = 'Plot')),
                                        box(height = 820, width = 550, align = "center",
                                            plotlyOutput(outputId = 'int_DR_plot', width = 1100, height = 800)))
                        ),
                        box(height = 900, width = 900, solidHeader = TRUE, status = "primary", title = 'Bar plots', align = 'center',
                            plotOutput(outputId = 'bar_plot', width = 1400, height = 800)
                        ),
                        box(height = 900, width = 900, solidHeader = TRUE, status = "primary", title = 'benchmarking metrics', align = 'center',
                            plotOutput(outputId = 'benchmark', width = 1400, height = 800)
                        )

                
        )),
        tabItem(tabName = "Features",
                column(12 ,
                        box(height = 900, width = 900, solidHeader = FALSE, status = "primary", title = 'Feature scatter plots',
                            splitLayout(cellWidths = c("25%", "75%"),
                                        box(height = 820, width = 300, title = 'Feature plot',
                                            uiOutput(outputId = 'assay_selector'),
                                            uiOutput(outputId = 'features_selector'),
                                            uiOutput(outputId = 'reduction_feature'),
                                            actionButton(inputId = 'plot_feature', label = 'Plot')
                                            ),
                                        box(height = 820, width = 550, align = "center",
                                            plotOutput(outputId = 'feature_plot', width = 800, height = 800)
                                            )
                         )
                     ),
                     box(height = 900, width = 900, solidHeader = FALSE, status = "primary", title = 'Feature violin plots',
                         splitLayout(cellWidths = c("25%", "75%"),
                                     box(height = 820, width = 300, title = 'Violin plot',
                                         selectInput(inputId = 'graph_type', label = 'Please select graph:', 
                                                     choices = c('Violin plot', 'Heatmap'), multiple = FALSE),
                                         uiOutput(outputId = 'assay_selector_vln'),
                                         uiOutput(outputId = 'features_selector_vln'),
                                         uiOutput(outputId = 'cluster_selector_vln'),
                                         uiOutput(outputId = 'cluster.column_vln'),
                                         actionButton(inputId = 'plot_vln', label = 'Plot')
                                     ),
                                     box(height = 820, width = 550, align = "center",
                                         plotOutput(outputId = 'vln_plot', width = 800, height = 800)
                                     )
                         )
                     )
            )
          )
        )
      )
    )
  ),
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
        metadata(x)[['clustering']][['metadata']] <- as.data.frame(colData(x))
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
    
    output$group_by <- renderUI({
      req(forout_reactive$object)
      obj <- forout_reactive$object
      dt <- colData(obj)
      choices <- colnames(dt)[grepl('factor|logical|character',sapply(dt,class))]
      selectInput(inputId = 'categorical_metadata', label = 'Select categorical grouping', choices = choices, multiple = FALSE)
    })
    
    output$bar_plot <- renderPlot({
      req(input$categorical_metadata)
      req(input$cluster_technique)
      req(input$cluster_column)
      x <- forout_reactive$object
      plot.barplot(object = x, x.value = x[[input$categorical_metadata]], metadata(x)[['clustering']][[input$cluster_technique]][[input$cluster_column]])
    })
    
    output$reduction_selector <- renderUI({
      choices <- unlist(forout_reactive$rednames)
      selectInput(inputId = 'reduction_technique',
                  label = 'Select reduction technique',
                  choices = choices, multiple = FALSE)
    })
    
    
    output$cluster.column <- renderUI({
      selectInput(inputId = 'cluster_column',
                  label = 'Select cluster column',
                  choices = unlist(names(forout_reactive$clustering_assignments[[as.character(input$cluster_technique)]])), multiple = FALSE)
    })

    observeEvent(input$plot_DR, {
      g <- forout_reactive$object
      print(input$reduction_technique)
      print(input$cluster_technique)
      print(input$cluster_column)
      if(input$cluster_dimensions == '3D') {
        p <- plot.reduced.dim(object = forout_reactive$object, 
                              reduction = input$reduction_technique, 
                              pt.size = input$pt_size, metadata.access = 'clustering', 
                              sub.access = input$cluster_technique, 
                              group.by = input$cluster_column, 
                              dimensions = 3)
      } else if (input$cluster_dimensions == '2D') {
        p <- plot.reduced.dim(object = forout_reactive$object, 
                              reduction = input$reduction_technique, 
                              pt.size = input$pt_size, metadata.access = 'clustering', 
                              sub.access = input$cluster_technique, 
                              group.by = input$cluster_column, 
                              dimensions = 2)
      }
      forout_reactive$DR_plot <- p
    })
    
    output$int_DR_plot <- renderPlotly({
      forout_reactive$DR_plot
    })
    
    output$benchmark <- renderPlot({
      req(forout_reactive$clustering_benchmarking)
      g <- forout_reactive$object
      if(input$cluster_technique == 'metadata') {
        return(NULL)
      } else {
        temp <- plot.benchmarking(object = g, clust.method = input$cluster_technique, ARI = TRUE)
      }
    })
    
    output$assay_selector <- renderUI({
      choices <- names(forout_reactive$assays)
      selectInput(inputId = 'assays',
                  label = 'Select assay',
                  choices = choices, multiple = FALSE)
    })
    
    output$features_selector <- renderUI({
      choices <- rownames(forout_reactive$object)
      
      selectizeInput(inputId = 'features', label = 'Identify features', 
                     choices = choices, selected = NULL, multiple = TRUE, 
                     options = list(create = TRUE))
    })
    
    output$reduction_feature <- renderUI({
      choices <- unlist(forout_reactive$rednames)
      selectInput(inputId = 'reduction_feature',
                  label = 'Select reduction technique',
                  choices = choices, multiple = FALSE)
    })
    
    output$vln_option <- renderUI({
      
    })
    
    observeEvent(input$plot_feature, {
      g <- forout_reactive$object
      print(g)
      print(input$assays)
      print(input$features)
      print(input$reduction_feature)
      p <- plot.features.multiple(object = g, 
                                  assay = input$assays, 
                                  reduction = input$reduction_feature, 
                                  lab.key = 'reduction', features = input$features)
      forout_reactive$feature_plot <- p
      
  })
    
    output$feature_plot <- renderPlot({
      forout_reactive$feature_plot
    })
    
    output$assay_selector_vln <- renderUI({
      choices <- names(forout_reactive$assays)
      selectInput(inputId = 'assays_vln',
                  label = 'Select assay',
                  choices = choices, multiple = FALSE)
    })
    
    output$features_selector_vln <- renderUI({
      choices <- rownames(forout_reactive$object)
      selectizeInput(inputId = 'features_vln', label = 'Identify features', 
                     choices = choices, selected = NULL, multiple = TRUE, 
                     options = list(create = TRUE))
    })
    #
    output$cluster_selector_vln <- renderUI({
      choices <- unlist(forout_reactive$alternative.cluster)
      selectInput(inputId = 'cluster_technique_vln',
                  label = 'Select cluster technique',
                  choices = choices, multiple = FALSE)
    })
    #
    output$cluster.column_vln <- renderUI({
      selectInput(inputId = 'cluster_column_vln',
                  label = 'Select cluster column',
                  choices = unlist(names(forout_reactive$clustering_assignments[[as.character(input$cluster_technique_vln)]])), multiple = FALSE)
    })
    
    observeEvent(input$plot_vln, {
      g <- forout_reactive$object
      print(g)
      print(input$graph_type)
      if(input$graph_type == 'Violin plot') {
        print(input$assays_vln)
        print(input$features_vln)
        print(input$cluster_technique_vln)
        print(input$cluster_column_vln)
        p <- plot.vln(object = g, 
                      assay = input$assays_vln, 
                      features = input$features_vln, 
                      group.by = metadata(g)[['clustering']][[input$cluster_technique_vln]][[input$cluster_column_vln]])
        forout_reactive$vln_plot <- p
      } else if(input$graph_type == 'Heatmap') {
        print(input$assays_vln)
        print(input$features_vln)
        p <- plot.heatmap(object = g, assay = input$assays_vln, features = input$features_vln)
        forout_reactive$vln_plot <- p
        }
    })
    
    output$vln_plot <- renderPlot({
      forout_reactive$vln_plot
    })
    
  }
)
