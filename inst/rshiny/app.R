library(shiny)
library(shinydashboard)
library(plotly)
library(cowplot)
library(SingleCellExperiment)
library(egg)
library(IBRAP)

shinyApp(
  ui <- dashboardPage(
    dashboardHeader(title = "IBRAP Exploration Centre"),
    dashboardSidebar(
      sidebarMenu(
        h4('Upload IBRAP item:'),
        fileInput(inputId = 'rds_file', label = 'RDS file upload', multiple = FALSE,
                  accept = c('.rds')),
        actionButton(inputId = 'generate_metadata', 'Activate'),
        uiOutput(outputId = 'active_assay'),
        uiOutput(outputId = 'reduction_selector'),
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
                          box(height = 900, width = 900, solidHeader = TRUE, status = "primary", title = 'Plotting Cell Labels',
                              splitLayout(cellWidths = c("25%", "75%"),
                                          box(height = 820, width = 300, title = 'Points Labelling',
                                              h5('Select between cell label dataframes:'),
                                              uiOutput(outputId = 'cluster_selector'),
                                              hr(),
                                              h5('Select which column to use within the dataframe:'),
                                              uiOutput(outputId = 'cluster.column'),
                                              hr(),
                                              h5('Specify the point size:'),
                                              numericInput(inputId = 'pt_size', label = 'Point size', value = 3),
                                              actionButton(inputId = 'plot_DR', label = 'Plot')
                                          ),
                                          box(align = "center", height = 820, width = 820, align = 'center',
                                              plotlyOutput(outputId = 'int_DR_plot', width = '100%', height = 800)))
                          ),
                          box(height = 1350, width = 900, solidHeader = TRUE, status = "primary", title = 'benchmarking metrics', align = 'center',
                              plotOutput(outputId = 'benchmark', width = '100%'),
                              hr(),
                              plotOutput(outputId = 'integration_benchmarking', width = '100%'),
                              hr(),
                              h3('Either 5 (ground truth available) or 3 (ground truth unavailable) benchmarking metrices are provided', align = "left"),
                              br(),
                              p(strong('No ground truth metrices:'), align = "left"),
                              p(' - ASW, Average Silhouette Width determines the separation of a cluster to its closest neighbour cluster', align = "left"),
                              p(' - Dunn Index evaluates the compactness of a cluster and its distance to its closest neighbour cluster', align = "left"),
                              p(' - connectivity determines how connecetd the cluster assignments points are to eachother', align = "left"),
                              br(),
                              p(strong('WARNING: these metrices are for guidance purposes and will not identify optimal cluster assignments without investigation, especially in respect to sub-populations of larger cell types.'), align = "left"),
                              br(),
                              p(strong('Ground truth metrices:'), align = "left"),
                              p(' - ARI, Adjusted Rand Index measures the agreement between the cluster assignments between the ground truth and novel assignments, this metric is adjusted for randomness', align = "left"),
                              p(' - NMI, Normalised Mutual Information functions similarly to ARI however it is adjusted for cluster sizes', align = "left"),
                              br(),
                              p('The boxplot (if present) describes the ASW between batches if batch correction was performed. A higher value indicates higher batch effects whilst a lower demonstrates less.', align = "left"),
                          )
                          
                          
                  )),
          tabItem(tabName = "Features",
                  column(12 ,
                         box(height = 900, width = 900, solidHeader = FALSE, status = "primary", title = 'Feature scatter plots',
                             splitLayout(cellWidths = c("25%", "75%"),
                                         box(height = 820, width = 300,
                                             uiOutput(outputId = 'assay_selector'),
                                             uiOutput(outputId = 'features_selector'),
                                             uiOutput(outputId = 'reduction_feature'),
                                             p('Please indicate the percetile range, default = 0-1 (0%-100%)'),
                                             numericInput(inputId = 'upper_percentile', value = 1, label = 'Upper percentile'),
                                             numericInput(inputId = 'lower_percentile', value = 0, label = 'Lower percentile'),
                                             numericInput(inputId = 'feature_ptsize', value = 3, label = 'Point size'),
                                             actionButton(inputId = 'plot_feature', label = 'Plot')
                                         ),
                                         box(height = 820, width = 550, align = "center",
                                             plotOutput(outputId = 'feature_plot', width = 800, height = 800)
                                         )
                             )
                         ),
                         box(height = 900, width = 900, solidHeader = FALSE, status = "primary",
                             splitLayout(cellWidths = c("25%", "75%"),
                                         box(height = 820, width = 300, title = 'Violin plot',
                                             uiOutput(outputId = 'features_selector_vln'),
                                             uiOutput(outputId = 'cluster_selector_vln'),
                                             uiOutput(outputId = 'cluster.column_vln'),
                                             actionButton(inputId = 'plot_vln', label = 'Plot')
                                         ),
                                         box(height = 820, width = 550, align = "center",
                                             plotOutput(outputId = 'vln_plot', width = 800, height = 800)
                                         )
                             )
                         ),
                         box(height = 900, width = 900, solidHeader = FALSE, status = "primary",
                             splitLayout(cellWidths = c("25%", "75%"),
                                         box(height = 820, width = 300, title = 'Dot plot',
                                             uiOutput(outputId = 'features_selector_heat'),
                                             uiOutput(outputId = 'cluster_selector_heat'),
                                             uiOutput(outputId = 'cluster.column_heat'),
                                             actionButton(inputId = 'plot_heat', label = 'Plot')
                                         ),
                                         box(height = 820, width = 550, align = "center",
                                             plotOutput(outputId = 'heat_plot', width = 800, height = 800)
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
        cat(crayon::cyan(paste0(Sys.time(), ': loading RDS file\n')))
        req(input$rds_file)
        direc <- input$rds_file
        x <- readRDS(file = as.character(direc$datapath))
        cat(crayon::cyan(paste0(Sys.time(), ': RDS file loaded \n')))
        
        for(i in names(x@methods)[2:length(names(x@methods))]) {
          
          x@methods[[i]]@cluster_assignments[['metadata']] <- x@sample_metadata
          
        }
        
        forout_reactive$obj <- x
        forout_reactive$assay.names <- names(x@methods)[2:length(names(x@methods))]
        forout_reactive$assays <- x@methods[forout_reactive$assay.names]
        
      })
    })
    
    output$active_assay <- renderUI({
      
      req(forout_reactive$assay.names)
      selectInput(inputId = 'assay', choices = forout_reactive$assay.names, label = 'Select assay:')
      
    })
    
    observeEvent(input$assay, {
      
      forout_reactive$active.assay <- forout_reactive$assays[[input$assay]]
      
    })
    
    observeEvent(forout_reactive$obj, {
      
      showNotification("Project now active", closeButton = TRUE)
      
    })
    
    output$reduction_selector <- renderUI({
      
      req(forout_reactive$active.assay)
      choices <- names(forout_reactive$active.assay@visualisation_reductions)
      selectInput(inputId = 'reduction_technique',
                  label = 'Select reduction technique',
                  choices = choices, multiple = FALSE)
      
    })
    
    observeEvent(eventExpr = input$reduction_technique, {
      
      temp <- function(x) {
        
        n <- sum(length(x)-1)
        r <- head(x, n=n)
        return(r)
        
      }
      
      red <- input$reduction_technique
      
      ff <- unlist(lapply(X = lapply(X = strsplit(x = names(forout_reactive$active.assay@cluster_assignments), split = '_'), FUN = temp), FUN = paste0, collapse = '_'))
      
      f <- unlist(lapply(X = lapply(X = strsplit(x = red, split = '_'), FUN = temp), FUN = paste0, collapse = '_'))
      
      forout_reactive$allowed_labels <- c(names(forout_reactive$active.assay@cluster_assignments)[f == ff], 'metadata')
      
    })
    
    output$cluster_selector <- renderUI({
      req(forout_reactive$active.assay)
      choices <- forout_reactive$allowed_labels
      selectInput(inputId = 'cluster_technique',
                  label = 'Select dataframe from @cluster_assignments',
                  choices = choices, multiple = FALSE)
    })
    
    output$cluster.column <- renderUI({
      req(forout_reactive$active.assay)
      selectInput(inputId = 'cluster_column',
                  label = 'Select dataframe column',
                  choices = suppressWarnings(colnames(forout_reactive$active.assay@cluster_assignments[[as.character(input$cluster_technique)]])), 
                  multiple = FALSE)
    })
    
    observeEvent(input$plot_DR, {
      g <- forout_reactive$object
      
      p <- IBRAP::plot.reduced.dim.interactive(object = forout_reactive$obj, 
                                               reduction = input$reduction_technique, 
                                               assay = input$assay, 
                                               pt.size = input$pt_size, 
                                               clust.method = input$cluster_technique, 
                                               column = input$cluster_column, 
                                               dimensions = 2)
      forout_reactive$DR_plot <- p
    })
    
    output$int_DR_plot <- renderPlotly({
      forout_reactive$DR_plot
    })
    
    output$benchmark <- renderPlot({
      req(input$cluster_technique != 'metadata')
      g <- forout_reactive$obj
      h <- names(forout_reactive$active.assay@benchmark_results$clustering[[input$cluster_technique]])
      if(length(h) > 3) {
        temp <- plot.cluster.benchmarking(object = g, 
                                          assay = input$assay, 
                                          clustering = input$cluster_technique, 
                                          ARI = TRUE)
      } else {
        temp <- plot.cluster.benchmarking(object = g, 
                                          assay = input$assay, 
                                          clustering = input$cluster_technique, 
                                          ARI = FALSE)
      }
      
    })
    
    output$integration_benchmarking <- renderPlot({
      req(forout_reactive$obj)
      
      if(is.null(pancreas_bench@methods$SCT@benchmark_results$integration)) {
        
        return(NULL)
        
      }
      
      p <- plot.integration.benchmarking(object = forout_reactive$obj, assay = input$assay)
      
    })
    
    output$features_selector <- renderUI({
      choices <- rownames(forout_reactive$obj@methods[[1]]@counts)
      suppressWarnings(selectizeInput(inputId = 'features', label = 'Identify features', 
                                      choices = choices, selected = NULL, multiple = TRUE, 
                                      options = list(create = TRUE)))
    })
    
    observeEvent(input$plot_feature, {
      g <- forout_reactive$obj
      
      p <- IBRAP::plot.features(object = g, assay = input$assay, slot = 'normalised', 
                                percentile = c(as.numeric(input$lower_percentile), as.numeric(input$upper_percentile)),
                                pt_size = input$feature_ptsize, reduction = input$reduction_technique, features = input$features)
      
      forout_reactive$feature_plot <- p
      
    })
    
    output$feature_plot <- renderPlot({
      forout_reactive$feature_plot
    })
    
    output$features_selector_vln <- renderUI({
      choices <- rownames(forout_reactive$obj@methods[[1]]@counts)
      suppressWarnings(selectizeInput(inputId = 'features_vln', label = 'Identify features',
                                      choices = choices, selected = NULL, multiple = TRUE, 
                                      options = list(create = TRUE)))
    })
    
    output$cluster_selector_vln <- renderUI({
      choices <- names(forout_reactive$active.assay@cluster_assignments)
      selectInput(inputId = 'cluster_technique_vln',
                  label = 'Select cluster technique',
                  choices = choices, multiple = FALSE)
    })
    
    output$cluster.column_vln <- renderUI({
      selectInput(inputId = 'cluster_column_vln',
                  label = 'Select cluster column',
                  choices = suppressWarnings(colnames(forout_reactive$active.assay@cluster_assignments[[input$cluster_technique_vln]])), multiple = FALSE)
    })
    
    observeEvent(input$plot_vln, {
      g <- forout_reactive$obj
      
      p <- plot.vln(object = g, assay = input$assay, 
                    slot = 'normalised', 
                    features = input$features_vln, 
                    group.by = forout_reactive$active.assay@cluster_assignments[[input$cluster_technique_vln]][[input$cluster_column_vln]])
      
      forout_reactive$vln_plot <- p
    })
    
    output$vln_plot <- renderPlot({
      forout_reactive$vln_plot
    })
    
    output$features_selector_heat <- renderUI({
      
      choices <- rownames(forout_reactive$obj@methods[[1]]@counts)
      suppressWarnings(selectizeInput(inputId = 'features_heat', label = 'Identify features', 
                                      choices = choices, selected = NULL, multiple = TRUE, 
                                      options = list(create = TRUE)))
      
    })
    
    output$cluster_selector_heat <- renderUI({
      choices <- names(forout_reactive$active.assay@cluster_assignments)
      selectInput(inputId = 'cluster_selector_heat',
                  label = 'Select cluster technique',
                  choices = choices, multiple = FALSE)
    })
    
    output$cluster.column_heat <- renderUI({
      
      selectInput(inputId = 'cluster.column_heat',
                  label = 'Select cluster column',
                  choices = suppressWarnings(colnames(forout_reactive$active.assay@cluster_assignments[[input$cluster_selector_heat]])), multiple = FALSE)
      
    })
    
    output$heat_plot <- renderPlot({
      
      req(input$plot_heat)
      
      plot.dot.plot(object = forout_reactive$obj,
                    assay = input$assay,
                    slot = 'normalised',
                    features = input$features_heat, 
                    clust.method = input$cluster_selector_heat, 
                    column = input$cluster.column_heat)
      
    })
    
    
  }
)
