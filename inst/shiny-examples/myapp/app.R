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
                          box(height = 900, width = 900, solidHeader = TRUE, status = "primary", title = 'Clustering plots',
                              splitLayout(cellWidths = c("25%", "75%"),
                                          box(height = 820, width = 300, title = 'Select metadata',
                                              h4('Select between methods here:'),
                                              uiOutput(outputId = 'cluster_selector'),
                                              hr(),
                                              h4('Generate a dimensionality reduced plot:'),
                                              uiOutput(outputId = 'cluster.column'),
                                              numericInput(inputId = 'pt_size', label = 'Point size', value = 5, min = 0.01, max = 10),
                                              br(),
                                              actionButton(inputId = 'plot_DR', label = 'Plot')),
                                          box(height = 820, width = 550, align = "center",
                                              plotlyOutput(outputId = 'int_DR_plot', width = 1075, height = 800)))
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
                         box(height = 900, width = 900, solidHeader = FALSE, status = "primary", title = 'Feature heatmap',
                             splitLayout(cellWidths = c("25%", "75%"),
                                         box(height = 820, width = 300, title = 'Violin plot',
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
    
    output$cluster_selector <- renderUI({
      req(forout_reactive$active.assay)
      choices <- names(forout_reactive$active.assay@cluster_assignments)
      selectInput(inputId = 'cluster_technique',
                  label = 'Select cell assignments ',
                  choices = choices, multiple = FALSE)
    })
    
    output$cluster.column <- renderUI({
      req(forout_reactive$active.assay)
      selectInput(inputId = 'cluster_column',
                  label = 'Select cell assignment column',
                  choices = unlist(names(forout_reactive$active.assay@cluster_assignments[[as.character(input$cluster_technique)]])), 
                  multiple = FALSE)
    })
    
    output$reduction_selector <- renderUI({
      req(forout_reactive$active.assay)
      choices <- names(forout_reactive$active.assay@visualisation_reductions)
      selectInput(inputId = 'reduction_technique',
                  label = 'Select reduction technique',
                  choices = choices, multiple = FALSE)
    })
    
    observeEvent(input$plot_DR, {
      g <- forout_reactive$object
      print(input$reduction_technique)
      print(input$cluster_technique)
      print(input$cluster_column)
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
      h <- names(forout_reactive$active.assay@benchmark_results[[input$cluster_technique]])
      if(length(h) > 3) {
        temp <- IBRAP::plot.benchmarking(object = g, 
                                         assay = input$assay, 
                                         clustering = input$cluster_technique, 
                                         ARI = TRUE)
      } else {
        temp <- IBRAP::plot.benchmarking(object = g, 
                                         assay = input$assay, 
                                         clustering = input$cluster_technique, 
                                         ARI = FALSE)
      }
      
    })
    
    output$features_selector <- renderUI({
      choices <- rownames(forout_reactive$obj@methods[[1]]@counts)
      selectizeInput(inputId = 'features', label = 'Identify features', 
                     choices = choices, selected = NULL, multiple = TRUE, 
                     options = list(create = TRUE))
    })
    
    observeEvent(input$plot_feature, {
      g <- forout_reactive$obj
      print(g)
      print(input$assay)
      print(input$features)
      print(input$reduction_technique)
      p <- IBRAP::plot.features.multiple(object = g, assay = input$assay, slot = 'normalised', 
                                         reduction = input$reduction_technique, features = input$features)
      forout_reactive$feature_plot <- p
      
    })
    
    output$feature_plot <- renderPlot({
      forout_reactive$feature_plot
    })
    
    output$features_selector_vln <- renderUI({
      choices <- rownames(forout_reactive$obj@methods[[1]]@counts)
      selectizeInput(inputId = 'features_vln', label = 'Identify features',
                     choices = choices, selected = NULL, multiple = TRUE, 
                     options = list(create = TRUE))
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
                  choices = colnames(forout_reactive$active.assay@cluster_assignments[[input$cluster_technique_vln]]), multiple = FALSE)
    })
    
    observeEvent(input$plot_vln, {
      g <- forout_reactive$obj
      print(g)
      p <- IBRAP::plot.vln(object = g, assay = input$assay, 
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
      selectizeInput(inputId = 'features_heat', label = 'Identify features', 
                     choices = choices, selected = NULL, multiple = TRUE, 
                     options = list(create = TRUE))
      
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
                  choices = colnames(forout_reactive$active.assay@cluster_assignments[[input$cluster_selector_heat]]), multiple = FALSE)
      
    })
    
    output$heat_plot <- renderPlot({
      
      req(input$plot_heat)
      
      IBRAP::plot.dot.plot(object = forout_reactive$obj,
                           assay = input$assay,
                           slot = 'normalised',
                           features = input$features_heat, 
                           clust.method = input$cluster_selector_heat, 
                           column = input$cluster.column_heat)
      
    })
    
    
  }
)
