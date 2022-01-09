#' @name plot.reduced.dim.interactive
#' @aliases plot.reduced.dim.interactive
#' 
#' @title Interactive plot of reduced dimensions and labels
#'
#' @param object An IBRAP S4 class object
#' @param reduction Character. Which reduction to use to display points
#' @param assay Character. Which assay within the object to access
#' @param clust.method Character. Which clustering method to access, supply the name of any dataframe contained in the cluster_assignments sections. If you wish to access metadata, just specify 'metadata'
#' @param column Character. Which column to access within the supplied clust.column
#' @param pt.size Numeric. What should the point size be
#' @param dimensions Numeric. How many dimensions should be displayed, can be either 2 or 3
#' @param cells Numeric. Which cells should be subset for plotting, Default = NULL
#' 
#' @return An interactive plot of reduced dimensions annotated with assignments
#'
#' @export plot.reduced.dim.interactive

plot.reduced.dim.interactive <- function(object,
                                         reduction,
                                         assay,
                                         clust.method,
                                         column,
                                         pt.size=5, 
                                         dimensions,
                                         cells = NULL,
                                         ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(reduction)) {
    
    stop('reduction must be character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('assay does not exist\n')
    
  }
  
  if(!reduction %in% names(c(object@methods[[assay]]@computational_reductions, 
                             object@methods[[assay]]@visualisation_reductions, 
                             object@methods[[assay]]@integration_reductions))) {
    
    stop('reduction does not exist\n')
    
  }
  
  if(!is.character(clust.method)) {
    
    stop('clust.method must be character string\n')
    
  }
  
  if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
    
    if(!clust.method == 'metadata') {
      
      stop('clust.method does not exist\n')
      
    }
    
  }
  
  if(!is.character(column)) {
    
    stop('column must be character string\n')
    
  }
  
  if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
    
    if(!column %in% colnames(object@sample_metadata)) {
      
      stop('column:', column, ', does not exist in clust.method: ', clust.method, '\n')
      
    }
    
  }
  
  if(!is.numeric(pt.size)) {
    
    stop('pt.size must be numerical \n')
    
  }
  
  if(!is.numeric(dimensions)) {
    
    stop('dimensions must be numerical: 2 or 3 \n')
    
  }
  
  reduction.list <- list()
  red.names <- c(names(object@methods[[assay]]@computational_reductions), 
                 names(object@methods[[assay]]@integration_reductions),
                 names(object@methods[[assay]]@visualisation_reductions))
  
  for(i in red.names) {
    
    if(i %in% names(object@methods[[assay]]@computational_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@computational_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@integration_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@integration_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@visualisation_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@visualisation_reductions[[i]]
      
    }
    
  }
  
  results <- as.data.frame(reduction.list[[reduction]])
  
  if(clust.method == 'metadata') {
    
    project.met <- object@sample_metadata
    
    if(dimensions == 2) {
      
      if(ncol(results) < 2) {
        
        stop('Not enough dimensions present in ', reduction, '\n')
        
      }
      
      results <- results[,1:2]
      
    } else if (dimensions == 3) {
      
      if(ncol(results) < 3) {
        
        stop('Not enough dimensions present in ', reduction, '\n')
        
      }
      
      results <- results[,1:3]

    }
    
    orig.names <- colnames(results)
    
    results <- cbind(results, project.met[,column])
    
    colnames(results) <- c(orig.names, 'variable')
    
    rownames(results) <- colnames(object)
    
    if(!is.null(cells)) {
      
      results <- results[cells,]
      
    }
    
  } else {
    
    if(dimensions == 2) {
      
      if(ncol(results) < 2) {
        
        stop('Not enough dimensions present in ', reduction, '\n')
        
      } else {
        
        results <- results[,1:2]
        
      }
      
    } else if (dimensions == 3) {
      
      if(ncol(results) < 3) {
        
        stop('Not enough dimensions present in ', reduction, '\n')
        
      } else {
        
        results <- results[,1:3]
        
      }
      
    }
    
    orig.names <- colnames(results)

    project.met <- object@sample_metadata

    assay.met <- object@methods[[assay]]@cluster_assignments[[clust.method]]

    assay.met <- assay.met[match(rownames(project.met), rownames(assay.met)),]

    meta <- cbind(project.met, assay.met)
    
    results <- cbind(results,meta[,column])

    colnames(results) <- c(orig.names, 'variable')

    rownames(results) <- colnames(object)
    
    if(!is.null(cells)) {
      
      results <- results[cells,]
      
    }
    
  }
  
  if(dimensions == 3) {
    
    p <- suppressMessages(suppressWarnings(plotly::plot_ly(data = results,
                                                           x = as.formula(paste0('~', colnames(results)[1])), 
                                                           y = as.formula(paste0('~', colnames(results)[2])),
                                                           z = as.formula(paste0('~', colnames(results)[3])), 
                                                           color = as.formula(paste0('~',colnames(results)[4])), 
                                                           colors = scales::hue_pal()(length(unique(results[,'variable']))),
                                                           mode = "markers", 
                                                           marker = list(size = pt.size, width=0.5), 
                                                           text=as.formula(paste0('~',colnames(results)[4])), 
                                                           hoverinfo="text", plot_bgcolor = 'black', ...)))
    
    print(p)
    
  } else if (dimensions == 2) {
    
    p <- suppressMessages(suppressWarnings(plotly::plot_ly(data = as.data.frame(results), 
                                                           x = as.formula(paste0('~', colnames(results)[1])), 
                                                           y = as.formula(paste0('~', colnames(results)[2])), 
                                                           color = as.formula(paste0('~',colnames(results)[3])),
                                                           colors = scales::hue_pal()(length(unique(results[,'variable']))), 
                                                           mode = "markers", 
                                                           marker = list(size = pt.size, width=0.5), 
                                                           text=as.formula(paste0('~',colnames(results)[3])), 
                                                           hoverinfo="text", plot_bgcolor = 'black', ...)))
    
    print(p)
    
  }
  
}
