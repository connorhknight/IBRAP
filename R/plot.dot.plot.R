#' @name plot.dot.plot
#' @aliases plot.dot.plot
#' 
#' @title Plots a dot plot of gene expression
#'
#' @description Creates a dot plot of gene expression
#'
#' @param object An object of IBRAP class
#' @param assay Character. Which assay to access
#' @param slot Character. Which expression matrix slot should be used. Default = 'normalised'
#' @param clust.method Character. Which cluster method should be used, metadata accesses objects metadata 
#' @param column Character. Which column to access in the defined clust.method
#' @param features Character. Which features to plot
#'
#' @export plot.dot.plot

plot.dot.plot <- function(object, assay, slot='normalised', clust.method, column, features, ...) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string \n')
    
  } else if (is.character(assay)) {
    
    if(!assay %in% names(object@methods)) {
      
      stop('assay is not contained within object@methods \n')
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be character string \n')
    
  } else if(is.character(slot)) {
    
    if(!slot %in% c('counts','normalised','norm.scaled')) {
      
      stop('slot must be either counts, normalised or norm.scaled \n')
      
    }
    
  }
  
  if(!is.character(clust.method)) {
    
    stop('clust.method must be character string \n')
    
  } 
  
  if(!is.character(column)) {
    
    stop('column must be character string \n')
    
  }
  
  if(clust.method != 'metadata') {
    
    if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
      
      stop('clust.method not contained within cluster assignments \n')
      
    } else {
      
      if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
        
        stop('column not contained within clust.method \n')
        
      }
      
    }
    
  }
  
  for(x in features) {
    
    if(!x %in% rownames(object@methods[[assay]][[slot]])) {
      
      stop(paste0(x, ' not contained within expression matrix \n'))
      
    }
    
  }

  if(clust.method != 'metadata') {
    
    assignment <- object@methods[[assay]]@cluster_assignments[[clust.method]]
    
  } else if (clust.method == 'metadata') {
    
    assignment <- object@sample_metadata
    
  }
  
  sce <- SingleCellExperiment::SingleCellExperiment(list(expression=object@methods[[assay]][[slot]]),
                                                    colData=data.frame(cell_assignment=assignment[,column]))

  return(dittoDotPlot(object = sce, vars = features, group.by = 'cell_assignment', ...))
  
}
