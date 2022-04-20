#' @name prepare.for.dittoSeq
#' @aliases prepare.for.dittoSeq
#' 
#' @title Create a SCE object to use for creating plots in dittoSeq
#'
#' @param object An IBRAP S4 class object
#' @param reduction Character. Which reduction to use to display points
#' @param assay Character. Which assay within the object to access
#' @param clust.method. Character. Whcih cluster assignments should be added to the metadata. the dataframe contained in `object@sample_metadata` is alread contained in the object. NULL will mean only sample_metadata will be in the returned object
#' @param pt.size Numeric. What should the point size be
#' @param cells Numeric. Which cells should be subset for plotting, Default = NULL
#' 
#' @return A SCE object containing 
#'
#' @export prepare.for.dittoSeq

prepare.for.dittoSeq <- function(object, assay, slot='normalised', clust.method, column, reduction = NULL, ...) {
  
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
  
  if(!is.null(reduction)) {
    
    if(!reduction %in% names(c(object@methods[[assay]]@computational_reductions, 
                               object@methods[[assay]]@visualisation_reductions, 
                               object@methods[[assay]]@integration_reductions))) {
      
      stop('reduction does not exist\n')
      
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
    
  }

  reduction <- reduction.list[reduction]

  if(!is.null(clust.method)) {
    
    assignment <- cbind(object@sample_metadata, object@methods[[assay]]@cluster_assignments[[clust.method]])
    
  } else {
    
    assignment <- object@sample_metadata
    
  }
  
  if(!is.null(reduction)) {
    
    sce <- SingleCellExperiment::SingleCellExperiment(list(expression=object@methods[[assay]][[slot]]), 
                                                      reducedDims = reduction,
                                                      colData=assignment)
    
  } else {
    
    sce <- SingleCellExperiment::SingleCellExperiment(list(expression=object@methods[[assay]][[slot]]), 
                                                      colData=assignment)
    
  }

  return(sce)
  
}