#' @name perform.tsne
#' @aliases perform.tsne
#' 
#' @title Performs t-SNE reduction
#'
#' @description Performs t-SNE reduction on defined method-assays and supplied reductions or graphs.
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the t-SNE algorithm. Default = NULL
#' @param reduction.name.suffix Character. What should be appended to the end of tsne as the reduction name.
#' @param n.dims Numerical. How many dimensions of the supplied reduction be used, Null equates to all dimensions. Default = NULL
#' @param n_components Numerical. How many tsne dimensions should be produced, if you are supplying graphs, only 2 dimensions can be produced. Default = 3
#' @param ... Numerical. Arguments to be passed to ProjectionBasedClustering::tSNE
#' 
#' @return UMAP reduction saved in the visualisation_reductions section in the supplied method-assays
#'
#' @export

perform.tsne <- function(object, 
                         assay,
                         reduction=NULL,
                         reduction.name.suffix='', 
                         n.dims=NULL,
                         n_components = 3, 
                         ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Must be IBRAP class object\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop(paste0('assay must be character string\n'))
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop('assay: ', x, ' does not exist \n')
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    stop(paste0('reduction must be character string\n'))
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        stop(paste0('reduction: ', x, ' does not exist\n'))
        
      }
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    stop(paste0('reduction.save must be character string\n'))
    
  }
  
  if(!is.numeric(n_components)) {
    
    stop(paste0('n_components must be numerical\n'))
    
  }
  
  if(!is.list(n.dims)) {
    
    stop('dimensions must be supplied in list format\n')
    
  }
  
  for(g in assay) {
    
    reduction.list <- list()
    
    red.names <- c(names(object@methods[[g]]@computational_reductions), 
                   names(object@methods[[g]]@integration_reductions),
                   names(object@methods[[g]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[g]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[g]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[g]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[g]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[g]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[g]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    count <- 1
    
    for(r in reduction) {
      
      red <- reduction.list[[r]]
      
      dim <- n.dims[[count]]
      
      cat(crayon::cyan(paste0(Sys.time(), ': processing ', r, 'for assay:', g,'\n')))
      
      if(dim != 0) {
        
        c <- ProjectionBasedClustering::tSNE(DataOrDistances = red[,1:dim], 
                                             OutputDimension = n_components, Iterations = 1000, 
                                             verbose = TRUE,...)$ProjectedPoints
        
      } else {
        
        c <- ProjectionBasedClustering::tSNE(DataOrDistances = red, 
                                             OutputDimension = n_components, Iterations = 1000, 
                                             verbose = TRUE, ...)$ProjectedPoints
        
      }
      
      cat(crayon::cyan(paste0(Sys.time(), ': t-SNE reduction completed\n')))
      
      dim.names <- list()
      
      for(t in 1:n_components) {
        
        dim.names[[t]] <- paste0('tsne_', t)
        
      }
      
      colnames(c) <- unlist(dim.names)
      
      rownames(c) <- colnames(object)
      
      if('_' %in% unlist(strsplit(x = reduction.name.suffix, split = ''))) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in reduction.name.suffix, replacing with - \n')))
        reduction.name.suffix <- sub(pattern = '_', replacement = '-', x = reduction.name.suffix)
        
      }
      
      object@methods[[g]]@visualisation_reductions[[paste0(r, '_tsne', reduction.name.suffix)]] <- c
      
      cat(crayon::cyan(paste0(Sys.time(), ': t-SNE data added\n')))
      
      count <- count + 1
      
    }
    
  }
  
  return(object)
  
}