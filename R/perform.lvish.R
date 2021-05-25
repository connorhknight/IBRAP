#' @name perform.lvish
#' @aliases perform.lvish
#' 
#' @title Performs LargeVis reduction
#'
#' @description Performs LargeVis reduction on defined method-assays and supplied reductions or graphs.
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the LargeVis algorithm. Default = NULL
#' @param reduction.save Character. What should the reduction be saved as. Default = 'lvish'
#' @param n.dims Numerical. How many dimensions of the supplied reduction be used, Null equates to all dimensions. Default = NULL
#' @param n_components Numerical. How many tsne dimensions should be produced, if you are supplying graphs, only 2 dimensions can be produced. Default = 3
#' @param ... Numerical. Arguments to be passed to uwot::lvish
#' 
#' @return LargeVis reduction saved in the visualisation_reductions section in the supplied method-assays
#'
#' @export

perform.lvish <- function(object, 
                          assay,
                          reduction='pca',
                          reduction.save='lvish',
                          n.dim=NULL, 
                          n_components = 3, 
                          ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan('reduction must be character string\n'))
    return(object)
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    cat(crayon::cyan('n_components must be numerical\n'))
    return(object)
    
  }
  
  if(!is.list(n.dim)) {
    
    cat(crayon::cyan('dimensions must be supplied in list format\n'))
    return(object)
    
  }
  
  for(u in assay) {
    
    reduction.list <- list()
    
    red.names <- c(names(object@methods[[u]]@computational_reductions), 
                   names(object@methods[[u]]@integration_reductions),
                   names(object@methods[[u]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[u]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[u]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[u]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[u]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[u]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[u]]@visualisation_reductions[[i]]
        
      }
    }
    
    count <- 1
    
    for(r in reduction) {
      
      red <- reduction.list[[r]]
      
      red.save <- reduction.save[count]
      
      dim <- n.dim[[count]]
      
      cat(crayon::cyan('Processing', r, 'for assay:', u,'\n'))
      
      if(!is.null(dim)) {
        
        c <- uwot::lvish(X = red[,dim], n_components = n_components, verbose = TRUE, ...)
        
      } else {
        
        c <- uwot::lvish(X = red, n_components = n_components, verbose = TRUE, ...)
        
      }
      
      dim.names <- list()
      
      for(l in 1:n_components) {
        
        dim.names[[l]] <- paste0('lvish_', l)
        
      }
      
      colnames(c) <- unlist(dim.names)
      
      object@methods[[u]]@visualisation_reductions[[red.save]] <- c
      
    }
    
  }
  
  return(object)
  
}