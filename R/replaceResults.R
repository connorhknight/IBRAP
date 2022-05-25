#' @name replaceResults
#' @aliases replaceResults
#' 
#' @title Shows the contents in your IBRAP object
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String showing which assay to access
#' @param item.to.replace Character. The name of the item to replace in the object.
#' @param replacement.item Supply the item to replace the item with. 
#' 
#' @return Prints out the contents of the supplied assays
#' 
#' @examples 
#'
#' @export replaceResults

replaceResults <- function(object, assay, item.to.replace, replacement.item) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string(s) \n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('assay: ', x, 'does not exist\n'))
      
    }
    
  }
  
  if(!is.character(item.to.replace)) {
    
    stop('item must be character string(s) \n')
    
  }
  
  if(item.to.replace %in% c('counts','normalised')) {
    
    if(!is(replacement.item, 'dgCMatrix')) {
      
      replacement.item <- as(replacement.item, 'dgCMatrix')
      
    }
    
    object@methods[[assay]][[item.to.replace]] <- replacement.item
    
  } else if (item.to.replace == 'norm.scaled') {
    
    if(!is.matrix(replacement.item)) {
      
      replacement.item <- as_matrix(replacement.item)
      
    }
    
    object@methods[[assay]]@norm.scaled <- replacement.item
    
  } else if (item.to.replace=='highly.variable.genes') {
    
    if(!is.character(replacement.item)) {
      
      replacement.item <- as.character(replacement.item)
      
    }
    
    object@methods[[assay]]@highly.variable.genes <- replacement.item
    
  } else if (item.to.replace %in% names(object@methods[[assay]]@computational_reductions)) {
    
    if(!is.data.frame(replacement.item)) {
      
      replacement.item <- as.data.frame(replacement.item)
      
    }
    
    object@methods[[assay]]@computational_reductions[[item.to.replace]] <- replacement.item
    
  } else if (item.to.replace %in% names(object@methods[[assay]]@integration_reductions)) {
    
    if(!is.data.frame(replacement.item)) {
      
      replacement.item <- as.data.frame(replacement.item)
      
    }
    
    object@methods[[assay]]@integration_reductions[[item.to.replace]] <- replacement.item
    
  } else if (item.to.replace %in% names(object@methods[[assay]]@visualisation_reductions)) {
    
    if(!is.data.frame(replacement.item)) {
      
      replacement.item <- as.data.frame(replacement.item)
      
    }
    
    object@methods[[assay]]@visualisation_reductions[[item.to.replace]] <- replacement.item
    
  } else if (item.to.replace %in% names(object@methods[[assay]]@neighbours)) {
    
    if(is.list(replacement.item)) {
      
      count <- 1
      
      for (x in 1:length(replacement.item)) {
        
        if(!is(replacement.item[[x]], 'dgCMatrix')) {
          
          replacement.item[[x]] <- as(replacement.item[[x]], 'dgCMatrix')
          
        }
        
      }
      
    } else {
      
      stop('neighbourhood graphs should be supplied as a list with the graph inside named connectivities \n')
      
    }
    
    object@methods[[assay]]@neighbours[[item.to.replace]] <- item.to.replace
    
  } else if (item.to.replace %in% names(object@methods[[assay]]@cluster_assignments)) {
    
    if(!is.data.frame(replacement.item)) {
      
      replacement.item <- as.data.frame(replacement.item)
      
    }
    
    object@methods[[assay]]@cluster_assignments[[item.to.replace]] <- replacement.item
    
  } else if ('clustering' %in% names(object@methods[[assay]]@benchmark_results)) {
    
    if(item.to.replace %in% names(object@methods[[assay]]@benchmark_results$clustering)) {
      
      if(!is.data.frame(replacement.item)) {
        
        replacement.item <- as.data.frame(replacement.item)
        
      }
      
      object@methods[[assay]]@benchmark_results$clustering[[item.to.replace]] <- replacement.item
      
    }

  }
  
  else if ('intergration' %in% names(object@methods[[assay]]@benchmark_results)) {
    
    if(!is.data.frame(replacement.item)) {
      
      replacement.item <- as.data.frame(replacement.item)
      
    }
    
    object@methods[[assay]]@benchmark_results$integration[[item.to.replace]] <- replacement.item
    
  } else {
    
    stop('item.to.replace cannot be found \n')
    
  }
  
  return(object)
  
}