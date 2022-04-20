fetchResults <- function(object, assay, item) {
  
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
  
  if(!is.character(item)) {
    
    stop('item must be character string(s) \n')
    
  }
  
  if(item %in% c('counts','normalised','norm.scaled')) {
    
    return(object@methods[[assay]][[item]])
    
  } else if (item=='highly.variable.genes') {
    
    return(object@methods[[assay]]@highly.variable.genes)
    
  } else if (item %in% names(object@methods[[assay]]@computational_reductions)) {
    
    return(object@methods[[assay]]@computational_reductions[[item]])
    
  } else if (item %in% names(object@methods[[assay]]@integration_reductions)) {
    
    return(object@methods[[assay]]@integration_reductions[[item]])
    
  } else if (item %in% names(object@methods[[assay]]@visualisation_reductions)) {
    
    return(object@methods[[assay]]@visualisation_reductions[[item]])
    
  } else if (item %in% names(object@methods[[assay]]@neighbours)) {
    
    return(object@methods[[assay]]@neighbours[[item]])
    
  } else if (item %in% names(object@methods[[assay]]@neighbours)) {
    
    return(object@methods[[assay]]@neighbours[[item]])
    
  } else if (item %in% names(object@methods[[assay]]@cluster_assignments)) {
    
    return(object@methods[[assay]]@cluster_assignments[[item]])
    
  } else if ('clustering' %in% names(object@methods[[assay]]@benchmark_results)) {
    
    return(object@methods[[assay]]@benchmark_results$clustering[[item]])
    
  }
  
  else if ('clustering' %in% names(object@methods[[assay]]@benchmark_results)) {
    
    return(object@methods[[assay]]@benchmark_results$clustering[[item]])
    
  } else {
    
    stop('item cannot be found \n')
    
  }
  
}