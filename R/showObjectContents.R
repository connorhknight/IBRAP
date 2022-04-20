#' @name showObjectContents
#' @aliases showObjectContents
#' 
#' @title Shows the contents in your IBRAP object
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing which assays to show
#' 
#' @return Prints out the contents of the supplied assays
#' 
#' @examples 
#'
#' @export showObjectContents

showObjectContents <- function(object, assay) {
  
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
  
  for(x in assay) {
    cat(crayon::cyan( x, ' contains:\n'))
    cat(crayon::cyan('counts\n'))
    cat(crayon::cyan('normalised\n'))
    cat(crayon::cyan('norm.scaled\n'))
    if(!is.null(object@methods[[x]]@highly.variable.genes)) {
      cat(crayon::cyan('HVGs: TRUE \n')) 
    } else {
      cat(crayon::cyan('HVGs: FALSE \n')) 
    }
    cat(crayon::cyan('computational_reductions: ')) 
    cat(crayon::cyan(paste0(names(object@methods[[x]]@computational_reductions), collapse = ', '), '\n'))
    cat(crayon::cyan('integration_reductions: '))
    cat(crayon::cyan(paste0(names(object@methods[[x]]@integration_reductions), collapse = ', '), '\n'))
    cat(crayon::cyan('visualisation_reductions: '))
    cat(crayon::cyan(paste0(names(object@methods[[x]]@visualisation_reductions), collapse = ', '), '\n'))
    cat(crayon::cyan('neighbourhood: '))
    cat(crayon::cyan(paste0(names(object@methods[[x]]@neighbours), collapse = ', '), '\n'))
    cat(crayon::cyan('cluster_assignments: '))
    cat(crayon::cyan(paste0(names(object@methods[[x]]@cluster_assignments), collapse = ', '), '\n'))
    cat(crayon::cyan('benchmarking_results: \n'))
    if('clustering' %in% names(object@methods[[x]]@benchmark_results)) {
      cat(crayon::cyan('(clustering) '))
      cat(crayon::cyan(paste0(names(object@methods[[x]]@benchmark_results$clustering), collapse = ', '), '\n'))
    }
    if('integration' %in% names(object@methods[[x]]@benchmark_results)) {
      cat(crayon::cyan('(integration) '))
      cat(crayon::cyan(paste0(names(object@methods[[x]]@benchmark_results$integration), collapse = ', '), '\n'))
    }
    cat(crayon::cyan('          \n'))
  }
}
