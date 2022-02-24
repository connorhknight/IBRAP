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
  
  for(x in assay) {
    
    cat(crayon::cyan( x, ' contains:\n'))
    cat(crayon::cyan('counts\n'))
    cat(crayon::cyan('normalised\n'))
    cat(crayon::cyan('norm.scaled\n'))
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
    cat(crayon::cyan('benchmarking_results: '))
    cat(crayon::cyan(paste0(names(object@methods[[x]]@benchmark_results), collapse = ', '), '\n'))
    cat(crayon::cyan('          \n'))
    
  }

}
