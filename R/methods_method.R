#' @name showmethod_methods
#' 
#' @title Method for methods S4 object
#' 
#' @exportMethod 

setMethod(f = 'show', signature = 'methods', definition = function(object) {
  cat(crayon::white(paste0('  Contains: ', 
                           nrow(object@counts), 
                           ' features by ', 
                           ncol(object@counts), 
                           ' samples\n')))
})