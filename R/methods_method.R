#' @title Method for methods S4 object
#'
#' @export

setMethod(f = 'show', signature = 'methods', definition = function(object) {
  cat(crayon::white(paste0(object@method.name, ', an IBRAP method\n')))
  cat(crayon::white(paste0('  Contains: ', 
                           nrow(object@counts), 
                           ' features by ', 
                           ncol(object@counts), 
                           ' samples\n')))
})