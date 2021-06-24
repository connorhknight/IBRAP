#' @name nrow_methods
#' 
#' @title Method override for nrow function
#' 
#' @exportMethod 

setMethod(f = 'nrow', 
          signature = 'methods', 
          function(x) {
            nrow(x@counts)
          })
