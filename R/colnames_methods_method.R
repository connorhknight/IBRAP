#' @name colnames_methods
#' 
#' @title Method override for colnames function
#' 
#' @exportMethod 

setMethod(f = 'colnames', 
          signature = 'methods', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@counts)
          })
