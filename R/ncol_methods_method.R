#' @name ncol_methods
#' 
#' @title Method override for ncol function
#' 
#' @exportMethod 

setMethod(f = 'ncol', 
          signature = 'methods', 
          function(x) {
            ncol(x@counts)
          })
