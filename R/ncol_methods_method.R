#' @name ncol
#' 
#' @title Method override for ncol function
#' 
#' @exportMethod ncol

setMethod(f = 'ncol', 
          signature = 'methods', 
          function(x) {
            ncol(x@counts)
          })

setMethod(f = 'ncol', 
          signature = 'IBRAP', 
          function(x) {
            ncol(x@methods[[1]]@counts)
          })