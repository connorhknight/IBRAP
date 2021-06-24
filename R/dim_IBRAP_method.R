#' @name dim_ibrap
#' 
#' @title Method override for dim function
#' 
#' @exportMethod 

setMethod(f = 'dim', 
          signature = 'IBRAP',
          function(x) {
            dim(x@methods[[1]]@counts)
          })