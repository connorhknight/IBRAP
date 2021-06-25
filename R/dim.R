#' @name dim
#' 
#' @title Method override for dim function
#' 
#' @exportMethod dim

setMethod(f = 'dim', 
          signature = 'IBRAP',
          function(x) {
            dim(x@methods[[1]]@counts)
          })