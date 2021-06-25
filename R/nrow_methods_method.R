#' @name nrow
#' 
#' @title Method override for nrow function
#' 
#' @exportMethod nrow

setMethod(f = 'nrow', 
          signature = 'methods', 
          function(x) {
            nrow(x@counts)
          })

setMethod(f = 'nrow', 
          signature = 'IBRAP', 
          function(x) {
            nrow(x@methods[[1]]@counts)
          })