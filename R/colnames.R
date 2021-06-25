#' @name colnames
#' 
#' @exportMethod colnames

setMethod(f = 'colnames', 
          signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@methods[[1]]@counts)
          })


setMethod(f = 'colnames', 
          signature = 'methods', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@counts)
          })