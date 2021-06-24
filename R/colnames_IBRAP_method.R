#' @name colnames_ibrap 
#' 
#' @title Method override for colnames function
#' 
#' @exportMethod 

setMethod(f = 'colnames', 
          signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@methods[[1]]@counts)
          })
