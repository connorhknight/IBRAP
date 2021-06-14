#' @title Method override for colnames function

setMethod(f = 'colnames', 
          signature = 'methods', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@counts)
          })
