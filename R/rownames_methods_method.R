#' @title Method override for rownames function
#' 
#' @exportMethod 

setMethod(f = 'rownames', 
          signature = 'methods', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            rownames(x@counts)
          })
