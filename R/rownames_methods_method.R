#' @title Method override for rownames function

setMethod(f = 'rownames', 
          signature = 'methods', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            rownames(x@counts)
          })
