#' @title Method override for rownames function

setMethod(f = 'rownames', signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            rownames(x@methods[[1]]@counts)
          })