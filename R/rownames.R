#' @name rownames
#' 
#' @title Method override for rownames function
#' 
#' @exportMethod rownames

setMethod(f = 'rownames', signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            rownames(x@methods[[1]]@counts)
          })

setMethod(f = 'rownames', 
          signature = 'methods', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            rownames(x@counts)
          })