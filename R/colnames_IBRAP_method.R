#' @title Method override for colnames function
#'
#' @export

setMethod(f = 'colnames', 
          signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@methods[[x@active.method]]@counts)
          })