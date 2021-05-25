#' @title Method override for dim function

setMethod(f = 'dim', 
          signature = 'IBRAP',
          function(x) {
            dim(x@methods[[1]]@counts)
          })