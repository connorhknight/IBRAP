#' @title Method override for dim function
#'
#' @export

setMethod(f = 'dim', 
          signature = 'IBRAP',
          function(x) {
            dim(x@methods[[x@active.method]]@counts)
          })