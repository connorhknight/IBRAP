#' @title Method override for ncol function

setMethod(f = 'ncol', 
          signature = 'methods', 
          function(x) {
            ncol(x@counts)
          })
