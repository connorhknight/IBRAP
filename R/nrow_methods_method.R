#' @title Method override for nrow function

setMethod(f = 'nrow', 
          signature = 'methods', 
          function(x) {
            nrow(x@counts)
          })
