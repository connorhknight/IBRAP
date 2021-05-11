#' @title Method override for subset method `'[['` for IBRAP object

setMethod(f = '[[', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   ...) {
            x@sample_metadata[[i, exact = TRUE]]
          })