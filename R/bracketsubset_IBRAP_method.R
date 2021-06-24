#' @name doublet_bracket_ibrap
#' 
#' @title Method override for subset method `'[['` for IBRAP object
#' 
#' @exportMethod 

setMethod(f = '[[', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   ...) {
            x@sample_metadata[[i, exact = TRUE]]
          })