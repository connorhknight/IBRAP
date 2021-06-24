#' @name double_bracket_replace_methods
#' 
#' @title Method override for `'[['` subset function regarding methods S4 object
#' 
#' @exportMethod 

setMethod(f = '[[', signature = 'methods',
          function(x, 
                   i, 
                   j, 
                   ...) {
            
            y <- as.list(x)
            y[[i]]
            
          })