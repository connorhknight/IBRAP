#' @title Method override for `'[['` subset function regarding methods S4 object
#'
#' @export

setMethod(f = '[[', signature = 'methods',
          function(x, 
                   i, 
                   j, 
                   ...) {
            
            y <- as.list(x)
            y[[i]]
            
          })