#' @name dollar
#' 
#' @title Method override for `'$'` subset function regarding IBRAP S4 object
#' 
#' @exportMethod `$`

setMethod(f = '$', signature = 'IBRAP',
          function(x, 
                   name){
            
            x@sample_metadata[[name]]
            
          })