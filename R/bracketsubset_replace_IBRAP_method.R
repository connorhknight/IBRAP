#' @name doublet_bracket_replace_ibrap
#' 
#' @title Method override for colnames function
#' 
#' @exportMethod 

setMethod(f = '[[<-', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   value) {
            x@sample_metadata[[i]] <- value
            return(x)
          })