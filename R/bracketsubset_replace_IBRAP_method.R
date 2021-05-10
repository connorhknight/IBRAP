

setMethod(f = '[[<-', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   value) {
            x@sample_metadata[[i]] <- value
            return(x)
          })