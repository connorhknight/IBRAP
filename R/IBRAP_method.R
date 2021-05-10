#' @title Method for IBRAP S4 object
#'
#' @export

setMethod(f = 'show', signature = 'IBRAP', definition = function(object) {
  cat(crayon::white(paste0('An object of class ', 
                           class(object), 
                           '\n')))
  
  cat(crayon::white(paste0('  ', 
                           nrow(object@methods[[object@active.method]]@counts), 
                           ' features by ', 
                           ncol(object@methods[[object@active.method]]@counts), 
                           ' samples\n')))
  
  cat(crayon::white(paste0('  ', 'Active method:', 
                           object@active.method, ' (features:', nrow(object@methods[[object@active.method]]@counts), 
                           ', samples:', paste0(ncol(object@methods[[object@active.method]]@counts),')'), '\n')))
  
  lol <- names(object@methods)[1]
  
  if(length(names(object@methods)) > 1) {
    
    for(x in names(object@methods)[2:length(names(object@methods))]) {
      
      lol <- paste0(lol, ', ', x)
      
    }
    
  }
  
  cat(crayon::white(paste0('  Available methods: ', lol, '\n')))
  
})