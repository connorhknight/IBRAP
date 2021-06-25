#' @name show
#' 
#' @exportMethod show

setMethod(f = 'show', signature = 'IBRAP', definition = function(object) {
  cat(crayon::white(paste0('An object of class ', 
                           class(object), 
                           '\n')))
  
  cat(crayon::white(paste0('  ', 
                           nrow(object@methods[[1]]@counts), 
                           ' features by ', 
                           ncol(object@methods[[1]]@counts), 
                           ' samples\n')))
  
  lol <- names(object@methods)[1]
  
  if(length(names(object@methods)) > 1) {
    
    for(x in names(object@methods)[2:length(names(object@methods))]) {
      
      lol <- paste0(lol, ', ', x)
      
    }
    
  }
  
  cat(crayon::white(paste0('  Available methods: ', lol, '\n')))
  
})

setMethod(f = 'show', signature = 'methods', definition = function(object) {
  cat(crayon::white(paste0('  Contains: ', 
                           nrow(object@counts), 
                           ' features by ', 
                           ncol(object@counts), 
                           ' samples\n')))
})