#' @name prepare.reticulate
#' @aliases prepare.reticulate
#' 
#' @title Installs or identifies if python modules are installed
#'
#' @description This function checks if IBRAPs dependent  
#' 
#' @export

prepare.reticulate <- function() {
  
  if(isFALSE(reticulate::py_module_available('scrublet'))){
    
    reticulate::py_install('scrublet', pip = T)
    
  }
  
  if(isFALSE(reticulate::py_module_available('scrublet'))) {
    
    cat(crayon::cyan('Scrublet is not installed, please try manually.\n'))
    
  } else {
    
    cat(crayon::cyan('scrublet installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(reticulate::py_module_available('scanpy'))){
    
    reticulate::py_install('scanpy', pip = T)
    
  }
  
  if(isFALSE(reticulate::py_module_available('scanpy'))) {
    
    cat(crayon::cyan('Scanpy is not install, please try manually.\n'))
    
  } else {
    
    cat(crayon::cyan('scanpy installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(reticulate::py_module_available('bbknn'))){
    
    reticulate::py_install('bbknn', pip = T)
    
  }
  
  if(isFALSE(reticulate::py_module_available('bbknn'))) {
    
    cat(crayon::cyan('BBKNN is not install, please try manually.\n'))
    
  } else {
    
    cat(crayon::cyan('BBKNN installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(reticulate::py_module_available('scanorama'))){
    
    reticulate::py_install('scanorama', pip = T)
    
  }
  
  if(isFALSE(reticulate::py_module_available('scanorama'))) {
    
    cat(crayon::cyan('Scanorama is not install, please try manually.\n'))
    
  } else {
    
    cat(crayon::cyan('Scanorama installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(reticulate::py_module_available('dbmap'))){
    
    reticulate::py_install('nmslib', pip = T)
    reticulate::py_install('dbmap', pip = T)
    
  }
  
  if(isFALSE(reticulate::py_module_available('dbmap'))) {
    
    cat(crayon::cyan('dbMAP is not install, please try manually.\n'))
    
  } else {
    
    cat(crayon::cyan('dbMAP installed.\n'))
    
  }
  
}


