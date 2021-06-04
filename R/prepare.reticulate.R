#' @name prepare.reticulate
#' @aliases prepare.reticulate
#' 
#' @title Installs or identifies if python modules are installed
#'
#' @description This function checks if IBRAPs dependent  
#' 
#' @export

prepare.reticulate <- function() {
  
  if(isFALSE(py_module_available('scrublet'))){
    
    py_install('scrublet', pip = T)
    
  }
  
  if(isFALSE(py_module_available('scrublet'))) {
    
    cat(cyan('Scrublet is not installed, please try manually.\n'))
    
  } else {
    
    cat(cyan('scrublet installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(py_module_available('scanpy'))){
    
    py_install('scanpy', pip = T)
    
  }
  
  if(isFALSE(py_module_available('scanpy'))) {
    
    cat(cyan('Scanpy is not install, please try manually.\n'))
    
  } else {
    
    cat(cyan('scanpy installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(py_module_available('bbknn'))){
    
    py_install('bbknn', pip = T)
    
  }
  
  if(isFALSE(py_module_available('bbknn'))) {
    
    cat(cyan('BBKNN is not install, please try manually.\n'))
    
  } else {
    
    cat(cyan('BBKNN installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(py_module_available('scanorama'))){
    
    py_install('scanorama', pip = T)
    
  }
  
  if(isFALSE(py_module_available('scanorama'))) {
    
    cat(cyan('Scanorama is not install, please try manually.\n'))
    
  } else {
    
    cat(cyan('Scanorama installed.\n'))
    
  }
  
  ####################################################
  
  if(isFALSE(py_module_available('dbmap'))){
    
    py_install('nmslib', pip = T)
    py_install('dbmap', pip = T)
    
  }
  
  if(isFALSE(py_module_available('dbmap'))) {
    
    cat(cyan('dbMAP is not install, please try manually.\n'))
    
  } else {
    
    cat(cyan('dbMAP installed.\n'))
    
  }
  
}


