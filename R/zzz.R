# Set up environment

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to IBRAP")
}

.onLoad <- function(libname, pkgname) {
  
  packageStartupMessage("Welcome to IBRAP")
  
  reticulate::py_install('scrublet==0.2.3', pip = T)
  
  if(isFALSE(reticulate::py_module_available('scrublet'))) {
    
    cat(crayon::cyan('scrublet is not installed, please try manually: scrublet\n'))
    
  } else {
    
    cat(crayon::cyan('scrublet installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('scanpy==1.7.2', pip=T)
  
  if(isFALSE(reticulate::py_module_available('scanpy'))) {
    
    cat(crayon::cyan('scanpy is not install, please try manually: scanpy==1.8.1\n'))
    
  } else {
    
    cat(crayon::cyan('scanpy installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('bbknn==1.5.1', pip = T)
  
  if(isFALSE(reticulate::py_module_available('bbknn'))) {
    
    cat(crayon::cyan('BBKNN is not install, please try manually: bbknn==1.5.1\n'))
    
  } else {
    
    cat(crayon::cyan('BBKNN installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('scanorama==1.7.1', pip = T)
  
  if(isFALSE(reticulate::py_module_available('scanorama'))) {
    
    cat(crayon::cyan('scanorama is not install, please try manually: scanorama==1.7.1\n'))
    
  } else {
    
    cat(crayon::cyan('scanorama installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('louvain', pip = T)
  
  if(isFALSE(reticulate::py_module_available('louvain'))) {
    
    cat(crayon::cyan('louvain is not install, please try manually: louvain\n'))
    
  } else {
    
    cat(crayon::cyan('louvain installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('leiden==1.0.1', pip = T)
  
  if(isFALSE(reticulate::py_module_available('leiden'))) {
    
    cat(crayon::cyan('leiden is not install, please try manually: leiden\n'))
    
  } else {
    
    cat(crayon::cyan('louvain installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('leidenalg==0.8.8', pip = T)
  
  if(isFALSE(reticulate::py_module_available('leidenalg'))) {
    
    cat(crayon::cyan('leidenalg is not install, please try manually: leidenalgn'))
    
  } else {
    
    cat(crayon::cyan('leidenalg installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('annoy==1.16.0', pip = T)
  
  if(isFALSE(reticulate::py_module_available('annoy'))) {
    
    cat(crayon::cyan('annoy is not install, please try manually: annoy==1.16.0\n'))
    
  } else {
    
    cat(crayon::cyan('annoy installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('umap-learn==0.5.1', pip = T)
  
  if(isFALSE(reticulate::py_module_available('umap'))) {
    
    cat(crayon::cyan('umap-learn is not install, please try manually: umap-learn==0.5.1\n'))
    
  } else {
    
    cat(crayon::cyan('umap-learn installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('pandas==1.1.5', pip = T)
  
  if(isFALSE(reticulate::py_module_available('pandas'))) {
    
    cat(crayon::cyan('pandas is not install, please try manually: pandas==1.1.5 \n'))
    
  } else {
    
    cat(crayon::cyan('pandas installed.\n'))
    
  }  
  
  ####################################################
  
  reticulate::py_install('numba==0.53.1', pip = T)
  
  if(isFALSE(reticulate::py_module_available('numba'))) {
    
    cat(crayon::cyan('numba is not install, please try manually: numba==0.53.1 \n'))
    
  } else {
    
    cat(crayon::cyan('numba installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('numpy==1.21.0', pip = T)
  
  if(isFALSE(reticulate::py_module_available('numpy'))) {
    
    cat(crayon::cyan('numpy is not install, please try manually: numpy==1.21.0 \n'))
    
  } else {
    
    cat(crayon::cyan('numpy installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('matplotlib==3.4.2', pip = T)
  
  if(isFALSE(reticulate::py_module_available('matplotlib'))) {
    
    cat(crayon::cyan('matplotlib is not install, please try manually: matplotlib==3.4.2 \n'))
    
  } else {
    
    cat(crayon::cyan('matplotlib installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('scipy==1.7.0', pip = T)
  
  if(isFALSE(reticulate::py_module_available('scipy'))) {
    
    cat(crayon::cyan('scipy is not install, please try manually: scipy==1.7.0 \n'))
    
  } else {
    
    cat(crayon::cyan('scipy installed.\n'))
    
  }
  
  ####################################################
  
  reticulate::py_install('scikit-learn==0.24.2', pip = T)
  
  if(isFALSE(reticulate::py_module_available('sklearn'))) {
    
    cat(crayon::cyan('scikit-learn is not install, please try manually: scikit-learn==0.24.2 \n'))
    
  } else {
    
    cat(crayon::cyan('scikit-learn installed.\n'))
    
  }
  
}