#' @name perform.dbmap
#' @aliases perform.dbmap
#' 
#' @title Performs dbMAP reduction
#'
#' @description Performs dbMAP reduction on defined method-assays. Unscaled but normalised data is automatically subset according to the identified HVGs for that method-assay
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param n_components Numerical. How many components should be produced, NULL enables automatic determination. Default = NULL
#' @param n_neighbors Numerical. How many neighbours are expected per cell. Default = 15
#' @param reduction.save Character. What should this reduction be saved as in computation_reduction. Default = 'dbmap'
#' @param save.object Boolean. Should the Anndata object use in python be saved under alt_objects. Default = TRUE
#' 
#' @return dbMAP reductions contained within the computational_reduction list in the defined assays
#'
#' @export

perform.dbmap <- function(object, 
                          assay, 
                          slot='normalised',
                          n_components = NULL, 
                          n_neighbors = 15, 
                          reduction.save='dbmap',
                          save.object = TRUE) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s)\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan('assay: ', x, ' does not exist\n'))
      return(object)
      
    }
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist \n'))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    if(!is.null(n_components)) {
      
      cat(crayon::cyan('n_components must be numerical or NULL\n'))
      return(object)
      
    }
    
  }
  
  if(!is.numeric(n_neighbors)) {
    
    cat(crayon::cyan('n_neighbors must be numerical \n'))
    return(object)
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string(s)\n'))
    return(object)
    
  }
  
  if(!is.logical(save.object)) {
    
    cat(crayon::cyan('save.object must be logical: TRUE/FALSE \n'))
    return(object)
    
  }
  
  if(is.null(reticulate::import('scipy.sparse', convert = FALSE))) {
    
    cat(crayon::cyan('scipy.sparse is not installed\n'))
    return(object)
    
  }
  
  if(is.null(reticulate::import('dbmap', convert = FALSE))) {
    
    cat(crayon::cyan('dbmap is not installed\n'))
    return(object)
    
  }
  
  scipy.sparse <- reticulate::import('scipy.sparse', convert = FALSE)
  
  dbmap <- reticulate::import('dbmap', convert = FALSE)
  
  for(o in assay) {
    
    cat(crayon::cyan(paste0('calculating dbmap for assay: ', o,'\n')))
    
    cellnames <- colnames(object@methods[[o]][[slot]])
    
    data <- scipy.sparse$csr_matrix(reticulate::r_to_py(t(as.matrix(object@methods[[o]][[slot]]))[,object@methods[[o]]@highly.variable.genes]))
    
    if(!is.null(n_components)) {
      
      diff <- dbmap$diffusion$Diffusor(n_components = as.integer(n_components), n_neighbors = as.integer(n_neighbors),
                                       transitions = as.logical(F),
                                       norm = as.logical(F), ann_dist = as.character('cosine'),
                                       n_jobs = as.integer(10), kernel_use = as.character('simple'))$fit(data)
      
    } else {
      
      diff <- dbmap$diffusion$Diffusor(n_neighbors = as.integer(n_neighbors),
                                       transitions = as.logical(F),
                                       norm = as.logical(F), ann_dist = as.character('cosine'),
                                       n_jobs = as.integer(10), kernel_use = as.character('simple'))$fit(data)
      
    }
    
    dbmap_components <- reticulate::py_to_r(diff$transform(data))
    
    res <- diff$return_dict()
    
    rownames(dbmap_components) <- cellnames
    
    dim.names <- list()
    for(t in 1:length(colnames(dbmap_components))) {
      dim.names[[t]] <- paste0('dbmap_', t)
    }
    colnames(dbmap_components) <- unlist(dim.names)
    
    object@methods[[o]]@computational_reductions[[reduction.save]] <- as.matrix(dbmap_components)
    
    if(isTRUE(save.object)) {
      
      object@methods[[o]]@alt_objects[[reduction.save]] <- diff
      
    }
    
  }
  
  return(object)
  
}