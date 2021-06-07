#' @name perform.scanorama
#' @aliases perform.scanorama
#' 
#' @title Performs Scanorama integration
#'
#' @description Performs Scanorama integration on defined method-assays and reductions contained within. This is performed on reductions. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String defining which slot in the assay to supply to Scanorama. Default = NULL
#' @param split.by Character. indicating the metadata column containing the batch to split the assay by. 
#' @param n.dims Numerical. The number of Scanorama dimensions to be produced. Default = 50
#' @param reduction.save Character. What should the Scanorama reduction be saved as. Default = 'scanorama'
#' @param batch_size Numerical. The batch size used in the alignment vector computation. Useful when integrating large datasets. Default = 5000
#' @param approx Boolean. Use appoximate nearest neighbours within python, speeds up runtime. Default = TRUE
#' @param sigma Numerical. Correction smoothing parameter on Gaussian kernel. Default = 15
#' @param alpha Numerical. Alignment score minimum cutoff. Default = 0.1
#' @param knn Numerical. Number of nearest neighbors to use for matching. Default = 20
#' @param union Boolean. Should genes between datasets be intersected or not. Default = FALSE
#' @param seed Numerical. The seed to use when integrating these datasets. Default = 12345
#' 
#' @return Scanorama reduction saved in the supplied method-assays
#'
#' @export

perform.scanorama <- function(object, 
                              assay,
                              slot = 'norm.scaled',
                              split.by, 
                              n.dims = 50, 
                              reduction.save='scanorama', 
                              batch_size = 5000, 
                              approx = TRUE, 
                              sigma = 15, 
                              alpha = 0.1, 
                              knn = 20,
                              union = FALSE,
                              seed = 12345) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be a character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(split.by)) {
    
    cat(crayon::cyan('split.by must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(n.dims)) {
    
    cat(crayon::cyan('n.dims must be numerical\n'))
    return(object)
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(batch_size)) {
    
    cat(crayon::cyan('batch_size must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(approx)) {
    
    cat(crayon::cyan('approx must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.numeric(sigma)) {
    
    cat(crayon::cyan('sigma must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(alpha)) {
    
    cat(crayon::cyan('alpha must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(knn)) {
    
    cat(crayon::cyan('knn must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(union)) {
    
    cat(crayon::cyan('union must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.numeric(seed)) {
    
    cat(crayon::cyan('seed must be numerical\n'))
    return(object)
    
  }
  
  cat(crayon::cyan('Loading python modules\n'))
  scanorama <- reticulate::import('scanorama', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  
  count <- 1
  
  for(p in assay) {
    
    cat(crayon::cyan(paste0('Initialising scanorama for assay: ', p, '\n')))
    
    list.matrix <- list()
    column.names <- list()
    sep <- unique(object@sample_metadata[,split.by])
    mat <- object@methods[[p]][[slot]]
    counter <- 1
    
    for(x in sep) {
      column.names[[counter]] <- colnames(mat[,object@sample_metadata[,split.by] == x])
      list.matrix[[counter]] <- t(mat[,object@sample_metadata[,split.by] == x])
      counter <- counter + 1
    }
    
    cat(crayon::cyan('Matrices isolated\n'))
    gene.list <- list()
    
    for(x in 1:length(sep)) {
      gene.list[[x]] <- rownames(mat[,object@sample_metadata[,split.by] == x])
    }
    
    cat(crayon::cyan('Genes identified\n'))
    cat(crayon::cyan('Corrections starting\n'))
    integrated.corrected.data <- scanorama$correct(datasets_full = reticulate::r_to_py(list.matrix), 
                                                   genes_list = reticulate::r_to_py(gene.list), 
                                                   dimred = as.integer(n.dims), 
                                                   return_dimred=TRUE, 
                                                   return_dense=FALSE, 
                                                   verbose = TRUE, 
                                                   batch_size = as.integer(batch_size), 
                                                   approx = approx, 
                                                   sigma = as.integer(sigma), 
                                                   alpha = as.numeric(alpha), 
                                                   knn = as.integer(knn),
                                                   union = as.logical(union),
                                                   seed = as.integer(seed))
    
    dims <- list()
    cat(crayon::cyan('Isolating scanorama reduced dimensions\n'))
    dim.names <- list()
    
    for(c in 1:n.dims) {
      dim.names[[c]] <- paste0('scanorama_', c)
    }
    
    dim.names <- unlist(dim.names)
    
    for(x in 1:length(sep)) {
      transposed <- t(reticulate::py_to_r(integrated.corrected.data)[[1]][[x]])
      colnames(transposed) <- column.names[[x]]
      rownames(transposed) <- dim.names
      dims[[x]] <- transposed
    }
    
    cat(crayon::cyan('Combining samples\n'))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan('Samples concatenated\n'))
    object@methods[[p]]@integration_reductions[[reduction.save]] <- t(combined)
    
  }
  
  return(object)
  
}