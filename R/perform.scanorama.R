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
#' @param reduction.save.suffix Character. Should a suffix be added to the end of scanorama, This cannot include underscores.
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
                              reduction.save.suffix='', 
                              batch_size = 5000, 
                              approx = TRUE, 
                              sigma = 15, 
                              alpha = 0.1, 
                              knn = 20,
                              union = FALSE,
                              seed = 12345) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('reduction: ', x, 'does not exist\n'))
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be a character string\n')
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    stop('slot does not exist\n')
    
  }
  
  if(!is.character(split.by)) {
    
    stop('split.by must be character string\n')
    
  }
  
  if(!is.numeric(n.dims)) {
    
    stop('n.dims must be numerical\n')
    
  }
  
  if(!is.character(reduction.save.suffix)) {
    
    stop('reduction.save.suffix must be character string\n')
    
  }
  
  if(!is.numeric(batch_size)) {
    
    stop('batch_size must be numerical\n')
    
  }
  
  if(!is.logical(approx)) {
    
    stop('approx must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.numeric(sigma)) {
    
    stop('sigma must be numerical\n')
    
  }
  
  if(!is.numeric(alpha)) {
    
    stop('alpha must be numerical\n')
    
  }
  
  if(!is.numeric(knn)) {
    
    stop('knn must be numerical\n')
    
  }
  
  if(!is.logical(union)) {
    
    stop('union must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed must be numerical\n')
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': loading python modules\n')))
  scanorama <- reticulate::import('scanorama', convert = FALSE)
  cat(crayon::cyan(paste0(Sys.time(), ': python modules loaded\n')))
  
  temp <- function(x) {
    
    return(paste(x, collapse = '_'))
    
  }
  
  count <- 1
  
  for(p in assay) {
    
    if(length(split.by) > 1) {
      
      df <- object@sample_metadata[,split.by]
      df <- as.data.frame(apply(X = df, MARGIN = 1, FUN = temp))
      
      cat(crayon::cyan(paste0(Sys.time(), ': initialising scanorama for assay: ', p, '\n')))
      
      list.matrix <- list()
      column.names <- list()
      sep <- unique(df[,1])
      mat <- object@methods[[p]][[slot]]
      counter <- 1
      
      for(x in sep) {
        
        column.names[[counter]] <- colnames(mat[,df[,1] == x])
        list.matrix[[counter]] <- t(mat[,df[,1] == x])
        counter <- counter + 1
        
      }
      
      cat(crayon::cyan(paste0(Sys.time(), ': matrices isolated\n')))
      gene.list <- list()
      
      for(x in 1:length(sep)) {
        gene.list[[x]] <- rownames(mat[,df[,1] == x])
      }
      
      cat(crayon::cyan(paste0(Sys.time(), ': genes identified\n')))
      cat(crayon::cyan(paste0(Sys.time(), ': corrections starting\n')))
      
    } else {
      
      cat(crayon::cyan(paste0(Sys.time(), ': initialising scanorama for assay: ', p, '\n')))
      
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
      
      cat(crayon::cyan(paste0(Sys.time(), ': matrices isolated\n')))
      gene.list <- list()
      
      for(x in 1:length(sep)) {
        gene.list[[x]] <- rownames(mat[,object@sample_metadata[,split.by] == x])
      }
      
      cat(crayon::cyan(paste0(Sys.time(), ': genes identified\n')))
      cat(crayon::cyan(paste0(Sys.time(), ': corrections starting\n')))
      
    }
    
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
    cat(crayon::cyan(paste0(Sys.time(), ': isolating scanorama reduced dimensions\n')))
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
    
    if('_' %in% unlist(strsplit(x = reduction.save.suffix, split = ''))) {
      
      reduction.save.suffix <- sub(pattern = '_', replacement = '-', x = reduction.save.suffix)
      
    }
    
    cat(crayon::cyan(paste0(Sys.time(), ': combining samples\n')))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan(paste0(Sys.time(), ': samples concatenated\n')))
    object@methods[[p]]@integration_reductions[[paste0('scanorama', reduction.save.suffix)]] <- t(combined)
    
  }
  
  return(object)
  
}