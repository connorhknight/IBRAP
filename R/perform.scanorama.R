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
#' @param verbose Logical Should function messages be printed?
#' @param seed Numerical What seed should be set. Default = 1234
#' 
#' @return Scanorama reduction saved in the supplied method-assays
#' 
#' @examples 
#' 
#' object <- perform.scanorama(object = object, 
#'                             assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                             slot = 'norm.scaled', 
#'                             split.by = 'original.project', 
#'                             n.dims = 50)
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
                              verbose = FALSE,
                              seed=1234) {
  
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
    
  } else if(is.character(split.by)) {
    
    for(x in split.by) {
      
      if(!x %in% names(object@sample_metadata)) {
        
        stop(paste0(x, ' is not contained within object@sample_metadata'))
        
      }
      
    }
    
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
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  } 
  
  if(!is.numeric(seed)) {
    
    stop('seed must be a numerical value \n')
    
  }
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': loading python modules\n')))
    
  }
  
  scanorama <- reticulate::import('scanorama', convert = FALSE)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': python modules loaded\n')))
    
  }
  
  temp <- function(x) {
    
    return(paste(x, collapse = '_'))
    
  }
  
  count <- 1
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  
  if(!'integration_method' %in% colnames(object@pipelines)) {
    
    tmp <- tibble::add_column(.data = object@pipelines, integration_method=NA, integration_time=NA)
    
  } else {
    
    tmp <- object@pipelines
    
  }
  
  for(p in assay) {
    
    start_time <- Sys.time()
    
    if(length(split.by) > 1) {
      
      df <- object@sample_metadata[,split.by]
      
      df <- as.data.frame(apply(X = df, MARGIN = 1, FUN = temp))
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': initialising scanorama for assay: ', p, '\n')))
        
      }
      
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
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': matrices isolated\n')))
        
      }
      
      gene.list <- list()
      
      for(x in 1:length(sep)) {
        gene.list[[x]] <- rownames(mat[,df[,1] == x])
      }
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': genes identified\n')))
        cat(crayon::cyan(paste0(Sys.time(), ': corrections starting\n')))
        
      }
      
    } else {
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': initialising scanorama for assay: ', p, '\n')))
        
      }
      
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
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': matrices isolated\n')))
        
      }
      
      gene.list <- list()
      
      for(x in 1:length(sep)) {
        gene.list[[x]] <- rownames(mat[,object@sample_metadata[,split.by] == x])
      }
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': genes identified\n')))
        cat(crayon::cyan(paste0(Sys.time(), ': corrections starting\n')))
        
      }
      
    }
    
    if(isTRUE(verbose)) {
      
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
      
    }
    
    if(isFALSE(verbose)) {
      
      integrated.corrected.data <- scanorama$correct(datasets_full = reticulate::r_to_py(list.matrix), 
                                                     genes_list = reticulate::r_to_py(gene.list), 
                                                     dimred = as.integer(n.dims), 
                                                     return_dimred=TRUE, 
                                                     return_dense=FALSE, 
                                                     verbose = FALSE, 
                                                     batch_size = as.integer(batch_size), 
                                                     approx = approx, 
                                                     sigma = as.integer(sigma), 
                                                     alpha = as.numeric(alpha), 
                                                     knn = as.integer(knn),
                                                     union = as.logical(union),
                                                     seed = as.integer(seed))
      
    }
    
    dims <- list()
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': isolating scanorama reduced dimensions\n')))
      
    }
    
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
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': combining samples\n')))
      
    }
    
    combined <- do.call('cbind', dims)
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': samples concatenated\n')))
      
    }
    
    object@methods[[p]]@integration_reductions[[paste0('scanorama', reduction.save.suffix)]] <- t(combined)
    
    end_time <- Sys.time()
    
    function_time <- end_time - start_time
    
    if(!'integration_method' %in% colnames(object@pipelines)) {
      
      tmp[which(x = tmp$normalisation_method==p),'integration_method'] <- paste0('SCANORAMA', reduction.save.suffix)
      
      tmp[which(x = tmp$normalisation_method==p),'integration_time'] <- as.difftime(function_time, units = 'secs')
      
    }
    
    if('integration_method' %in% colnames(object@pipelines)) {
      
      if(paste0('SCANORAMA', reduction.save.suffix) %in% tmp$integration_method) {
        
        tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('SCANORAMA', reduction.save.suffix)),] <- c(tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('SCANORAMA', reduction.save.suffix)),c('normalisation_method','normalisation_time')], paste0('SCANORAMA', reduction.save.suffix), as.difftime(function_time, units = 'secs'))  
        
      }
      
      if(!paste0('SCANORAMA', reduction.save.suffix) %in% object@pipelines$integration_method) {
        
        df <- tmp[which(tmp$normalisation_method==p),]
        
        df <- df[!duplicated(df$normalisation_method),]
        
        df[,'integration_method'] <- paste0('SCANORAMA', reduction.save.suffix)
        
        df[,'integration_time'] <- function_time
        
        tmp <- rbind(tmp, df)
        
      }
      
    }
    
  }
  
  if(!'integration_method' %in% colnames(object@pipelines)) {
    
    tmp$integration_time <- as.difftime(tim = tmp$integration_time, units = 'secs')
    
    rownames(tmp) <- 1:nrow(tmp)
    
    object@pipelines <- tmp
    
  } else if ('integration_method' %in% colnames(object@pipelines)) {
    
    tmp$integration_time <- as.difftime(tim = tmp$integration_time, units = 'secs')
    
    rownames(tmp) <- 1:nrow(tmp)
    
    object@pipelines <- tmp
    
  }
  
  return(object)
  
}