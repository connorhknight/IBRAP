#' @name perform.harmony
#' @aliases perform.harmony
#' 
#' @title Performs Harmony integration
#'
#' @description Performs Harmony integration on defined method-assays and reductions contained within. This is performed on reductions. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param vars.use Character. A string of the column nmae that contains variables to regress. 
#' @param reduction Character. String defining the name of the reduction to provide for BBKNN. Default = 'pca'
#' @param reduction.save.suffix Character. What should be appended to the end of harmony as the reduction name
#' @param dims.use Numerical. Number of dimensions of the provided reduction to input into harmony, NULL equates to all dimensions. Default = NULL
#' @param ... Arguments to be passed to harmony::HarmonyMatrix
#' 
#' @return PCA reductions contained within the computational_reduction list in the defined assays
#' 
#' @examples
#' 
#' object <- perform.harmony(object = object, 
#'                           assay = c('SCRAN', 'SCT', 'SCANPY'), 
#'                           vars.use = 'original.project', 
#'                           reduction = c('pca'),  
#'                           max.iter.harmony = 100,
#'                           dims.use = list(NULL))
#'
#' @export

perform.harmony <- function(object, 
                            assay, 
                            vars.use, 
                            reduction = 'pca', 
                            reduction.save.suffix = '',
                            dims.use = NULL, 
                            ...) {
  
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
  
  if(!is.character(vars.use)) {
    
    stop('vars.use must be character string(s)\n')
    
  }
  
  if(!vars.use %in% names(object@sample_metadata)) {
    
    stop('vars.use does not exist\n')
    
  }
  
  if(!is.character(reduction)) {
    
    stop('reduction must be character string(s)\n')
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        stop(paste0('reduction: ', x, ' does not exist\n'))
        
      }
      
    }
    
  }
  
  if(!is.character(reduction.save.suffix)) {
    
    stop('reduction.save.suffix must be character string(s)\n')
    
  }
  
  if(is.list(dims.use)) {
    
    for(x in dims.use) {
      
      if(!is.null(x)) {
        
        if(!is.numeric(x)) {
          
          stop('dims must either be numeric or NULL\n')
          
        } 
        
      } else if (!is.numeric(x)) {
        
        if(!is.null(x)) {
          
          stop('dims must either be numeric or NULL\n')
          
        }
        
      }
      
    }
    
  }
  
  for(p in assay) {
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[p]]@computational_reductions), 
                   names(object@methods[[p]]@integration_reductions),
                   names(object@methods[[p]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[p]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        stop('reductions could not be found\n')
        
      }
      
    }
    
    count <- 1
    
    for(g in reduction) {
      
      red <- reduction.list[[g]]
      
      dims <- dims.use[[count]]
      
      if(is.null(dims)) {
        
        dims <- 1:ncol(red)
        
      }
      
      if('_' %in% unlist(strsplit(x = reduction.save.suffix, split = ''))) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in reduction.save.suffix, replacing with - \n')))
        reduction.save.suffix <- sub(pattern = '_', replacement = '-', x = reduction.save.suffix)
        
      }
      
      cat(crayon::cyan(paste0(Sys.time(), ': initialising harmony for assay: ', p, ', reduction: ', g, '\n')))
      
      harm <- harmony::HarmonyMatrix(data_mat = red[,dims], meta_data = object@sample_metadata, vars_use = vars.use, do_pca = FALSE, verbose = TRUE, ...)
      
      object@methods[[p]]@integration_reductions[[paste0(r, '_harmony', reduction.save.suffix)]] <- harm
      
      count <- count + 1
      
    }
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': harmony completed\n')))
  return(object)
}