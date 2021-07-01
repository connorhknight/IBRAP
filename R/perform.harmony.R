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
#' @param reduction.save Character. What should the harmony reduction be saved as. Default = 'harmony'
#' @param dims.use Numerical. Number of dimensions of the provided reduction to input into harmony, NULL equates to all dimensions. Default = NULL
#' @param ... Arguments to be passed to PCAtools::pca
#' 
#' @return PCA reductions contained within the computational_reduction list in the defined assays
#'
#' @export

perform.harmony <- function(object, 
                            assay, 
                            vars.use, 
                            reduction = 'pca', 
                            reduction.save = 'harmony',
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
  
  if(!is.character(reduction.save)) {
    
    stop('reduction.save must be character string(s)\n')
    
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
      
      red.save <- reduction.save[[count]]
      
      dims <- dims.use[[count]]
      
      if(is.null(dims)) {
        
        dims <- 1:ncol(red)
        
      }
      
      cat(crayon::cyan(paste0(Sys.time(), ': initialising harmony for assay: ', p, ', reduction: ', g, '\n')))
      
      harm <- harmony::HarmonyMatrix(data_mat = red[,dims], meta_data = object@sample_metadata, vars_use = vars.use, do_pca = FALSE, verbose = TRUE, ...)
      
      object@methods[[p]]@integration_reductions[[red.save]] <- harm
      
      count <- count + 1
      
    }
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': harmony completed\n')))
  return(object)
}