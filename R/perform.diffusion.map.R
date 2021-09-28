#' @name perform.diffusion.map
#' @aliases perform.diffusion.map
#' 
#' @title Diffusion maps
#' 
#' @description Produces diffusion maps from previous reductions, i.e. PCA. Diffusion maps are known to better represent cellular trajectories in non-linear space 
#'
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm.
#' @param dims Numerical list. The number of dimensions to use for each reduction. This is supplied as a list respective to the order of reductions. 
#' @param n.dcs Numerical. The number of diffusion components to produce. Default = 15
#' @param k Numerical. How many neighbours should be found per cell. A higher value tends to be more accurate. Default = 15
#' @param diffmap.name.suffix Character. What should be used as a suffix for diffmap
#' 
#' @examples 
#' 
#' samp <- perform.diffusion.map(object = samp, 
#'                               assay = c('SCT','SCRAN','SCANPY'), 
#'                               reduction = 'pca', 
#'                               dims = list(20))                     
#'
#' @export

perform.diffusion.map <- function(object, 
                                  assay, 
                                  reduction, 
                                  dims, 
                                  n.dcs = 15,
                                  k = 15, 
                                  diffmap.name.suffix='',
                                  ...) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0(x, ' not in object@methods\n'))
      
    }
    
  }
  
  if(!is.list(dims)) {
    
    stop('dims must be supplied in list format \n')
    
  } else if(is.list(dims)) {
    
    for(x in dims) {
      
      if(!is.numeric(x)) {
        
        stop('dims items must be numerical \n')
        
      }
      
    }
    
  }
  
  if(!is.numeric(n.dcs)) {
    
    stop('n.dcs must be numerical \n')
    
  }
  
  if(!is.numeric(n.dcs)) {
    
    stop('n.dcs must be numerical \n')
    
  }
  
  if(!is.numeric(k)) {
    
    stop('k must be numerical \n')
    
  }
  
  if(!is.character(diffmap.name.suffix)) {
    
    stop('diffmap name msut be numerical \n')
    
  }
  
  if('_' %in% unlist(strsplit(x = diffmap.name.suffix, split = ''))) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in diffmap.name.suffix, replacing with - \n')))
    diffmap.name.suffix <- sub(pattern = '_', replacement = '-', x = diffmap.name.suffix)
    
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
    
    if(!is.null(reduction)) {
      
      for(r in reduction) {
        
        if(!r %in% names(reduction.list)) {
          
          stop('reductions could not be found\n')
          
        }
        
      }
      
    }
    
    count <- 1
    
    for(r in reduction) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': calculating diffusion map for assay: ', p, ' reduction: ', r, '\n')))
      
      dim <- dims[[count]]
      
      if(dim == 0) {
        
        dim <- ncol(reduction.list[[r]])
        
      }
      
      object@methods[[p]]@computational_reductions[[paste0(r, ':diffmap', diffmap.name.suffix)]] <- destiny::DiffusionMap(data = reduction.list[[r]][,1:dim], k = k, n_eigs = n.dcs, verbose = T, ...)@eigenvectors
      
    }
    
    count <- count + 1
    
  }
  
  return(object)
  
}

