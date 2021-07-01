#' @name perform.scanpy.neighbours
#'
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm.
#' @param neighbour.name Character. String defining the names to store the neighbourhood graph results under.
#' @param n_neighbors Numerical. How many neighbours should be found per cell, a higher value typically achieves more accurate results. Default = 15
#' @param dims Numerical. How many components of the reduction should be used, 0 means that all will be used. Default = 0
#' @param random_state Numerical. The seed value to use. Default = 0
#' @param method Character. String indicating which methodology to use including: ‘umap’, ‘gauss’ or ‘rapids’
#' @param metric Character. String indicating which distance metric to use, including: ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’, ‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’ or ‘manhattan’
#' @param generate.diffmap Boolean. Should diffusion maps be generated from the neighourhood graphs, these will be stored in computational_reductions and can be used for umap generation and further neighbourhood generation. Default = TRUE
#' @param n_comps Numerical. How many components should be generated for the diffusion maps. Default = 15
#' @param diffmap.name Character. What should the diffusion maps be named. 
#'
#' @export

perform.scanpy.neighbours <- function(object, 
                                      assay,
                                      reduction,
                                      neighbour.name,
                                      n_neighbors = 15, 
                                      dims = 0, 
                                      random_state = 0, 
                                      method = 'umap', 
                                      metric='euclidean',
                                      generate.diffmap = FALSE,
                                      n_comps = 15,
                                      diffmap.name = NULL
                                      ) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0(x, ' not in object@methods\n'))
      
    }
    
  }
  
  
  if(length(neighbour.name) != length(reduction)) {
    
    stop('n of neighbour.name is not equal to the n of supplied reductions \n')
    
  }
  
  if(!is.numeric(n_neighbors)) {
    
    stop('n_neighbors must be numerical \n')
    
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
  
  if(!is.numeric(random_state)) {
    
    stop('random_state must be numerical \n')
    
  }
  
  if(!is.character(method)) {
    
    stop('method must be character string \n')
    
  } else if (is.character(method)) {
    
    if(!method %in% c('umap', 'gauss', 'rapids')) {
      
      stop('method must be either ‘umap’, ‘gauss’ or ‘rapids’ \n')
      
    }
    
  }
  
  if(!is.character(metric)) {
    
    stop('metric must be character string \n')
    
  } else if (is.character(metric)) {
    
    if(!metric %in% c('braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule', 'cityblock', 'cosine', 'euclidean' , 'l1', 'l2', 'manhattan')) {
      
      stop('metric must be either ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’, ‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’ or ‘manhattan’ \n')
      
    }
    
  }
  
  if(!is.logical(generate.diffmap)) {
    
    stop('generate.diffmap must be boolean, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(n_comps)) {
    
    stop('n_comps must be numerical \n')
    
  }
  
  if(!is.null(diffmap.name)) {
    
    if(!is.character(diffmap.name)) {
      
      stop('diffmap.name must be character \n')
      
    } else if (is.character(diffmap.name)) {
      
      if(length(diffmap.name) != length(neighbour.name)) {
        
        stop('n of diffmap.name must be equal to n of neighbour.name \n')
        
      }
      
    }
    
  }
  
  sc <- reticulate::import('scanpy')
  pd <- reticulate::import('pandas')
  
  cat(crayon::cyan(paste0(Sys.time(), ': importing python packages \n')))
  
  
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
    
    for(g in reduction) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': processing assay: ', p, ', reduction: ', g, '\n')))
      
      scobj <- sc$AnnData(X = t(as.matrix(object@methods[[2]]@norm.scaled)))
      scobj$obs_names <- as.factor(colnames(object@methods[[2]]@norm.scaled))
      scobj$var_names <- as.factor(rownames(object@methods[[2]]@norm.scaled))
      
      if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
        
        scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))
        
      }
      
      red <- reduction.list[[g]]
      red.key <- reduction[count]
      dim <-dims[[count]]
      
      scobj$obsm$update(X_pca = as.matrix(red))
      
      sc$pp$neighbors(adata = scobj, 
                      n_neighbors = as.integer(n_neighbors), 
                      n_pcs = as.integer(dim),
                      random_state = as.integer(random_state), 
                      method = as.character(method), 
                      metric = as.character(metric))
      
      if(isTRUE(generate.diffmap)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': calcualting diffusion map\n')))
        
        sc$tl$diffmap(adata = scobj, n_comps = as.integer(n_comps))
        
        cat(crayon::cyan(paste0(Sys.time(), ': diffusion map calculated\n')))
        
        diffmap <- as.matrix(scobj$obsm[['X_diffmap']])

        DC_names <- list()
        
        counter <- 1
        
        for(x in 1:ncol(diffmap)) {
          
          DC_names[[counter]] <- paste0('DC_', counter)
          
          counter <- counter + 1
          
        }
        
        DC_names <- unlist(DC_names)
        
        colnames(diffmap) <- DC_names
        rownames(diffmap) <- colnames(object)
        
      }

      cat(crayon::cyan(paste0(Sys.time(), ': neighbours calculated\n')))
      
      conn <- as(as.matrix(scobj$obsp[['connectivities']]), 'dgCMatrix')
      colnames(conn) <- colnames(object)
      rownames(conn) <- colnames(object)
      
      object@methods[[p]]@neighbours[[neighbour.name[count]]][['connectivities']] <- conn
      
      dis <- as(as.matrix(scobj$obsp[['distances']]), 'dgCMatrix')
      colnames(dis) <- colnames(object)
      rownames(dis) <- colnames(object)
      
      object@methods[[p]]@neighbours[[neighbour.name[count]]][['distances']] <- dis
      
      cat(crayon::cyan(paste0(Sys.time(), ': results saved as ', neighbour.name[count], '\n')))
      
      if(isTRUE(generate.diffmap)) {
        
        object@methods[[p]]@computational_reductions[[diffmap.name[count]]] <- diffmap
        
        cat(crayon::cyan(paste0(Sys.time(), ': diffmap saved as ', diffmap.name[count], '\n')))
        
      }
      
      count <- count + 1
      
    }
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': completed \n')))
  
  return(object)
  
}
