#' @name perform.nn
#' @aliases perform.nn
#' 
#' @title Finds the shared nearest neighbourhood for the cells. This supplies a graph. 
#' 
#' @description Neighbourhood graph generator utilised. According to the assay name a different method will be used. If the assay is derived from scanpy then the scanpy algorithm will be applied. Otherwise, the seurat implementation will be applied. 
#'
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm.
#' @param neighbour.name.suffix Character. String defining the names to store the neighbourhood graph results under.
#' 
#' @param k.param Numerical. The number of k when calcuating k-nearest neighbour. Default = 20
#' @param compute.SNN Boolean. Should the shared nearest neighbour graph be calculated. Default = TRUE 
#' @param prune.SNN Numerical. Setas acceptance cutoff for jaccard index whilst computing neighbourhood overlap for SNN construction. Any edges with a value less than this parameter will be removed. 0 = no pruning and 1 = prune everything. Default = 0
#' @param nn.method Character. Nearest neighbour method, either 'rann' or 'annoy'. Default = 'annoy'
#' @param n.trees Numerical. More trees facilitates hgiher precision when using 'annoy' method. Default = 20
#' @param nn.eps Numerical. Margin of error when performing nearest neighbour search whilst using rann method. 0 would imply an exact search. Default = 0.0
#' @param annoy.metric Character. Distance metric for annoy method. Options: 'euclidean', 'cosine', 'manhattan', 'hamming'. Default = 'euclidean'
#'
#' @param n_neighbors Numerical. (scanpy only) How many neighbours should be found per cell, a higher value typically achieves more accurate results. Default = 15
#' @param dims.use Numerical. (scanpy only) How many components of the reduction should be used, 0 means that all will be used. Default = 0
#' @param random_state Numerical. (scanpy only) The seed value to use. Default = 0
#' @param method Character. (scanpy only) String indicating which methodology to use including: ‘umap’, ‘gauss’ or ‘rapids’
#' @param metric Character. (scanpy only) String indicating which distance metric to use, including: ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘correlation’, ‘dice’, ‘hamming’, ‘jaccard’, ‘kulsinski’, ‘mahalanobis’, ‘minkowski’, ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’, ‘sokalsneath’, ‘sqeuclidean’, ‘yule’, ‘cityblock’, ‘cosine’, ‘euclidean’, ‘l1’, ‘l2’ or ‘manhattan’
#' @param generate.diffmap Boolean. (scanpy only) Should diffusion maps be generated from the neighourhood graphs, these will be stored in computational_reductions and can be used for umap generation and further neighbourhood generation. Default = TRUE
#' @param n_comps Numerical. (scanpy only) How many components should be generated for the diffusion maps. Default = 15
#' @param diffmap.name Character. (scanpy only) What should the diffusion maps be named. 
#' 
#' @examples 
#' 
#' # generates a diffusion map from the scanpy assay
#' object <- perform.nn(object = object, assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                      reduction = c('pca'),
#'                      dims = list(0,0))
#'                         
#'
#' @export

perform.nn <- function(object, 
                       assay, 
                       reduction, 
                       
                       neighbour.name.suffix = '',
                       dims.use=NULL, 
                       
                       k.param=20,
                       prune.SNN=1/15, 
                       nn.method='annoy',
                       n.trees = 50,
                       nn.eps=0.0, 
                       annoy.metric='euclidean',
                       
                       n_neighbors = 15,
                       random_state = 0, 
                       method = 'umap', 
                       metric='euclidean',
                       generate.diffmap = FALSE,
                       n_comps = 15,
                       diffmap.name.suffix = '',
                       
                       verbose=FALSE,
                       seed=1234) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0(x, ' not in object@methods\n'))
      
    }
    
  }
  
  if(!is.numeric(k.param)) {
    
    stop('k.param must be numerical \n')
    
  }
  
  if(!is.numeric(prune.SNN)) {
    
    stop('prune.SNN must be numerical \n')
    
  }
  
  if(!is.numeric(n.trees)) {
    
    stop('n.trees must be numerical \n')
    
  }
  
  if(!is.numeric(nn.eps)) {
    
    stop('nn.eps must be numerical \n')
    
  }
  
  if(!is.character(nn.method)) {
    
    stop('nn.method must be a character string \n')
    
  } else if(is.character(nn.method)) {
    
    if(!nn.method %in% c('rann', 'annoy')) {
      
      stop('nn.method must be either: ‘rann’ or ‘annoy’ \n')
      
    }
    
  }
  
  if(!is.character(annoy.metric)) {
    
    stop('annoy.metric must be a character string \n')
    
  } else if(is.character(annoy.metric)) {
    
    if(!annoy.metric %in% c('euclidean', 'cosine', 'manhattan', 'hamming')) {
      
      stop('annoy.metric must be either: ‘euclidean’, ‘cosine’, ‘manhattan’ or ‘hamming’ \n')
      
    }
    
  }
  
  if(!is.character(neighbour.name.suffix)) {
    
    stop('neighbour.name.suffix must be character string \n')
    
  }
  
  if(!is.character(diffmap.name.suffix)) {
    
    stop('diffmap.name.suffix must be character string \n')
    
  }
  
  
  if(!is.numeric(n_neighbors)) {
    
    stop('n_neighbors must be numerical \n')
    
  }
  
  if(is.null(dims.use)) {
    
    dims.use <- list()
    
    count <- 1
    
    for(x in 1:length(reduction)) {
      
      dims.use[[count]] <- 0
      
      count <- count + 1
      
    }
    
  } else if(is.list(dims.use)) {
    
    for(x in dims.use) {
      
      if(!is.numeric(x)) {
        
        stop('dims.use items must be numerical \n')
        
      }
      
    }
    
  } else {
    
    stop('dims.use must be either NULL or a list of numerical values \n')
    
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
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical \n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  
  # if(!'integration_method' %in% colnames(object@pipelines)) {
  #   
  #   tmp <- tibble::add_column(.data = object@pipelines, integration_method=NA, integration_time=NA)
  #   
  # } else {
  #   
  #   tmp <- object@pipelines
  #   
  # }

  for(p in assay) {
    
    ass <- strsplit(x = names(object@methods)[which(names(object@methods)==p)], split = '_')[[1]][1]
    
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
          
          stop('reductions could not be found \n')
          
        }
        
      }
      
    }
    
    count <- 1
    
    if(ass %in% c('SCT','SCRAN','TPM')){
      
      for(g in reduction) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': processing assay: ', p, ', reduction: ', g, '\n')))
          
        }
        
        
        
        seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[1]]@counts)
        
        red <- reduction.list[[g]]
        red.key <- reduction[count]
        dim <-dims.use[[count]]
        
        if(!is.matrix(red)) {
          
          red <- as.matrix(red)
          
        }
        
        seuobj@reductions[[red.key]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = red, key = paste0(red.key, '_')))
        
        if(dim != 0) {
          
          if(isTRUE(verbose)) {
            
            seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                             k.param = k.param,
                                                             reduction = red.key, 
                                                             verbose = TRUE, 
                                                             dims = 1:dim, 
                                                             prune.SNN = prune.SNN,
                                                             nn.method = nn.method, 
                                                             n.trees = n.trees,
                                                             annoy.metric = annoy.metric, 
                                                             nn.eps = nn.eps))
            
          } else if (isFALSE(verbose)) {
            
            seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                             k.param = k.param,
                                                             reduction = red.key, 
                                                             verbose = FALSE, 
                                                             dims = 1:dim, 
                                                             prune.SNN = prune.SNN,
                                                             nn.method = nn.method, 
                                                             n.trees = n.trees,
                                                             annoy.metric = annoy.metric, 
                                                             nn.eps = nn.eps))
            
          }

          
        } else if(dim == 0) {
          
          if(isTRUE(verbose)) {
            
            seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                             reduction = red.key, 
                                                             verbose = TRUE, 
                                                             dims = 1:ncol(red), 
                                                             k.param = k.param,
                                                             prune.SNN = prune.SNN,
                                                             nn.method = nn.method, 
                                                             n.trees = n.trees,
                                                             annoy.metric = annoy.metric, 
                                                             nn.eps = nn.eps))
            
          } else if (isFALSE(verbose)) {
            
            seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                             reduction = red.key, 
                                                             verbose = FALSE, 
                                                             dims = 1:ncol(red), 
                                                             k.param = k.param,
                                                             prune.SNN = prune.SNN,
                                                             nn.method = nn.method, 
                                                             n.trees = n.trees,
                                                             annoy.metric = annoy.metric, 
                                                             nn.eps = nn.eps))
            
          }

        }
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': neighbours calculated\n')))
          
        }

        conn <- as(object = seuobj@graphs[[2]], Class = 'dgCMatrix')
        
        if('_' %in% unlist(strsplit(x = neighbour.name.suffix, split = ''))) {
          
          if(isTRUE(verbose)) {
            
            cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in neighbour.save.suffix, replacing with - \n')))
            
          }
          
          neighbour.name.suffix <- sub(pattern = '_', replacement = '-', x = neighbour.name.suffix)
          
        }
        
        object@methods[[p]]@neighbours[[paste0(g, '_NN', neighbour.name.suffix)]][['connectivities']] <- conn
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': results saved as ', paste0(g, '_NN', neighbour.name.suffix), '\n')))
          
        }

        count <- count + 1
        
      }
      
    } else if (ass == 'SCANPY') { 
      
      for(g in reduction) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': processing assay: ', p, ', reduction: ', g, '\n')))
          
        }
        
        sc <- reticulate::import('scanpy')
        pd <- reticulate::import('pandas')
        
        scobj <- sc$AnnData(X = t(object@methods[[2]]@norm.scaled))
        scobj$obs_names <- as.factor(colnames(object@methods[[2]]@norm.scaled))
        scobj$var_names <- as.factor(rownames(object@methods[[2]]@norm.scaled))
        
        if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
          
          scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))
          
        }
        
        red <- reduction.list[[g]]
        red.key <- reduction[count]
        dim <-dims.use[[count]]
        
        if(!is.matrix(red)) {
          
          red <- as.matrix(red)
          
        }
        
        scobj$obsm$update(X_pca = red)
        
        sc$pp$neighbors(adata = scobj, 
                        n_neighbors = as.integer(n_neighbors), 
                        n_pcs = as.integer(dim),
                        random_state = as.integer(random_state), 
                        method = as.character(method), 
                        metric = as.character(metric))
        
        if(isTRUE(generate.diffmap)) {
          
          if(isTRUE(verbose)) {
            
            cat(crayon::cyan(paste0(Sys.time(), ': calcualting diffusion map\n')))
            
          }

          sc$tl$diffmap(adata = scobj, n_comps = as.integer(n_comps))
          
          if(isTRUE(verbose)) {
            
            cat(crayon::cyan(paste0(Sys.time(), ': diffusion map calculated\n')))
            
          }

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
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': neighbours calculated\n')))
          
        }

        conn <- as.csc.matrix(scobj$obsp[['connectivities']])
        colnames(conn) <- colnames(object@methods[[2]]@norm.scaled)
        rownames(conn) <- colnames(object@methods[[2]]@norm.scaled)
        
        if('_' %in% unlist(strsplit(x = neighbour.name.suffix, split = ''))) {
          
          if(isTRUE(verbose)) {
            
            cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in neighbour.name.suffix, replacing with - \n')))
            
          }

          neighbour.name.suffix <- sub(pattern = '_', replacement = '-', x = neighbour.name.suffix)
          
        }
        
        object@methods[[p]]@neighbours[[paste0(g, '_NN', neighbour.name.suffix)]][['connectivities']] <- conn
        
        dis <- as.csc.matrix(scobj$obsp[['distances']])
        colnames(dis) <- colnames(object@methods[[2]]@norm.scaled)
        rownames(dis) <- colnames(object@methods[[2]]@norm.scaled)
        
        object@methods[[p]]@neighbours[[paste0(g, '_NN', neighbour.name.suffix)]][['distances']] <- dis
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': results saved as ', paste0(g, '_NN', neighbour.name.suffix), '\n')))
          
        }

        if(isTRUE(generate.diffmap)) {
          
          if('_' %in% unlist(strsplit(x = diffmap.name.suffix, split = ''))) {
            
            if(isTRUE(verbose)) {
              
              cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in diffmap.name.suffix, replacing with - \n')))
              
            }
            
            
            
            diffmap.name.suffix <- sub(pattern = '_', replacement = '-', x = diffmap.name.suffix)
            
          }
          
          object@methods[[p]]@computational_reductions[[paste0(g, '_NN:DIFFUSIONMAP', diffmap.name.suffix)]] <- diffmap
          
          if(isTRUE(verbose)) {
            
            cat(crayon::cyan(paste0(Sys.time(), ': diffmap saved as ', paste0(g, '_NN', diffmap.name.suffix), '\n')))
            
          }

        }
        
        count <- count + 1
        
      }
      
    } else {
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': processing assay: ', p, ', reduction: ', g, '\n')))
        
      }

      seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[1]]@counts)
      
      red <- reduction.list[[g]]
      red.key <- reduction[count]
      dim <-dims.use[[count]]
      
      if(!is.matrix(red)) {
        
        red <- as.matrix(red)
        
      }
      
      seuobj@reductions[[red.key]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = red, key = paste0(red.key, '_')))
      
      if(dim != 0) {
        
        seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                         k.param = k.param,
                                                         reduction = red.key, 
                                                         verbose = FALSE, 
                                                         dims = dim, 
                                                         prune.SNN = prune.SNN,
                                                         nn.method = nn.method, 
                                                         n.trees = n.trees,
                                                         annoy.metric = annoy.metric, 
                                                         nn.eps = nn.eps))
        
      } else if(dim == 0) {
        
        seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                         reduction = red.key, 
                                                         verbose = FALSE, 
                                                         dims = 1:ncol(red), 
                                                         k.param = k.param,
                                                         prune.SNN = prune.SNN,
                                                         nn.method = nn.method, 
                                                         n.trees = n.trees,
                                                         annoy.metric = annoy.metric, 
                                                         nn.eps = nn.eps))
        
      }
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': neighbours calculated\n')))
        
      }

      conn <- seuobj@graphs[[2]]
      
      if('_' %in% unlist(strsplit(x = neighbour.name.suffix, split = ''))) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in neighbour.save.suffix, replacing with - \n')))
          
        }

        neighbour.name.suffix <- sub(pattern = '_', replacement = '-', x = neighbour.name.suffix)
        
      }
      
      object@methods[[p]]@neighbours[[paste0(g, '_NN', neighbour.name.suffix)]][['connectivities']] <- conn
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': results saved as ', paste0(g, '_NN', neighbour.name.suffix), '\n')))
        
      }

      count <- count + 1
      
    }
    
  }
  
  return(object)
  
}