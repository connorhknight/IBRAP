#' @name perform.umap
#' @aliases perform.umap
#' 
#' @title Performs UMAP reduction
#'
#' @description Performs UMAP reduction on defined method-assays and supplied reductions or graphs.
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the UMAP algorithm. Default = NULL
#' @param graph Character. If you wish to UMAP project a previously created connectivity graph (i.e. BBKNN output) supply the graph name here with reductions set to NULL. Default = NULL
#' @param reduction.name.suffix Character. What should be appended to the end of umap as the reduction name. 
#' @param n.dims Numerical. The number of UMAP dimensions to be produced. Must be supplied in list format relative to the order of reductions
#' @param n_components Numerical. How many UMAP dimensions should be produced, if you are supplying graphs, only 2 dimensions can be produced. Default = 3
#' @param n_neighbors Numerical. How many neighbours should be identified per cell. A higher value typically returns more accurate results. Default = 
#' @param verbose Logical. Should system information be printed. Default = FALSE
#' @param seed Numeric. What should the seed be set as. Default = 1234
#' @param ... Numerical. Arguments to be passed to Seurat::RunUMAP
#' 
#' @return UMAP reduction saved in the visualisation_reductions section in the supplied method-assays
#' 
#' @examples 
#' 
#' # How to perform umap reduction on a previous reduction 
#' object <- perform.umap(object = object, 
#'                        assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                        reduction = c('pca'), 
#'                        n_components = 2, 
#'                        n.dims = list(1:12))
#'                        
#' # How to perfrom umap reduction on neighbourhood graphs
#' object <- perform.umap(object = object, 
#'                        assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                        graph = 'pca_bbknn_bbknn')
#'                       
#'                       
#'
#' @export

perform.umap <- function(object, 
                         assay,
                         reduction=NULL,
                         graph=NULL,
                         reduction.name.suffix='',
                         dims.use=NULL, 
                         n_components = 2, 
                         n_neighbors = 30,
                         metric = 'cosine',
                         min_dist = 0.3,
                         verbose=FALSE,
                         seed=1234,
                         ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string \n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop('assay: ', x, 'does not exist\n')
      
      
    }
    
  }
  
  if(!is.character(reduction.name.suffix)) {
    
    stop('reduction.name.suffix must be character string \n')
    
  }
  
  if(!is.numeric(n_components)) {
    
    stop('n_components must be numerical\n')
    
  }
  
  if(!is.null(reduction) & !is.null(graph)) {
    
    stop('only graphs OR reductions can be provided\n')
    
  }
  
  if('_' %in% unlist(strsplit(x = reduction.name.suffix, split = ''))) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in reduction.name.suffix, replacing with - \n')))
      
    }
    
    reduction.name.suffix <- sub(pattern = '_', replacement = '-', x = reduction.name.suffix)
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical \n')
    
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
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  
  for(u in assay) {
    
    if(is.null(reduction) & !is.null(graph)) {
      
      count <- 1
      
      graphs <- object@methods[[u]]@neighbours
      
      for (g in graph) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': processing ', g, ' for assay: ', u,'\n')))
          
        }

        seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[u]]@counts))
        
        seuobj@graphs[['temp']] <- suppressWarnings(Seurat::as.Graph(object@methods[[u]]@neighbours[[g]]$connectivities))
        
        if(isTRUE(verbose)) {
          
          seuobj <- suppressWarnings(Seurat::RunUMAP(object = seuobj, 
                                                     graph = 'temp', 
                                                     n_components = n_components, 
                                                     verbose = TRUE,
                                                     ...))
          
        } else if(isFALSE(verbose)) {
          
          seuobj <- suppressWarnings(Seurat::RunUMAP(object = seuobj, 
                                                     graph = 'temp', 
                                                     n_components = n_components, 
                                                     verbose = FALSE,
                                                     ...))
          
        }
        

        
        red.iso <- Seurat::Embeddings(object = seuobj, 
                                      reduction = 'umap')
        
        rownames(red.iso) <- colnames(object)
        
        dim.names <- list()
        
        for(l in 1:2) {
          
          dim.names[[l]] <- paste0('umap_', l)
          
        }
        
        colnames(red.iso) <- unlist(dim.names)
        
        object@methods[[u]]@visualisation_reductions[[paste0(g, ':UMAP', reduction.name.suffix)]] <- red.iso
        
        count <- count + 1
        
      }
      
    } else {
      
      if(!is.character(reduction)) {
        
        stop('reduction must be character string \n')
        
      }
      
      for(r in reduction) {
        
        if(!r %in% c(names(object@methods[[u]]@computational_reductions), 
                     names(object@methods[[u]]@integration_reductions),
                     names(object@methods[[u]]@visualisation_reductions))) {
          
          stop(paste0('reduction:', r, ' does not exist\n'))
          
        }
        
      }
      
      reduction.list <- list()
      
      red.names <- c(names(object@methods[[u]]@computational_reductions), 
                     names(object@methods[[u]]@integration_reductions),
                     names(object@methods[[u]]@visualisation_reductions))
      
      for(i in red.names) {
        
        if(i %in% names(object@methods[[u]]@computational_reductions)) {
          
          reduction.list[[i]] <- object@methods[[u]]@computational_reductions[[i]]
          
        }
        
        if(i %in% names(object@methods[[u]]@integration_reductions)) {
          
          reduction.list[[i]] <- object@methods[[u]]@integration_reductions[[i]]
          
        }
        
        if(i %in% names(object@methods[[u]]@visualisation_reductions)) {
          
          reduction.list[[i]] <- object@methods[[u]]@visualisation_reductions[[i]]
          
        }
        
      }
      
      count <- 1
      
      for(i in reduction) {
        
        dim <- dims.use[[count]]
        
        red <- reduction.list[[i]]
        
        if(dim == 0) {
          
          dim <- ncol(red)
          
        }
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': processing ', i, ' for assay: ', u,'\n')))
          
        }

        seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[u]]@counts))
        
        seuobj@reductions$pca <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = red, assay = 'RNA', key = paste0(i, '_')))
        
        if(isTRUE(verbose)) {
          
          seuobj <- suppressWarnings(Seurat::RunUMAP(object = seuobj, 
                                                     dims = 1:dim,
                                                     n_components = n_components, 
                                                     reduction = 'pca',
                                                     verbose = TRUE,
                                                     ...))
          
        } else if(isFALSE(verbose)) {
          
          seuobj <- suppressWarnings(Seurat::RunUMAP(object = seuobj, 
                                                     dims = 1:dim,
                                                     n_components = n_components, 
                                                     reduction = 'pca',
                                                     verbose = FALSE,
                                                     ...))
          
        }
        
        red.iso <- Seurat::Embeddings(object = seuobj, 
                                      reduction = 'umap')
        
        rownames(red.iso) <- colnames(object)
        
        dim.names <- list()
        
        for(l in 1:2) {
          
          dim.names[[l]] <- paste0('umap_', l)
          
        }
        
        colnames(red.iso) <- unlist(dim.names)
        
        object@methods[[u]]@visualisation_reductions[[paste0(i, '_UMAP', reduction.name.suffix)]] <- red.iso
        
        count <- count + 1
        
      }
      
    }
    
  }
  
  return(object)
  
}
