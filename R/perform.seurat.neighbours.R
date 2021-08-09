#' @name perform.seurat.neighbours
#'
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm.
#' @param neighbour.name Character. String defining the names to store the neighbourhood graph results under.
#' @param k.param Numerical. The number of k when calcuating k-nearest neighbour. Default = 20
#' @param compute.SNN Boolean. Should the shared nearest neighbour graph be calculated. Default = TRUE 
#' @param prune.SNN Numerical. Setas acceptance cutoff for jaccard index whilst computing neighbourhood overlap for SNN construction. Any edges with a value less than this parameter will be removed. 0 = no pruning and 1 = prune everything. Default = 0
#' @param nn.method Character. Nearest neighbour method, either 'rann' or 'annoy'. Default = 'annoy'
#' @param n.trees Numerical. More trees facilitates hgiher precision when using 'annoy' method. Default = 20
#' @param nn.eps Numerical. Margin of error when performing nearest neighbour search whilst using rann method. 0 would imply an exact search. Default = 0.0
#' @param annoy.metric Character. Distance metric for annoy method. Options: 'euclidean', 'cosine', 'manhattan', 'hamming'. Default = 'euclidean'
#' 
#' @export

perform.nn.v2 <- function(object, 
                          assay, 
                          reduction, 
                          neighbour.name.suffix = '',
                          dims=NULL, 
                          k.param=20,
                          prune.SNN=1/15, 
                          nn.method='annoy', 
                          n.trees = 50,
                          nn.eps=0.0, 
                          annoy.metric='euclidean') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0(x, ' is not contained within methods \n'))
      
    }
    
  }
  
  if(!is.character(neighbour.name.suffix)) {
    
    stop('neighbour.name.suffix must be a character string \n')
    
  }
  
  if(!is.list(dims)) {
    
    stop('dims must be a list \n')
    
  } else if (is.list(dims)) {
    
    for (x in dims) {
      
      if(!is.numeric(x)) {
        
        stop('dims items must be numerical \n')
        
      }
      
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
      
      seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[1]]@counts)
      
      red <- reduction.list[[g]]
      red.key <- reduction[count]
      dim <-dims[[count]]
      
      seuobj@reductions[[red.key]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = as.matrix(red), key = paste0(red.key, '_')))
      
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
      
      cat(crayon::cyan(paste0(Sys.time(), ': neighbours calculated\n')))
      
      conn <- seuobj@graphs[[2]]
      
      if('_' %in% unlist(strsplit(x = neighbour.name.suffix, split = ''))) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in neighbour.save.suffix, replacing with - \n')))
        neighbour.name.suffix <- sub(pattern = '_', replacement = '-', x = neighbour.name.suffix)
        
      }
      
      object@methods[[p]]@neighbours[[paste0(g, '_nn.v2', neighbour.name.suffix)]][['connectivities']] <- conn
      
      cat(crayon::cyan(paste0(Sys.time(), ': results saved as ', paste0(g, '_nn.v2', neighbour.name.suffix), '\n')))
      
      count <- count + 1
      
    }
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': completed \n')))
  
  return(object)
  
}
