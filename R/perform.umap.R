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
#' @param n.dims Numerical. The number of Scanorama dimensions to be produced. Default = 50
#' @param n_components Numerical. How many UMAP dimensions should be produced, if you are supplying graphs, only 2 dimensions can be produced. Default = 3
#' @param ... Numerical. Arguments to be passed to Seurat::RunUMAP
#' 
#' @return UMAP reduction saved in the visualisation_reductions section in the supplied method-assays
#'
#' @export

perform.umap <- function(object, 
                         assay,
                         reduction=NULL,
                         graph=NULL,
                         reduction.save='umap',
                         n.dims=NULL, 
                         n_components = 3, 
                         ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string \n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan('assay: ', x, 'does not exist\n'))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string \n'))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    cat(crayon::cyan('n_components must be numerical\n'))
    return(object)
    
  }
  
  if(!is.null(reduction) & !is.null(graph)) {
    
    cat(crayon::cyan('only graphs OR reductions can be provided\n'))
    return(object)
    
  }
  
  for(u in assay) {
    
    if(is.null(reduction) & !is.null(graph)) {
      
      count <- 1
      
      graphs <- object@methods[[u]]@graphs
      
      for (g in graph) {
        
        cat(crayon::cyan('Processing', g, 'for assay:', u,'\n'))
        
        red.save <- reduction.save[count]
        
        seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[u]]@counts))
        
        seuobj[['temp']] <- suppressWarnings(Seurat::as.Graph(object@methods[[u]]@graphs[[g]]$connectivities))
        
        seuobj <- suppressWarnings(Seurat::RunUMAP(object = seuobj, 
                                                   graph = 'temp', 
                                                   n_components = n_components, 
                                                   verbose = TRUE,
                                                   ...))
        
        red.iso <- Seurat::Embeddings(object = seuobj, 
                                      reduction = 'umap')
        
        rownames(red.iso) <- colnames(object)
        
        dim.names <- list()
        
        for(l in 1:2) {
          
          dim.names[[l]] <- paste0('umap_', l)
          
        }
        
        colnames(red.iso) <- unlist(dim.names)
        
        object@methods[[u]]@visualisation_reductions[[red.save]] <- red.iso
        
        count <- count + 1
        
      }
      
    } else {
      
      if(!is.list(n.dims)) {
        
        cat(crayon::cyan('dimensions must be supplied in list format\n'))
        return(object)
        
      }
      
      if(!is.character(reduction)) {
        
        cat(crayon::cyan('reduction must be character string \n'))
        return(object)
        
      }
      
      for(r in reduction) {
        
        if(!r %in% c(names(object@methods[[u]]@computational_reductions), 
                     names(object@methods[[u]]@integration_reductions),
                     names(object@methods[[u]]@visualisation_reductions))) {
          
          cat(crayon::cyan(paste0('reduction:', r, ' does not exist\n')))
          return(object)
          
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
        
        dim <- n.dims[[count]]
        
        red.save <- reduction.save[count]
        
        red <- reduction.list[[i]]
        
        cat(crayon::cyan('Processing', i, 'for assay:', u,'\n'))
        
        if(!is.null(dim)) {
          
          c <- uwot::umap(X = red[,dim], n_components = n_components, verbose = TRUE, ...)
          
        } else {
          
          c <- uwot::umap(X = red, n_components = n_components, verbose = TRUE, ...)
          
        }
        
        dim.names <- list()
        
        for(l in 1:n_components) {
          
          dim.names[[l]] <- paste0('umap_', l)
          
        }
        
        colnames(c) <- unlist(dim.names)
        
        rownames(c) <- colnames(object)
        
        object@methods[[u]]@visualisation_reductions[[red.save]] <- c
        
        count <- count + 1
        
      }
      
    }
    
  }
  
  return(object)
  
}
