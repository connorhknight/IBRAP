#' @name perform.seurat.cluster
#' @aliases perform.seurat.cluster
#' 
#' @title Performs Seurat clustering
#'
#' @description Performs Seurat clustering on defined method-assays and supplied reductions or graphs. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm. Default = NULL
#' @param dims Numerical. How many dimensions of the reduciton should be supplied, NULL equates to all/. Default = NULL
#' @param graph Character. Name of the graph to use for clustering (i.e. BBKNN integrated graph or previously calculated nearest neighbour graphs), only either graphs or reductions can be supplied in one instance. Default = NULL
#' @param k.param Numerical. The number of k when calcuating k-nearest neighbour. Default = 20
#' @param compute.SNN Boolean. Should the shared nearest neighbour graph be calculated. Default = TRUE 
#' @param prune.SNN Numerical. Setas acceptance cutoff for jaccard index whilst computing neighbourhood overlap for SNN construction. Any edges with a value less than this parameter will be removed. 0 = no pruning and 1 = prune everything. Default = 0
#' @param nn.method Character. Nearest neighbour method, either 'rann' or 'annoy'. Default = 'annoy'
#' @param n.trees Numerical. More trees facilitates hgiher precision when using 'annoy' method. Default = 20
#' @param nn.eps Numerical. Margin of error when performing nearest neighbour search whilst using rann method. 0 would imply an exact search. Default = 0.0
#' @param annoy.metric Character. Distance metric for annoy method. Options: 'euclidean', 'cosine', 'manhattan', 'hamming'. Default = 'euclidean'
#' @param algorithm Numerical. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default = 1 Default = NULL
#' @param assignment.df.name Character. What to call the df contained in clusters. Default = 'seurat
#' @param res Numerical vector. Which resolution to run the clusterign algorithm at, a smaller and larger value identified less and more clusters, respectively. Default = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
#' @param ... arguments to be passed to Seurat::FindClusters
#' 
#' @return Cluster assignments using the list of resolutions provided contained within cluster_assignments under assignment.df.name
#'
#' @export

perform.seurat.cluster <- function(object, 
                                   assay,
                                   reduction=NULL,
                                   dims=NULL,
                                   graph=NULL,
                                   k.param=20,
                                   compute.SNN=TRUE,
                                   prune.SNN=1/15, 
                                   nn.method='annoy', 
                                   n.trees = 50,
                                   nn.eps=0.0, 
                                   annoy.metric='euclidean',
                                   algorithm=1,
                                   assignment.df.name='seurat',
                                   res=c(0.1,0.2,0.3,0.4,0.5,
                                         0.6,0.7,0.8,0.9,1,
                                         1.1,1.2,1.3,1.4,1.5), 
                                   ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  
  
  if(!is.null(reduction)) {
    
    if(is.list(dims)) {
      
      if(length(dims) != length(reduction)) {
        
        cat(crayon::cyan('number of dims must match number of reduction(s)\n'))
        return(object)
        
      }
      
      if(length(dims) != length(assignment.df.name)) {
        
        cat(crayon::cyan('number of dims must match number of assignment.df.name(s)\n'))
        return(object)
        
      }
      
      if(!is.character(assay)) {
        
        cat(crayon::cyan('assay must be character string\n'))
        return(object)
        
      }
      
    }
    
    if(!is.list(dims)) {
      
      cat(crayon::cyan('dims must be provided in a list\n'))
      return(object)
      
    }
    
    if(!is.character(reduction)) {
      
      cat(crayon::cyan('reduction must be character string\n'))
      return(object)
      
    }
    
    for(x in reduction) {
      
      for(i in assay) {
        
        if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                           object@methods[[i]]@visualisation_reductions, 
                           object@methods[[i]]@integration_reductions))) {
          
          cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
          return(object)
          
        }
        
      }
      
    }
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('assay: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  
  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan('assignment.df.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(res)) {
    
    cat(crayon::cyan('res must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(algorithm)) {
    
    cat(crayon::cyan('method must either be numeric\n'))
    return(object)
    
  }
  
  if(!algorithm %in% c(1,2,3,4)) {
    
    cat(crayon::cyan('method must either: 1 (original Louvain), 2 (Louvain with multilevel refinement), 3 (SLM) or 4 (Leiden)\n'))
    return(object)
    
  }
  
  if(!is.numeric(prune.SNN)) {
    
    cat(crayon::cyan('prune.SNN must either be numeric\n'))
    return(object)
    
  }
  
  if(!is.character(nn.method)) {
    
    cat(crayon::cyan('nn.method must either be character string\n'))
    return(object)
    
  }
  
  if(!is.character(annoy.metric)) {
    
    cat(crayon::cyan('annoy.metric must either be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(nn.eps)) {
    
    cat(crayon::cyan('nn.eps must either be numeric\n'))
    return(object)
    
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
          
          cat(crayon::cyan('reductions could not be found\n'))
          return(object)
          
        }
        
      }
      
    }
    
    if(!is.null(graph)) {
      
      for(g in graph) {
        
        if(!g %in% names(object@methods[[p]]@graphs)) {
          
          cat(crayon::cyan('graphs could not be found\n'))
          return(object)
          
        }
        
      }
      
    }
    
    count <- 1
    
    if(is.null(reduction) & !is.null(graph)) {
      
      for(t in graph) {
        
        graph.iso <- object@methods[[p]]@graphs[[t]]$connectivities
        
        cat(crayon::cyan(paste0('Calculating Seurat clusters for assay: ', p, ' using graph: ', t, '\n')))
        
        tmp <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
        
        orig.name <- names(tmp@meta.data)
        
        tmp[['temp']] <- Seurat::as.Graph(x = graph.iso)
        
        cat(crayon::cyan('Graph added\n'))
        
        tmp <- suppressWarnings(Seurat::FindClusters(object = tmp, resolution = res, graph.name = 'temp', algorithm = algorithm, verbose = FALSE, ...))
        
        cat(crayon::cyan('Clusters identified\n'))
        
        post.name <- names(tmp@meta.data)
        
        post.name <- post.name[!post.name %in% orig.name]
        
        post.name <- post.name[1:length(post.name)-1]
        
        new.clusters <- tmp@meta.data[,post.name]
        
        new.clusters <- new.clusters[match(rownames(object@sample_metadata), rownames(new.clusters)),]
        
        object@methods[[p]]@cluster_assignments[[assignment.df.name[count]]] <- new.clusters
        
        cat(crayon::cyan('Seurat clusters added\n'))
        
        count <- count + 1
        
      }
      
      
      
    } else {
      
      for(g in reduction) {
        
        dim <- dims[[count]]
        red <- reduction.list[[g]]
        red.key <- reduction[count]
        
        cat(crayon::cyan(paste0('Calculating Seurat clusters for assay: ', p, ' using reduction: ', g, '\n')))
        
        tmp <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
        
        tmp@reductions[[red.key]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = as.matrix(red), key = red.key))
        
        orig.name <- names(tmp@meta.data)
        
        if(!is.null(dim)) {
          
          tmp <- suppressWarnings(Seurat::FindNeighbors(object = tmp, 
                                                        k.param = k.param,
                                                        reduction = red.key, 
                                                        verbose = FALSE, 
                                                        dims = dim, 
                                                        compute.SNN = compute.SNN, 
                                                        prune.SNN = prune.SNN,
                                                        nn.method = nn.method, 
                                                        n.trees = n.trees,
                                                        annoy.metric = annoy.metric, 
                                                        nn.eps = nn.eps))
          
        } else {
          
          tmp <- suppressWarnings(Seurat::FindNeighbors(object = tmp, 
                                                        reduction = red.key, 
                                                        verbose = FALSE, 
                                                        dims = 1:ncol(red), 
                                                        k.param = k.param,
                                                        compute.SNN = compute.SNN, 
                                                        prune.SNN = prune.SNN,
                                                        nn.method = nn.method, 
                                                        n.trees = n.trees,
                                                        annoy.metric = annoy.metric, 
                                                        nn.eps = nn.eps))
          
        }
        
        cat(crayon::cyan('Neighbours identified\n'))
        
        tmp <- suppressWarnings(Seurat::FindClusters(object = tmp, resolution = res, algorithm = algorithm, verbose = FALSE, ...))
        
        cat(crayon::cyan('Clusters identified\n'))
        
        post.name <- names(tmp@meta.data)
        
        post.name <- post.name[!post.name %in% orig.name]
        
        post.name <- post.name[1:length(post.name)-1]
        
        new.clusters <- tmp@meta.data[,post.name]
        
        new.clusters <- new.clusters[match(rownames(object@sample_metadata), rownames(new.clusters)),]
        
        object@methods[[p]]@cluster_assignments[[assignment.df.name[count]]] <- new.clusters
        
        cat(crayon::cyan('Seurat clusters added\n'))
        
        count <- count + 1
        
      }
      
    }
    
  }
  
  return(object)
  
}
