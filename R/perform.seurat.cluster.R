#' @name perform.seurat.cluster
#' @aliases perform.seurat.cluster
#' 
#' @title Performs Seurat clustering
#'
#' @description Performs Seurat clustering on defined method-assays and supplied reductions or graphs. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param neighbours Character. String indicating which neighbourhood graphs should be used. 
#' @param algorithm Numerical. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default = 1 Default = NULL
#' @param cluster.df.name Character. What to call the df contained in clusters. Default = 'seurat
#' @param res Numerical vector. Which resolution to run the clusterign algorithm at, a smaller and larger value identified less and more clusters, respectively. Default = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
#' @param ... arguments to be passed to Seurat::FindClusters
#' 
#' @return Cluster assignments using the list of resolutions provided contained within cluster_assignments under cluster.df.name
#'
#' @export

perform.seurat.cluster <- function(object, 
                                   assay,
                                   neighbours, 
                                   algorithm=1,
                                   cluster.df.name,
                                   res=c(0.1,0.2,0.3,0.4,0.5,
                                         0.6,0.7,0.8,0.9,1,
                                         1.1,1.2,1.3,1.4,1.5), 
                                   ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string(s)')
    
  } else if (is.character(assay)) {
    
    for(x in assay) {
      
      if(!x %in% names(object@methods)) {
        
        stop(paste0('assay: ', x, ' is not contained within object@methods \n'))
        
      }
      
    }
    
  }
  
  for(p in assay) {
    
    if(!is.character(neighbours)) {
      
      stop('neighbours must be character string(s)')
      
    } else if (is.character(neighbours)) {
      
      for(x in neighbours) {
        
        if(!x %in% names(object@methods[[p]]@neighbours)) {
          
          stop(paste0('neighbours: ', x, ' is not contained within assay: ', p,  ' \n'))
          
        }
        
      }
      
    }
    
  }
  
  if(!is.numeric(algorithm)) {
    
    stop('method must either be numeric\n')
    
  }
  
  if(!algorithm %in% c(1,2,3,4)) {
    
    stop('method must either: 1 (original Louvain), 2 (Louvain with multilevel refinement), 3 (SLM) or 4 (Leiden)\n')
    
  }
  
  for(p in assay) {
    
    count <- 1

      for(t in neighbours) {

        cat(crayon::cyan(paste0(Sys.time(), ': calculating Seurat clusters for assay: ', p, ' using graph: ', t, '\n')))
        
        tmp <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
        
        orig.name <- names(tmp@meta.data)
        
        tmp@graphs[['neighbourhood_graph']] <- suppressWarnings(Seurat::as.Graph(x = object@methods[[p]]@neighbours[[t]]$connectivities))
        
        cat(crayon::cyan(paste0(Sys.time(), ': graph added\n')))
        
        tmp <- suppressWarnings(Seurat::FindClusters(object = tmp, resolution = res, graph.name = 'neighbourhood_graph', 
                                                     algorithm = algorithm, verbose = FALSE, ...))
        
        cat(crayon::cyan(paste0(Sys.time(), ': clusters identified\n')))
        
        post.name <- names(tmp@meta.data)
        
        post.name <- post.name[!post.name %in% orig.name]
        
        post.name <- post.name[1:length(post.name)-1]
        
        new.clusters <- tmp@meta.data[,post.name]
        
        new.clusters <- new.clusters[match(rownames(object@sample_metadata), rownames(new.clusters)),]
        
        object@methods[[p]]@cluster_assignments[[cluster.df.name[count]]] <- new.clusters
        
        cat(crayon::cyan(paste0(Sys.time(), ': seurat clusters added\n')))
        
        count <- count + 1
        
      }
    
  }
  
  return(object)
  
}
