#' @name perform.graph.subclustering
#' @aliases perform.graph.subclustering
#' 
#' @title Seurat subclustering
#'
#' @param object An IBRAP S4 class object
#' @param assay Character. Which assay within the object to access
#' @param clust.method Character. Which cluster_assignments dataframe to access
#' @param column Character. Which column to access within the cluster_assignment dataframe
#' @param clusters Which cluster(s) would you like to subcluster
#' @param neighbours Character. String indicating which neighbourhood graphs should be used. 
#' @param algorithm Numerical. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default = 1 Default = NULL
#' @param cluster.df.name Character. What to call the df contained in clusters. Default = 'seurat
#' @param res Numerical vector. Which resolution to run the clusterign algorithm at, a smaller and larger value identified less and more clusters, respectively. Default = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
#' @param ... arguments to be passed to Seurat::FindClusters
#' 
#' @return A new column within the defined cluster_assignment dataframe containing original and new subclusters
#'
#' @export

perform.graph.subclustering <- function(object, 
                                        assay, 
                                        clust.method, 
                                        column, 
                                        clusters, 
                                        neighbours, 
                                        algorithm, 
                                        res = 0.6, ...) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string \n')
    
  } else {
    
    if(!assay %in% names(object@methods)) {
      
      stop('assay is not contained within object@methods')
      
    }
    
  }
  
  if(!is.character(clust.method)) {
    
    stop('clust.method must be a character string \n')
    
  } else {
    
    if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
      
      stop('clust.method is not contained within cluster_assignments \n')
      
    }
    
  }
  
  if(!is.character(column)) {
    
    stop('column must be character string \n')
    
  } else {
    
    if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
      
      stop('column is not contained within the defined cluster.method dataframe \n')
      
    }
    
  }
  
  for(x in clusters) {
    
    if(!x %in% object@methods[[assay]]@cluster_assignments[[clust.method]][,column]) {
      
      stop('defined clusters are not contained within the dataframe \n')
      
    }
    
  }
  
  cell_subset <- object[,object@methods[[assay]]@cluster_assignments[[clust.method]][,column]==clusters]
  
  cell_subset <- perform.seurat.cluster(object = cell_subset, assay = assay, neighbours = neighbours, res = res, cluster.df.name = clust.method, ...)
  
  
  subclusters <- as.numeric(cell_subset@methods[[assay]]@cluster_assignments[[clust.method]][,1])
  
  cat(crayon::cyan(paste0(Sys.time(), ': identified ', length(unique(subclusters)), ' subclusters\n')))
  
  if(length(unique(subclusters)) > 26) {
    
    let <- paste0(LETTERS, LETTERS)
    
  } else {
    
    let <- LETTERS
    
  }
  
  minimum <- min(subclusters)
  maximum <- max(subclusters)
  
  for(x in minimum:maximum) {
    
    subclusters[subclusters == x] <- let[x]
    
  }
  
  subclusters <- as.data.frame(as.factor(subclusters))
  rownames(subclusters) <- colnames(cell_subset)
  
  object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subres_', res)] <- as.character(object@methods[[assay]]@cluster_assignments[[clust.method]][,column])
  
  sub_ids <- object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subres_', res)]==as.character(clusters)
  
  object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subres_', res)][sub_ids] <- as.character(subclusters[,1])
  
  cat(crayon::cyan(paste0(Sys.time(), ': subclusters added under column name: ', paste0(column, '_subres_', res), '\n')))
  
  return(object)
  
}