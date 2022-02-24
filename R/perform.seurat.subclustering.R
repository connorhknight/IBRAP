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
#' @param verbose Logical Should function messages be printed?
#' @param seed Numeric. What should the seed be set as. Default = 1234
#' @param ... arguments to be passed to Seurat::FindClusters
#' 
#' @return A new column within the defined cluster_assignment dataframe containing original and new subclusters
#' 
#' @examples 
#' 
#' object <- perform.graph.subclustering(object = object, assay = 'SCT', 
#'                                       clust.method = 'pca', 
#'                                       column = 'neighbourhood_graph_res.0.7', clusters = c(1,5,9), 
#'                                       neighbours = 'pca_nn:', algorithm = 1)
#'
#' @export

perform.graph.subclustering <- function(object, 
                                        assay, 
                                        clust.method, 
                                        column, 
                                        clusters, 
                                        neighbours, 
                                        algorithm = 1, 
                                        res = 0.6,
                                        verbose=FALSE,
                                        seed=1234, 
                                        ...) {
  
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
    
  } else if (clust.method != 'metadata') {
    
    if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
      
      stop('clust.method is not contained within cluster_assignments \n')
      
    }
    
  }
  
  if(!is.character(column)) {
    
    stop('column must be character string \n')
    
  } else if(clust.method != 'metadata') {
    
    if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
      
      stop('column is not contained within the defined cluster.method dataframe \n')
      
    }
    
  } else if (clust.method == 'metadata') {
    
    if(!column %in% colnames(object@sample_metadata)) {
      
      stop('column is not contained within the defined metadata \n')
      
    }
    
  }
  
  for(x in clusters) {
    
    if(clust.method != 'metadata') {
      
      if(!x %in% object@methods[[assay]]@cluster_assignments[[clust.method]][,column]) {
        
        stop('defined clusters are not contained within the dataframe column \n')
        
      }
      
    } else if (clust.method == 'metadata') {
      
      if(!x %in% object@sample_metadata[,column]) {
        
        stop('defined clusters are not contained within the metadata column \n')
        
      }
      
    }
    
  }
  
  if(algorithm==1) {
    
    algo.name <- 'louvain'
    
  }
  
  if(algorithm==2) {
    
    algo.name <- 'louvainMLR'
    
  }
  
  if(algorithm==3) {
    
    algo.name <- 'SLM'
    
  }
  
  if(algorithm==4) {
    
    algo.name <- 'leiden'
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical\n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  
  if(clust.method != 'metadata') {
    
    cell_subset <- object[,object@methods[[assay]]@cluster_assignments[[clust.method]][,column] %in% clusters]
    
    cell_subset <- perform.graph.cluster(object = cell_subset, assay = assay, neighbours = neighbours, res = res, algorithm = 1, ...)
    
    subclusters <- as.character(cell_subset@methods[[assay]]@cluster_assignments[[paste0(neighbours, ':', algo.name)]][,1])
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': identified ', length(unique(subclusters)), ' subclusters\n')))
      
    }

    subclusters <- as.data.frame(as.factor(subclusters))
    
    rownames(subclusters) <- colnames(cell_subset)
    
    object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subcluster_', res)] <- as.character(object@methods[[assay]]@cluster_assignments[[clust.method]][,column])
    
    sub_ids <- which(object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subcluster_', res)] %in% as.character(clusters))
    
    orig.names <- as.data.frame(object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subcluster_', res)][sub_ids])
    
    rownames(orig.names) <- colnames(cell_subset)
    
    orig.names[,2] <- subclusters[match(x = rownames(subclusters), table = rownames(orig.names)),]
    
    orig.names[,3] <- paste0(orig.names[,1], '_', orig.names[,2])
    
    object@methods[[assay]]@cluster_assignments[[clust.method]][,paste0(column, '_subcluster_', res)][sub_ids] <- orig.names[,3]
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': subclusters added under column name: ', paste0(column, '_subcluster_', res), '\n')))
      
    }

  } else if (clust.method == 'metadata') {
    
    object2 <- object
    
    object2@methods[[assay]]@cluster_assignments[['metadata']] <- object@sample_metadata
    
    cell_subset <- object2[,object2@methods[[assay]]@cluster_assignments[['metadata']][,column] %in% clusters]
    
    print(dim(cell_subset))
    
    cell_subset <- perform.graph.cluster(object = cell_subset, assay = assay, neighbours = neighbours, res = res, algorithm = 1, ...)

    subclusters <- as.character(cell_subset@methods[[assay]]@cluster_assignments[[paste0(neighbours, ':', algo.name)]][,1])

    subclusters <- as.data.frame(as.factor(subclusters))
 
    rownames(subclusters) <- colnames(cell_subset)

    object2@methods[[assay]]@cluster_assignments[['metadata']][,paste0(column, '_subcluster_', res)] <- as.character(object2@methods[[assay]]@cluster_assignments[['metadata']][,column])

    sub_ids <- which(object2@methods[[assay]]@cluster_assignments[['metadata']][,paste0(column, '_subcluster_', res)] %in% as.character(clusters))

    orig.names <- as.data.frame(object2@methods[[assay]]@cluster_assignments[['metadata']][,paste0(column, '_subcluster_', res)][sub_ids])

    rownames(orig.names) <- colnames(cell_subset)

    orig.names[,2] <- subclusters[match(x = rownames(subclusters), table = rownames(orig.names)),]

    orig.names[,3] <- paste0(orig.names[,1], '_', orig.names[,2])

    object2@methods[[assay]]@cluster_assignments[['metadata']][,paste0(column, '_subcluster_', res)][sub_ids] <- as.character(orig.names[,3])
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': identified ', length(unique(subclusters)), ' subclusters\n')))
      
    }
    
    object@sample_metadata[,paste0(column, '_subcluster_', res)] <- object2@methods[[assay]]@cluster_assignments[['metadata']][,paste0(column, '_subcluster_', res)]
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': subclusters added under column name: ', paste0(column, '_subcluster_', res), '\n')))
      
    }

  }
  
  return(object)
  
}
