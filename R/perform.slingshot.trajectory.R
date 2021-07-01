#' @name perform.slingshot.trajectory
#' @aliases perform.slingshot.trajectory
#' 
#' @title Performs Slingshot trajectory inference
#'
#' @description Generates slingshot trajectory inference on the defined reduction and clustering
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm. Default = NULL
#' @param clust.method Character. Which cluster method should be used utilised from clustering results, if `'metadata'` is supplied, you will access the metadata.
#' @param column Character. Which column within the isolated clust.method should be used to define cell type annoation. 
#' @param start.clus Character. Which cluster should start the trajectory, if NULL then slingshot will attempt to discover this. Default = NULL
#' @param end.clus Character. Which cluster should end the trajectory, if NULL then slingshot will attempt to discover this. Default = NULL
#' @param ... arguments to be passed to slingshot::slingshot
#' 
#' @return A SingleshotDataSet class results object containing cellular lineages/curves
#'
#' @export

perform.slingshot.trajectory <- function(object, 
                                         reduction, 
                                         assay,
                                         clust.method,
                                         column, 
                                         start.clus = NULL, 
                                         end.clus = NULL,
                                         ...) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('Object must be of class IBRAP \n')
    
  }
  
  if(!is.character(reduction)) {
    
    stop('Reduction must be character string \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('Assay must be character string \n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('Assay does not exist \n')
    
  }
  
  if(!reduction %in% c(names(object@methods[[assay]]@computational_reductions), 
                       names(object@methods[[assay]]@integration_reductions),
                       names(object@methods[[assay]]@visualisation_reductions))) {
    
    stop('Reduction not present in assay \n')
    
  }
  
  reduction.list <- list()
  red.names <- c(names(object@methods[[assay]]@computational_reductions), 
                 names(object@methods[[assay]]@integration_reductions),
                 names(object@methods[[assay]]@visualisation_reductions))
  
  for(i in red.names) {
    
    if(i %in% names(object@methods[[assay]]@computational_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@computational_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@integration_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@integration_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@visualisation_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@visualisation_reductions[[i]]
      
    }
  }
  
  red <- reduction.list[[reduction]]
  
  if(!is.character(clust.method)) {
    
    stop('clust.method should be a character string \n')
    
  }
  
  if(clust.method != 'metadata') {
    
    if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
      
      stop('clust.method should either be metadata or cluster assignment data.frame name \n')
      
    }
    
    if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
      
      stop(paste0(column, ' does not exist in the defined clust.method dataframe \n'))
      
    } else if (column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
      
      clusters <- object@methods[[assay]]@cluster_assignments[[clust.method]][,column]
      
      if(is.null(clusters)) {
        
        stop(paste0('error, defined column is null \n'))
        
      }
      
    }
    
  } else if (clust.method == 'metadata') {
    
    if(!column %in% colnames(object@sample_metadata)) {
      
      stop(paste0('error, defined column is null \n'))
      
    } else if (column %in% colnames(object@sample_metadata)) {
      
      clusters <- object@sample_metadata[,column]
      
      if(is.null(clusters)) {
        
        stop(paste0('error, defined column is null \n'))
        
      }
      
    }
    
  }
  
  if(!is.null(start.clus)) {
    
    if(!start.clus %in% clusters) {
      
      stop(paste0('start cluster is not present within the defined clusters \n'))
      
    }
    
  }
  
  if(!is.null(end.clus)) {
    
    if(!end.clus %in% clusters) {
      
      stop(paste0('end cluster is not present within the defined clusters \n'))
      
    }
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': initiating slingshot \n')))
  
  res <- slingshot::slingshot(data = red, clusterLabels = clusters, start.clus = start.clus, end.clus = end.clus, ...)
  
  pt <- slingshot::slingPseudotime(res)
  
  cat(crayon::cyan(paste0(Sys.time(), ': initiating slingshot \n')))
  
  slingres <- list(assignments = res, pseudotimes = pt)
  
  return(slingres)
  
}
