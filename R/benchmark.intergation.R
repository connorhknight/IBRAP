#' @name benchmark.integration
#' @aliases benchmark.integration
#' 
#' @title Benchmarks the cluster assignmentws
#'
#' @description Adopting scones method of Average Silhouette Width determination between batches. A depleted value is superior. 
#' 
#' @param object IBRAP S4 class object
#' @param batch Character. Whcih column in the metadata dataframe contains the batches.
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. Which reduction(s) should be used to observe integration.
#' @param result.names Character. A vector of names for the benchmarking results. 
#' @param n.components Numerical. How many components of the reduced embeddings should be used. Default = 2
#' 
#' @return Benchmarking scores for the supplied integrations
#' 
#' @examples object <- benchmark.integration(object = object, 
#'                                           batch = 'original.project', 
#'                                           assays = c('SCT','SCRAN','SCANPY'), 
#'                                           reduction = c('pca_umap', 'pca_harmony_umap', 
#'                                                         'scanorama_umap', 'pca_bbknn_bbknn:umap',
#'                                                         'CCA_pca_umap'), 
#'                                           result.names = c('uncorrected', 'harmony', 
#'                                                            'scanorama', 'bbknn', 'cca'), 
#'                                           n.components = 2)
#'
#' @export

benchmark.integration <- function(object, batch, assays, reduction, result.names, n.components = 2) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(batch)) {
    
    stop('batch must be character string \n')
    
  } else {
    
    if(!batch %in% colnames(object@sample_metadata)) {
      
      stop('batch cannot be found in object@sample_metadata \n')
      
    }
    
  }
  
  if(!is.character(assays)) {
    
    stop('batch must be character string \n')
    
  } else {
    
    for(x in assays) {
      
      if(!x %in% names(object@methods)) {
        
        stop(paste0(x, ' cannot be found in object@methods \n'))
        
      }
      
    }
    
  }
  
  if(!is.numeric(n.components)) {
    
    stop('n.components must be numerical \n')
    
  }
  
  if(!is.character(result.names)) {
    
    stop('result.names must be charcter string(s) \n')
    
  } else {
    
    if(length(result.names) != length(reduction)) {
      
      stop('not enough result.names provided \n')
      
    }
    
  }
  
  for(p in assays) {
    
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
    
    counter <- 1
    
    for(t in reduction) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': benchmarking assay: ', p, ', reduction: ', t, '\n')))
      
      dist <- as.matrix(dist(reduction.list[[t]][, seq_len(n.components)]))
      
      if(isTRUE(anyNA(as.numeric(object@sample_metadata[,batch])))) {
        
        count <- 1
        
        batch_tmp <- object@sample_metadata[,batch]
        
        for(d in unique(object@sample_metadata[,batch])) {
          
          batch_tmp[batch_tmp==d] <- count
          
          count <- count + 1
          
        }
        
        batch_tmp <- as.numeric(batch_tmp)
        
      } else {
        
        batch_tmp <- as.numeric(object@sample_metadata[,batch])
        
      }
      
      object@methods[[p]]@benchmark_results[['integration']][[paste0(p, '_', result.names[counter])]] <- summary(cluster::silhouette(as.numeric(batch_tmp), dist))$avg.width

      counter <- counter + 1
      
    }
    
  }
  
  return(object)
  
}
