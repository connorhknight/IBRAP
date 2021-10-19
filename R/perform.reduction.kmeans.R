#' @name perform.reduction.kmeans
#' @aliases perform.reduction.kmeans
#' 
#' @title Performs Kmeans/PAM clustering on a reduction
#'
#' @description Performs Kmeans/PAM clustering on defined method-assays a defined reduction.
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. Which reduction(s) within the assay should be supplied. Default = NULL
#' @param dims Numerical. How many dimensions of the reduction should be supplied.Default = NULL
#' @param k Numerical. How many clusters should the algorithm identify, this can be a range. Default = NULL
#' @param assignment.df.name Character. What to call the df contained in clusters.
#' @param method Which algorithm should be used, options: 'kmeans', 'pam'. Default = 'kmeans
#' @param ... Arguments to be passed to either base::kmeans or cluster::pam, depending on the defined method
#' 
#' @return Cluster assignments using the list of resolutions provided contained within cluster_assignments under assignment.df.name
#'
#' @export

perform.reduction.kmeans <- function(object, 
                                     assay,
                                     reduction=NULL, 
                                     dims = NULL,
                                     k=NULL,
                                     cluster.df.name.suffix='',
                                     method='kmeans',
                                     ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string(s) \n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('reduction: ', x, 'does not exist\n'))
      
    }
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        stop(paste0('reduction: ', x, ' does not exist\n'))
        
      }
      
    }
    
  }
  
  if(!is.numeric(k)) {
    
    stop(paste0('k must be numerical\n'))
    
  }
  
  if(!is.null(cluster.df.name.suffix)) {
    
    if(!is.character(cluster.df.name.suffix)) {
      
      stop(paste0('cluster.df.name.suffix must be character string(s)\n'))
      
    }
    
  }
  
  if(!method %in% c('kmeans', 'pam')) {
    
    stop(paste0('method must be kmeans or pam\n'))
    
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
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
      stop('reductions could not be found\n')
        
      }
      
    }
    
    reduction.list <- reduction.list[reduction]
    
    count <- 1
    
    for(o in reduction) {
      
      red <- reduction.list[[o]] 
      
      dimen <- dims[[count]]
      
      if(is.null(dimen)) {
        
        dimen <- 1:ncol(red)
        
      }
      
      clusters <- data.frame(kmeans=numeric(length(colnames(object))))
      rownames(clusters) <- colnames(object)
      
      if(is.null(k)) {
        stop('specify number of clusters\n')
      }
      if(is.null(reduction)) {
        stop('provide assay\n')
      }
      if(method == 'pam') {
        for(i in k) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': calculating PAM for k = ', i, '\n')))

          clusters[,paste0('k_', i)] <- as.factor(cluster::pam(x = red[,dimen], k = i, ...)$clustering)
          
        }
      }
      if(method == 'kmeans') {
        for(i in k) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': calculating kmeans for k = ', i, '\n')))

          clusters[[paste0('k_', i)]] <- as.factor(kmeans(x = red[,dimen], centers = i, ...)$cluster)
        }
      } 
      
      clusters <- clusters[,2:length(colnames(clusters))]
      
      if('_' %in% unlist(strsplit(x = cluster.df.name.suffix, split = ''))) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in cluster.df.name.suffix, replacing with - \n')))
        cluster.df.name.suffix <- sub(pattern = '_', replacement = '-', x = cluster.df.name.suffix)
        
      }
      
      if(method == 'pam') {
        
        object@methods[[p]]@cluster_assignments[[paste0(o, ':pam', cluster.df.name.suffix)]] <- clusters
        
      }
      
      if(method == 'kmeans') {
        
        object@methods[[p]]@cluster_assignments[[paste0(o, ':kmeans', cluster.df.name.suffix)]] <- clusters
        
      }
      
      count <- count + 1
    }
    
  }
  
  if(method == 'pam') {
    
    cat(crayon::cyan(paste0(Sys.time(), ': finished PAM clustering \n')))
    
  }
  
  if(method == 'kmeans') {
    
    cat(crayon::cyan(paste0(Sys.time(), ': finished kmeans clustering \n')))
    
  }
  
  return(object)
}