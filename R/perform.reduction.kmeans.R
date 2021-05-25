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
                                     assignment.df.name,
                                     method='kmeans',
                                     ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(NULL)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s) \n'))
    return(NULL)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
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
  
  if(!is.numeric(k)) {
    
    cat(crayon::cyan(paste0('k must be numerical\n')))
    return(object)
    
  }
  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan(paste0('assignment.df.name must be character string(s)\n')))
    return(object)
    
  }
  
  if(!method %in% c('kmeans', 'pam')) {
    
    cat(crayon::cyan(paste0('method must be kmeans or pam\n')))
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
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        cat(crayon::cyan('reductions could not be found\n'))
        return(object)
        
      }
      
    }
    
    print(reduction)
    print(names(reduction.list))
    
    reduction.list <- reduction.list[reduction]
    
    count <- 1
    
    for(o in reduction) {
      
      if(length(reduction) == 1) {
        
        red <- reduction.list
        
      } else {
        
        red <- reduction.list[[o]] 
        
      }
      
      dimen <- dims[[count]]
      
      if(is.null(dimen)) {
        
        dimen <- 1:length(colnames(red))
        
      }
      
      clusters <- data.frame(kmeans=numeric(length(colnames(object))))
      rownames(clusters) <- colnames(object)
      
      if(is.null(k)) {
        cat(crayon::cyan('Specify number of clusters\n'))
      }
      if(is.null(reduction)) {
        cat(crayon::cyan('Provide assay\n'))
      }
      if(method == 'pam') {
        for(i in k) {
          clusters[,paste0('pam_clustering_K_', i)] <- as.factor(cluster::pam(x = red[,dimen], k = i, ...)$clustering)
        }
      }
      if(method == 'kmeans') {
        for(i in k) {
          clusters[[paste0('kmeans_clustering_K_', i)]] <- as.factor(kmeans(x = red[,dimen], centers = i, ...)$cluster)
        }
      } else {
        cat(crayon::cyan('Please specify method: pam or kmeans\n'))
      }
      
      clusters <- clusters[,2:length(colnames(clusters))]
      object@methods[[p]]@cluster_assignments[[assignment.df.name[[count]]]] <- clusters
      
      count <- count + 1
    }
    
  }
  
  return(object)
}