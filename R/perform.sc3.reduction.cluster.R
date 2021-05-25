#' @name perform.sc3.reduction.cluster
#' @aliases perform.sc3.reduction.cluster
#' 
#' @title Performs SC3 clustering on reduced embeddings
#'
#' @description Performs SC3 clustering on defined method-assays and supplied reductions. SC3 is designed to function on count matrices, reduction may not always work well.
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining which reduction to supply to the clustering algorithm.
#' @param dims Numerical. How many dimensions of the reduciton should be supplied, NULL equates to all.
#' @param assignment.df.name Character. What to call the df contained in clusters.
#' @param ks Numerical range. Number of clusters to identify, this can be a range, i.e. 5:10.
#' @param n.core Numerical. How many cores should be used to run SC3. Default = 3
#' 
#' @return Cluster assignments using the list of resolutions provided contained within cluster_assignments under assignment.df.name
#'
#' @export

perform.sc3.reduction.cluster <- function(object, 
                                          assay,
                                          reduction,
                                          dims,
                                          assignment.df.name,
                                          ks, 
                                          n.core = 3) {
  
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
  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan(paste0('assignment.df.name must be character string(s)\n')))
    return(object)
    
  }
  
  if(!is.numeric(ks)) {
    
    cat(crayon::cyan(paste0('ks must be numerical\n')))
    return(object)
    
  }
  
  if(!is.numeric(n.core)) {
    
    cat(crayon::cyan(paste0('n.core must be numerical\n')))
    return(object)
    
  }
  
  cat(crayon::cyan('Initialising SC3 clustering\n'))
  
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
    
    reduction.list <- reduction.list[reduction]
    
    count <- 1
    
    for(r in reduction) {
      
      red <- reduction.list[[r]]
      
      dimen <- dims[[count]]
      
      if(is.null(dimen)) {
        
        dimen <- 1:ncol(red)
        
      }
      
      temp.2 <- SingleCellExperiment::SingleCellExperiment(list('logcounts' = t(red)[dimen,]))
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- SC3::sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
      temp.2 <- SC3::sc3_calc_dists(temp.2)
      temp.2 <- SC3::sc3_calc_transfs(temp.2)
      temp.2 <- SC3::sc3_kmeans(temp.2, ks = ks)
      temp.2 <- SC3::sc3_calc_consens(temp.2)
      object@methods[[p]]@cluster_assignments[[assignment.df.name[[count]]]] <- as.data.frame(colData(temp.2))
      cat(crayon::cyan('SC3 clustering completed\n'))
      count <- count + 1
      
    }
    
  }
  
  return(object)
  
}