#' @name perform.sc3.slot.cluster
#' @aliases perform.sc3.slot.cluster
#' 
#' @title Performs SC3 clustering on matrix slot
#'
#' @description Performs SC3 clustering on defined method-assays slot.
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. Which slot within the assay should be supplied.
#' @param HVGs Boolean. Should previously the count matrix be subset to highly variable genes. 
#' @param assignment.df.name Character. What to call the df contained in clusters.
#' @param ks Numerical range. Number of clusters to identify, this can be a range, i.e. 5:10.
#' @param n.core Numerical. How many cores should be used to run SC3. Default = 3
#' 
#' @return Cluster assignments using the list of resolutions provided contained within cluster_assignments under assignment.df.name
#'
#' @export

perform.sc3.slot.cluster <- function(object, 
                                     assay,
                                     slot,
                                     HVGs,
                                     assignment.df.name,
                                     ks, 
                                     n.core = 3) {
  
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
  
  if(!is.character(slot)) {
    
    stop(paste0('slot must be character string\n'))
    
  }
  
  if(!is.logical(HVGs)) {
    
    cat(crayon::cyan(paste0('HVGs must be logical: TRUE/FALSE\n')))
    
  }
  
  if(!is.character(assignment.df.name)) {
    
    stop(paste0('assignment.df.name must be character string(s)\n'))
    
  }
  
  if(!is.numeric(ks)) {
    
    stop(paste0('ks must be numerical\n'))
    
  }
  
  if(!is.numeric(n.core)) {
    
    stop(paste0('n.core must be numerical\n'))
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': initialising SC3 clustering\n')))
  
  for(p in assay) {
    
    mat <- object@methods[[p]][[slot]]
    
    if(!is.null(HVGs)) {
      
      mat <- mat[object@methods[[p]]@highly.variable.genes,]
      
    }
    
    
    if(is.null(HVGs)) {
      
      temp.2 <- SingleCellExperiment(list('counts' = as.matrix(object@methods[[p]]@counts[rownames(mat),]), 'logcounts' = as.matrix(mat)))
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- sc3_prepare(temp.2, gene_filter = TRUE, n_cores = n.core)
      
    } else {
      
      temp.2 <- SingleCellExperiment(list('logcounts' = as.matrix(mat)))
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
      
    }
    
    temp.2 <- sc3_calc_dists(temp.2)
    temp.2 <- sc3_calc_transfs(temp.2)
    temp.2 <- sc3_kmeans(temp.2, ks = ks)
    temp.2 <- sc3_calc_consens(temp.2)
    object@methods[[p]]@cluster_assignments[[assignment.df.name]] <- as.data.frame(colData(temp.2))
    cat(crayon::cyan(paste0(Sys.time(), ': SC3 clustering completed\n')))
    
  }
  
  return(object)
  
}