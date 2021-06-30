#' @name perform.scanpy.neighbours
#'
#' @export

perform.scanpy.neighbours <- function(object, 
                                      assay,
                                      reduction,
                                      neighbour.name,
                                      n_neighbors = 15, 
                                      n_pcs = 0, 
                                      random_state = 0, 
                                      method = 'umap', 
                                      metric='euclidean') {
  
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
    
    sc <- reticulate::import('scanpy')
    
    count <- 1
    
    for(g in reduction) {
      
      scobj <- sc$AnnData(X = t(as.matrix(object@methods[[2]]@norm.scaled)))
      scobj$obs_names <- as.factor(colnames(object@methods[[2]]@norm.scaled))
      scobj$var_names <- as.factor(rownames(object@methods[[2]]@norm.scaled))
      
      if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
        
        pd <- reticulate::import('pandas')
        
        scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))
        
      }
      
      red <- reduction.list[[g]]
      red.key <- reduction[count]
      
      scobj$obsm$update(X_pca = as.matrix(red))
      
      sc$pp$neighbors(adata = scobj, 
                      n_neighbors = as.integer(n_neighbors), 
                      n_pcs = as.integer(n_pcs),
                      random_state = as.integer(random_state), 
                      method = as.character(method), 
                      metric = as.character(metric))
      
      object@methods[[p]]@neighbours[[neighbour.name]][['connectivities']] <- scobj$obsp[['connectivities']]
      
      object@methods[[p]]@neighbours[[neighbour.name]][['distances']] <- scobj$obsp[['distances']]
     
      count <- count + 1
       
    }
    
  }
  
  return(object)
  
}
