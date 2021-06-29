#' @name perform.bbknn
#' @aliases perform.bbknn
#' 
#' @title Performs BBKNN integration
#'
#' @description Performs BBKNN integration on defined method-assays and reductions contained within. This is performed on reductions. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining the name of the reduction to provide for BBKNN. Default = NULL
#' @param graph.name Character. What should the BBKNN graph be named. Default = 'bbknn
#' @param batch Character. Column name in metadata indicating batch. Can be multiple.
#' @param approx Character. Employs annoy's approximate neighbour finding. Useful for large datasets but may increase correction. 
#' @param metric. Character. Which distance metric to use when approx is TRUE, options: 'angular', 'euclidean', 'manhattan' or 'hamming'. Default = 'euclidean'
#' @param neighbors_within_batch Numerical. How many neighbours to report per batch. Default = 3
#' @param n_pcs Numerical. Range of principal components to use. Default = NULL
#' @param trim Numerical. Trims the n of neighbours per cell to this value. Helps with population independence. Default = NULL
#' @param n_trees Numerical. Number of trees to generate in annoy forest. More trees provides higher precision at the cost of increased resource demand and run time. Default = 10
#' @param use_faiss Boolean. Uses faiss package to compute nearest neighbour, this improves run time at the cost of precision. Default = TRUE
#' @param set_op_mix_ratio Numerical. UMAP connectivity parameter between 0 and 1. controls the blen d between a connectivity matrix formed exclusively from mutual nearest neighbour pairs (0) and a union of all observed neighbour relationships with the mutual pairs emphasised (1). Default = 1.0
#' @param local_connectivity Numerical. How many nearest neighbours of each cell are assumed to be fully connected. Default = 1
#' 
#' @return BBKNN connectivity graph contained in graphs in the indicated method-assays
#'
#' @export

perform.bbknn <- function(object,
                          assay,
                          reduction,
                          graph.name = 'bbknn',
                          batch,
                          approx = FALSE,
                          metric = 'euclidean',
                          neighbors_within_batch = 3,
                          n_pcs = NULL,
                          trim = NULL,
                          n_trees = 10,
                          use_faiss = TRUE,
                          set_op_mix_ratio = 1.0,
                          local_connectivity= 1,
                          save.object = TRUE) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan('reduction must be character string\n'))
    return(object)
    
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
  
  if(!is.character(graph.name)) {
    
    cat(crayon::cyan('graph.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.character(batch)) {
    
    cat(crayon::cyan('batch must be character string\n'))
    return(object)
    
  }
  
  if(!batch %in% colnames(object@sample_metadata)) {
    
    cat(crayon::cyan('batch does not exist\n'))
    return(object)
    
  }
  
  if(!is.logical(approx)) {
    
    cat(crayon::cyan('approx must be logical: TRUE/FALSE \n'))
    return(object)
    
  }
  
  if(!is.character(metric)) {
    
    cat(crayon::cyan('metric must be character \n'))
    return(object)
    
  }
  
  if(!is.numeric(neighbors_within_batch)) {
    
    cat(crayon::cyan('neighbors_within_batch must be numerical \n'))
    return(object)
    
  }
  
  
  
  if(!is.null(trim)) {
    
    if(!is.numeric(trim)) {
      
      cat(crayon::cyan('trim must be character string\n'))
      return(object)
      
    }
    
  }
  
  if(!is.numeric(n_trees)) {
    
    cat(crayon::cyan('n_trees must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(use_faiss)) {
    
    cat(crayon::cyan('use_faiss must be logical: TRUE/FALSE\n'))
    
  }
  
  if(!is.numeric(set_op_mix_ratio)) {
    
    cat(crayon::cyan('set_op_mix_ratio must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(local_connectivity)) {
    
    cat(crayon::cyan('local_connectivity must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(save.object)) {
    
    cat(crayon::cyan('save.object must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  for(p in assay) {
    
    count <- 1
    
    for(r in reduction) {
      
      sc <- reticulate::import('scanpy')
      
      scobj <- sc$AnnData(X = object@methods[[p]]@computational_reductions[[reduction[count]]])
      
      scobj$obs_names <- as.factor(colnames(object))
      
      scobj$var_names <- as.factor(colnames(object@methods[[p]]@computational_reductions[[reduction[count]]]))
      
      scobj$obsm$update(X_pca = object@methods[[p]]@computational_reductions[[reduction[count]]])
      
      if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
        
        pd <- reticulate::import('pandas')
        
        scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))
        
      }
      
      if(is.null(n_pcs)) {
        
        cat(crayon::cyan('npcs calculated\n'))
        
        n_pcs <- as.integer(length(colnames(object@methods[[p]]@computational_reductions[[reduction[count]]])))
        
      }
      
      if(is.null(trim)) {
        
        cat(crayon::cyan(paste0('initialising BBKNN for assay: ', p,  ', reduction: ', r, '\n')))
        
        sc$external$pp$bbknn(scobj,
                             batch_key = as.character(batch),
                             approx = as.logical(FALSE),
                             metric = as.character(metric),
                             neighbors_within_batch = as.integer(neighbors_within_batch),
                             n_pcs = n_pcs,
                             n_trees = as.integer(n_trees),
                             use_faiss = as.logical(use_faiss),
                             set_op_mix_ratio = set_op_mix_ratio,
                             local_connectivity = local_connectivity)
        
        cat(crayon::cyan('BBKNN complete \n'))
        
      } else if (!is.null(trim)) {
        
        cat(crayon::cyan(paste0('initialising BBKNN for assay: ', p,  ', reduction: ', r, '\n')))
        
        sc$external$pp$bbknn(scobj,
                             batch_key= as.character(batch),
                             approx = as.logical(FALSE),
                             metric = as.character(metric),
                             neighbors_within_batch = as.integer(neighbors_within_batch),
                             n_pcs = n_pcs,
                             trim = as.integer(trim),
                             n_trees = as.integer(n_trees),
                             use_faiss = as.logical(use_faiss),
                             set_op_mix_ratio = set_op_mix_ratio,
                             local_connectivity = local_connectivity)
        
        cat(crayon::cyan('BBKNN complete \n'))
        
      }
      
      graph.list <- list()
      connectivities <- scobj$obsp[['connectivities']]
      colnames(connectivities) <- colnames(object)
      rownames(connectivities) <- colnames(object)
      distances <- scobj$obsp[['distances']]
      colnames(distances) <- colnames(object)
      rownames(distances) <- colnames(object)
      connectivities <- as(object = as.matrix(connectivities), Class = 'dgCMatrix')
      distances <- as(object = as.matrix(distances), Class = 'dgCMatrix')
      
      graph.list[['connectivities']] <- connectivities
      graph.list[['distances']] <- distances
      
      object@methods[[p]]@neighbours[[graph.name[count]]] <- graph.list

      count <- count + 1
      
    }
    
  }
  
  return(object)
  
}
