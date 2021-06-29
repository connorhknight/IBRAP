

perform.seurat.neighbours <- function(object, 
                                      assay, 
                                      reduction, 
                                      dims=NULL,
                                      neighbour.name=NULL,
                                      k.param=20,
                                      compute.SNN=TRUE,
                                      prune.SNN=1/15, 
                                      nn.method='annoy', 
                                      n.trees = 50,
                                      nn.eps=0.0, 
                                      annoy.metric='euclidean') {
  
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

    seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
    
    count <- 1
    
    for(g in reduction) {
      
      dim <- dims[[count]]
      red <- reduction.list[[g]]
      red.key <- reduction[count]
      
      seuobj@reductions[[red.key]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(red), key = red.key)
      
      
      seuobj <- suppressWarnings(Seurat::FindNeighbors(object = seuobj, 
                                                       k.param = k.param,
                                                       reduction = red.key, 
                                                       verbose = T, 
                                                       dims = dim, 
                                                       compute.SNN = compute.SNN, 
                                                       prune.SNN = prune.SNN,
                                                       nn.method = nn.method, 
                                                       n.trees = n.trees,
                                                       annoy.metric = annoy.metric, 
                                                       nn.eps = nn.eps))
      
      object@methods[[p]]@neighbours[[neighbour.name]][['connectivities']] <- seuobj@graphs$RNA_snn
      
      count <- count + 1
      
    }
    
  }
  
  return(object)

}
