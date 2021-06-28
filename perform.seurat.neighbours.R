

perform.seurat.neighbours <- function(object, assay, reduction, dims=NULL) {
  
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
      
      
    }
    
  }

}
