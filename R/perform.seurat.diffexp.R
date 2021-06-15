#' @name perform.seurat.diffexp
#' @aliases perform.seurat.diffexp
#' 
#' @title Plots slingshot results
#'
#' @description Plots the results of the slingshot analysis using ggplots2
#' 
#' @return A ggplot of the reduced cellular embbedings and trajectories. 
#'
#' @export
#' 
perform.seurat.diffexp <- function(object, 
                                   assay,
                                   test, 
                                   identity,
                                   latent.vars = NULL,
                                   ...) {
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
  seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
  
  seuobj$clusters <- identity
  Seurat::Idents(seuobj) <- 'clusters'
  
  if(!is.null(latent.vars)) {
    
    if(!is.character(latent.vars)) {
      
      cat(crayon::cyan('latent.vars must be character(s)\n'))
      return(NULL)
      
    }
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, latent.vars = latent.vars, ...)
    
  } else {
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, ...)
    
  }
  
  return(results)
  
}
