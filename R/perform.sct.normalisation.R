#' @name perform.sct.normalisation
#' @aliases perform.sct.normalisation
#' 
#' @title Performs SCTransform
#'
#' @description A new method-assay is produced. Raw counts are normalised and HVGs identified using SCTransform 
#' 
#' @param object IBRAP S4 class object
#' @param assay A character string containing indicating which assay to use
#' @param slot String indicating which slot within the assay should be sourced
#' @param new.assay.name Character. Name of new method-assay to be produced. Default = 'SCT'
#' @param do.scale Whether to scale residuals to have unit variance; default is FALSE
#' @param do.center Whether to center residuals to have mean zero; default is TRUE
#' @param n.genes Numerical value of how many highly variable genes should be retained. Default = 1500
#' @param min_cells Numerical value of minimum cells required for a gene to not be filtered. Default = 3
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#'
#' @export

perform.sct.normalisation <- function(object, 
                                      assay,
                                      slot,
                                      new.assay.name = 'SCT',
                                      do.scale = TRUE,
                                      do.center = TRUE,
                                      n.genes = 1500,
                                      min_cells = 3,
                                      save.seuratobject = TRUE,
                                      vars.to.regress = NULL,
                                      ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be of class IBRAP\n'))
    return(object)
    
  }
  if(!is.character(assay)) {
    
    cat(crayon::cyan('Assay must be a character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('Slot must be a character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(new.assay.name)) {
    
    cat(crayon::cyan('new.assay.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.logical(do.scale)) {
    
    cat(crayon::cyan('do.scale must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.logical(do.center)) {
    
    cat(crayon::cyan('do.center must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.numeric(n.genes)) {
    
    cat(crayon::cyan('n.genes must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(min_cells)) {
    
    cat(crayon::cyan('min_cells must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(save.seuratobject)) {
    
    cat(crayon::cyan('save.seuratobject must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  cat(crayon::cyan('Converting to Seurat object\n'))
  seuratobj <- Seurat::CreateSeuratObject(counts = as.matrix(object@methods[[assay]][[slot]]), project = 'NA')
  seuratobj@meta.data <- cbind(seuratobj@meta.data, object@sample_metadata)
  cat(crayon::cyan('Initiating SCTransform\n'))
  seuratobj <- Seurat::SCTransform(object = seuratobj, do.scale = do.scale, do.center = do.center, vars.to.regress = vars.to.regress, min_cells = min_cells, variable.features.n = n.genes, ...)
  cat(crayon::cyan('SCTransform completed!\n'))
  .highly.variable.genes <- as.character(seuratobj@assays$SCT@var.features)
  .counts <- as(object = as.matrix(seuratobj@assays$SCT@counts), Class = 'dgCMatrix')
  .normalised <- as(as.matrix(seuratobj@assays$SCT@data), Class = 'dgCMatrix')
  .norm.scaled <- as.matrix(seuratobj@assays$SCT@scale.data)
  feat.meta <- feature_metadata(assay = as.matrix(.counts), col.prefix = new.assay.name)
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = .counts, 
                                          normalised = .normalised, 
                                          norm.scaled = .norm.scaled,
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feat.meta)
  if(isTRUE(save.seuratobject)) {
    
    object@methods[[new.assay.name]]@alt_objects[['seurat']] <- seuratobj
    
  }
  cat(crayon::cyan('Populated IBRAP object\n'))
  return(object)
}
