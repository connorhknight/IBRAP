#' @name perform.sct
#' @aliases perform.sct
#' 
#' @title Performs SCTransform
#'
#' @description Performs SCTransform normalisation, hvg selection, scaling and variance stabilisation and regression.  
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

perform.sct <- function(object, 
                        assay,
                        slot,
                        new.assay.name = 'SCT',
                        ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Object must be of class IBRAP\n')
    
  }
  if(!is.character(assay)) {
    
    stop('Assay must be a character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('assay does not exist\n')
    
  }
  
  if(!is.character(slot)) {
    
    stop('Slot must be a character string\n')
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    stop('slot does not exist\n')
    
  }
  
  if(!is.character(new.assay.name)) {
    
    stop('new.assay.name must be character string\n')
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': converting to Seurat object\n')))
  seuratobj <- Seurat::CreateSeuratObject(counts = as.matrix(object@methods[[assay]][[slot]]), project = 'NA')
  seuratobj@meta.data <- cbind(seuratobj@meta.data, object@sample_metadata)
  cat(crayon::cyan(paste0(Sys.time(), ': initiating SCTransform\n')))
  seuratobj <- Seurat::SCTransform(object = seuratobj, ...)
  
  .highly.variable.genes <- as.character(seuratobj@assays$SCT@var.features)
  .counts <- as(object = as.matrix(seuratobj@assays$SCT@counts), Class = 'dgCMatrix')
  .normalised <- as(as.matrix(seuratobj@assays$SCT@data), Class = 'dgCMatrix')
  .norm.scaled <- as.matrix(seuratobj@assays$SCT@scale.data)
  feat.meta <- feature_metadata(assay = as.matrix(.counts), col.prefix = new.assay.name)
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = as.matrix(.counts), col.prefix = new.assay.name))
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = .counts, 
                                          normalised = .normalised, 
                                          norm.scaled = .norm.scaled,
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feat.meta)
  cat(crayon::cyan(paste0(Sys.time(), ': SCT normalisation completed\n')))
  return(object)
}
