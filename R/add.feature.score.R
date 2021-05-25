#' @name add.feature.score
#' @aliases add.feature.score
#' 
#' @title Provides scores for a given vector of features
#'
#' @description Produces a module score per cell for the supplied genes
#' 
#' @param object IBRAP S4 class object
#' @param assay A character string containing indicating which assay to use
#' @param slot String indicating which slot within the assay should be sourced
#' @param transform Boolean. If raw counts are supplied, this must be TRUE to normalise data
#' @param features A character vector of genes to be scored
#' @param column.name Character naming the column containing the scores in the metadata dataframe
#' 
#' @usage add.cell.cycle(object = obj, assay = 'RAW', slot = 'counts')
#' 
#' @return IBRAP S4 class object containing module scores for each cell in the metadata
#'
#' @export

add.feature.score <- function(object, 
                              assay, 
                              slot,
                              transform, 
                              features, 
                              column.name,
                              ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.logical(transform)) {
    
    cat(crayon::cyan('transform must be logical: TRUE/FALSEt\n'))
    return(object)
    
  }
  
  if(!is.character(features)) {
    
    cat(crayon::cyan('features must be character string(s)\n'))
    return(object)
    
  }
  
  if(!is.character(column.name)) {
    
    cat(crayon::cyan('column.name must be character string\n'))
    return(object)
    
  }
  
  genes <- rownames(object)
  genes <- list(genes[genes %in% features])
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]])
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = genes, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    object@sample_metadata[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][['counts']])
    seuobj@assays$RNA@data <- object@methods[[assay]][[slot]]
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    object@sample_metadata[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}