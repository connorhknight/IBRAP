#' @name add.cell.cycle
#' @aliases add.cell.cycle
#' 
#' @title Scores cell cycle phases
#'
#' @description Scores each cell in which cell cycle stage it is currently in. 
#' 
#' @param object IBRAP S4 class object
#' @param assay A character string containing indicating which assay to use
#' @param slot String indicating which slot within the assay should be sourced
#' @param transform Boolean. If raw counts are supplied, this must be TRUE to normalise data
#' 
#' @usage add.cell.cycle(object = obj, assay = 'RAW', slot = 'counts')
#' 
#' @return IBRAP S4 class object containing cell cycle assignments and scores for each cell in the metadata
#'
#' @export

add.cell.cycle <- function(object, 
                           assay,
                           slot,
                           transform, ...) {
  
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
    
    cat(crayon::cyan('slot must be character string\n'))
    return(object)
    
  }
  
  if(!is.logical(transform)) {
    
    cat(crayon::cyan('transform must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  r <- read.table(text = as.character(Homo_sapiens$phase.geneID.GeneName), sep = ',')
  colnames(r) <- c('phase', 'geneID', 'geneName')
  
  cat(crayon::cyan('Cell cycle genes loaded\n'))
  
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]])
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Cell cycle scores identified\n'))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    df <- seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))]
    object@sample_metadata <- cbind(object@sample_metadata, df)
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][['counts']])
    seuobj@assays$RNA@data <- object@methods[[assay]][[slot]]
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Data transformed\n'))
    
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    
    df <- seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))]
    object@sample_metadata <- cbind(object@sample_metadata, df)
    cat(crayon::cyan('New metadata added\n'))
    
  }
  return(object)
}