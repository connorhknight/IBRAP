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
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('assay does not exist\n')
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be character string\n')
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    stop('slot must be character string\n')
    
  }
  
  if(!is.logical(transform)) {
    
    stop('transform must be logical: TRUE/FALSE\n')
    
  }
  
  r <- read.table(text = as.character(IBRAP::Homo_sapiens$phase.geneID.GeneName), sep = ',')
  colnames(r) <- c('phase', 'geneID', 'geneName')
  
  cat(crayon::cyan(paste0(Sys.time(), ': cell cycle genes loaded\n')))
  
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]])
    
    
    cat(crayon::cyan(paste0(Sys.time(), ':converted to Seurat object\n')))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan(paste0(Sys.time(), ': data transformed\n')))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan(paste0(Sys.time(), ': cell cycle scores identified\n')))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    df <- seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))]
    object@sample_metadata <- cbind(object@sample_metadata, df)
    cat(crayon::cyan(paste0(Sys.time(), ': new metadata added\n')))
  } else {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][['counts']])
    seuobj@assays$RNA@data <- object@methods[[assay]][[slot]]
    cat(crayon::cyan(paste0(Sys.time(), ': converted to Seurat object\n')))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan(paste0(Sys.time(), ': data transformed\n')))
    
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    
    df <- seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))]
    object@sample_metadata <- cbind(object@sample_metadata, df)
    cat(crayon::cyan(paste0(Sys.time(), ': new metadata added\n')))
    
  }
  return(object)
}