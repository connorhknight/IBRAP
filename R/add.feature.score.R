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
#' @param verbose Logical should function messages be printed?
#' @param seed Numerical What seed should be set. Default = 1234
#' 
#' @usage add.cell.cycle(object = obj, assay = 'RAW', slot = 'counts')
#' 
#' @return IBRAP S4 class object containing module scores for each cell in the metadata
#' 
#' @examples # object <- add.feature.score(object = object, 
#                                          assay = 'RAW', 
#                                          slot = 'counts',
#                                          transform = TRUE, 
#                                          features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
#                                                      'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
#                                                      'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
#                                                      'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
#                                                      'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
#                                          column.name = 'StressScore')
#'
#' @export

add.feature.score <- function(object, 
                              assay, 
                              slot,
                              transform, 
                              features, 
                              column.name,
                              verbose=FALSE,
                              seed=1234,
                              ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop(crayon::cyan('object must be of class IBRAP\n'))
    
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
    
    stop('slot does not exist\n')
    
  }
  
  if(!is.logical(transform)) {
    
    stop('transform must be logical: TRUE/FALSEt\n')
    
  }
  
  if(!is.character(features)) {
    
    stop('features must be character string(s)\n')
    
  }
  
  if(!is.character(column.name)) {
    
    stop('column.name must be character string\n')
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose must be logical: TRUE/FALSE\n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  genes <- rownames(object)
  
  genes <- list(genes[genes %in% features])
  
  if(transform == TRUE) {
    
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]])
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': converted to Seurat object\n')))
      
    }

    seuobj <- Seurat::NormalizeData(object = seuobj, verbose = verbose)
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': data transformed\n')))
      
    }

    seuobj <- Seurat::AddModuleScore(object = seuobj, features = genes, ...)
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': seurat gene score calculated\n')))
      
    }

    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
          
        }

        object@sample_metadata[,o] <- NULL
        
        
      }
      
    }
    
    object@sample_metadata[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': new metadata added\n')))
      
    }

  } else {
    
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][['counts']])
    
    seuobj@assays$RNA@data <- object@methods[[assay]][[slot]]
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': converted to Seurat object\n')))
      
    }

    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, verbose = verbose, ...)
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': seurat gene score calculated\n')))
      
    }

    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': found duplicated column name: ',o, ' removing old column names.\n')))
          
        }
        
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    
    object@sample_metadata[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': new metadata added\n')))
      
    }

  }
  
  return(object)
}