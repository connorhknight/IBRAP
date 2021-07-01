#' @name perform.singleR.annotation
#' @aliases perform.singleR.annotation
#' 
#' @title Performs automated cell annotation on query datasets using reference data
#'
#' @description SingleR iterates through singular cells and iterates through probabilitiy comparisons to identify which cell type the query cell is likely to be. If a probably cell type cannot be discovered then 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param ref AnyMatrix. A matrix of the reference datasets, if data is end-bias then it should be log normalised, if it is full-length then it requires tpm normalisation. Both can be completed within this function. 
#' @param log.transform Boolean. Should the reference data be log transformed. Default = TRUE
#' @param tpm.transform Boolean. Should the reference data be tpm normalised. Default = FALSE
#' @param ref.labels Vector. The cluster assignments for the reference data. Default = NULL
#' @param column.suffix Character. A suffix to append the end of the new metadata columns if this functiuons is to be used multiple times. Default = '1'
#' @param ... arguments to be passed to singleR::SingleR
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#'
#' @export

perform.singleR.annotation <- function(object, 
                                       assay = 'RAW', 
                                       slot = 'counts', 
                                       ref, 
                                       log.transform = TRUE, 
                                       tpm.transform = FALSE, 
                                       ref.labels, 
                                       column.suffix='1', 
                                       ...) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object should be IBRAP class \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be a character string \n'
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be a character string \n')
    
  }
  
  if(!is(ref, 'matrix')) {
    
    stop('reference matrix must be a matrix class \n')
    
  }
  
  if(!is.logical(log.transform)) {
    
    stop('log.transform must be logical, TRUE/FALSE \n'
    
  }
  
  
  if(!is.logical(tpm.transform)) {
    
    stop('tpm.transform must be logical, TRUE/FALSE \n')
    
  }
  
  if(!is.vector(ref.labels)) {
    
    stop('ref.labels must be vector, TRUE/FALSE \n')
    
  }
  if(!is.character(column.suffix)) {
    
    stop('column.suffix must be character string \n')
    
  }
  
  if(isTRUE(tpm.transform)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': tpm transforming reference data \n')))
    
    temp <- createIBRAPobject(counts = ref, original.project = 'tpm_transform')
    temp <- perform.tpm.normalisation(object = temp)
    ref <- temp@methods$TPM@normalised
    
  } else if (isTRUE(log.transform) && isFALSE(tpm.transform)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': log2 transforming reference data \n')))
    
    ref <- log2(ref+1)
    
  }
  
  query <- object@methods[[assay]][[slot]]
  
  cat(crayon::cyan(paste0(Sys.time(), ': initiating singleR automated labelling \n')))
  
  result <- SingleR::SingleR(test = query, ref = ref, labels = ref.labels, ...)
  
  cat(crayon::cyan(paste0(Sys.time(), ': completed singleR automated labelling \n')))
  
  result$pruned.labels[is.na(result$pruned.labels)] <- 'X'
  result$labels[is.na(result$labels)] <- 'X'
  
  object@sample_metadata[,paste0('singleR_pruned_labels_', column.suffix)] <- result$pruned.labels
  object@sample_metadata[,paste0('singleR_labels_', column.suffix)] <- result$labels
  
  cat(crayon::cyan(paste0(Sys.time(), ': appending annotations to metadata \n')))
  cat(crayon::cyan(paste0(Sys.time(), ': complete \n')))
  
  return(object)
  
}
