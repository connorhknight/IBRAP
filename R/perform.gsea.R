#' @name perform.gsea
#' @aliases perform.gsea
#' 
#' @title Returnas gene set enrichment scores for indivdiual cells
#'
#' @description calculates gene set enrichment scores
#' 
#' @param object IBRAP S4 class object
#' @param gene.sets Either a GeneSetCollection derived from escape::getGeneSets() or a list of genes with their grouping as the list item name. Default = escape::getGeneSets(library = 'H') 
#' @param groups Numerical. The number of groups the cells should be split into. Default = 1000.
#' @param cores Numerical. The number of cores to commit to the calculations.
#' @param return_object Logical. Whether to return an updated IBRAP object 
#' 
#' @usage perform.gsea(object = obj)
#' 
#' @return IBRAP S4 class object containing module scores for each cell in the metadata
#'
#' @export

perform.gsea <- function(object, assay='RAW', slot='counts', gene.sets=escape::getGeneSets(library = 'H'), groups = 1000, cores = 1, return_object=T) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('Assay must be character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('Assay not contained in the supplied IBRAP object \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('Slot must be character string\n')
    
  }
  
  if(!assay %in% c('counts','normalised','norm.scaled')) {
    
    stop('Slot must either be: counts, normalised or norm.scaled. However, counts is highly recommended! \n')
    
  }
  
  if(isFALSE(is(gene.sets, 'GeneSetCollection')) && isFALSE(is.list(gene.sets))) {
    
    stop('gene.sets must be either GeneSetCollection or list object \n')
    
  }
  
  if(!is.numeric(cores)) {
    
    stop('cores must be a numerical value \n')
    
  }
  
  if(!is.numeric(cores)) {
    
    stop('groups must be a numerical value \n')
    
  }
  
  if(!is.logical(return_object)) {
    
    stop('return_object must be a logical value \n')
    
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': initiating gene set enrichment analysis \n')))

  enriched <- enrichIt(obj = as_matrix(object@methods[[1]]@counts), gene.sets = gene.sets, groups = groups, cores = cores)
  
  if(isTRUE(return_object)) {
    
    for(o in names(enriched)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
          
        }
        
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    
    object@sample_metadata <- cbind(object@sample_metadata, enriched)
    
    cat(crayon::cyan(paste0(Sys.time(), ': gene set enrichment analysis complete \n')))
    
    return(object)
    
  } else {
    
    return(enriched)
    
  }
  
}