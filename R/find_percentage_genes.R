#' @name find_percentage_genes
#' @aliases find_percentage_genes
#' 
#' @title Calculates the fraction of counts from genes matching a pattern string
#'
#' @description Subsets gene names that match the pattern supplied. The percentage fraction of this gene group is then calculated.
#' 
#' @param object IBRAP S4 class object
#' @param pattern A character string containing a pattern to identfiy in rownames
#' @param assay String indicating which assay to source the raw counts from
#' @param slot String indicating which slot within the assay should be sourced
#' @param column.name String naming the column name in the metadata
#' 
#' @usage find_percentage_genes(object = obj, pattern = '^MT-', assay = 'RAW', slot = 'RAW', column.name = 'RAW_percent.mt')
#' 
#' @return IBRAP S4 class object containing highly varaible genes within the source assay
#' 
#' @examples 
#' 
#' object <- find_percentage_genes(object = object, pattern = '^MT-',
#                                  assay = 'RAW', slot = 'counts',
#                                  column.name = 'RAW_percent.mt')
#'
#' @export

find_percentage_genes <- function(object, 
                                  pattern='^MT-', 
                                  assay='RAW', 
                                  slot='counts',
                                  column.name = 'RAW_percent.mt') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(pattern)) {
    
    stop('pattern must be character string\n')
    
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

  cat(crayon::cyan(paste0(Sys.time(), ': calculating percentage\n')))
  mat <- as.matrix(object@methods[[assay]][[slot]])
  subbed <- mat[grep(pattern = pattern, x = rownames(mat)),]
  temp <- Matrix::colSums(subbed) / Matrix::colSums(mat) * 100
  
  cat(crayon::cyan(paste0(Sys.time(), ': percentage calculated\n')))
  temp <- as.data.frame(temp)
  colnames(temp) <- column.name
  if(column.name %in% colnames(object@sample_metadata)) {
    cat(crayon::cyan(paste0(Sys.time(), ': removing old metadata column\n')))
    object@sample_metadata <- object@sample_metadata[,colnames(object@sample_metadata) != column.name]
  }
  cat(crayon::cyan(paste0(Sys.time(), ': appending new column\n')))
  colnames(temp) <- column.name
  temp <- apply(temp, 2, function(x) as.numeric(x))
  object@sample_metadata <- cbind(object@sample_metadata, temp)
  return(object)
}