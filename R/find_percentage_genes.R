#' @name find_percentage_genes
#' @aliases find_percentage_genes
#' 
#' @title Finds the most variable genes in a dataset
#'
#' @description Find an x amount of top variable genes for your dataset
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
#' @export

find_percentage_genes <- function(object, 
                                  pattern='^MT-', 
                                  assay='RAW', 
                                  slot='counts',
                                  column.name = 'RAW_percent.mt') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(pattern)) {
    
    cat(crayon::cyan('pattern must be character string\n'))
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
  
  cat(crayon::cyan('Calculating percentage\n'))
  mat <- as.matrix(object@methods[[assay]][[slot]])
  subbed <- mat[grep(pattern = pattern, x = rownames(mat)),]
  temp <- Matrix::colSums(subbed) / Matrix::colSums(mat) * 100
  
  cat(crayon::cyan('Percentage calculated\n'))
  temp <- as.data.frame(temp)
  colnames(temp) <- column.name
  if(column.name %in% colnames(object@sample_metadata)) {
    cat(crayon::cyan('Removing old metadata column\n'))
    object@sample_metadata <- object@sample_metadata[,colnames(object@sample_metadata) != column.name]
  }
  cat(crayon::cyan('Appending new column\n'))
  colnames(temp) <- column.name
  temp <- apply(temp, 2, function(x) as.numeric(x))
  object@sample_metadata <- cbind(object@sample_metadata, temp)
  return(object)
}