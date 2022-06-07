#' @name Read10X_output
#' @aliases Read10X_output
#' 
#' @title Produce counts matrix from CellRanger output
#'
#' @description Concatenates the three output files of CellRanger alignment.
#' 
#' @param directory Absolute pathway to directory containing the three files
#' @param matrix.file Specify .mtx file name. Default = matrix.mtx
#' @param genes.file Specify file name for genes. Default = genes.tsv
#' @param barcodes.file Specify file name for barcodes. Default = barcodes.tsv
#' @param make.feat.unique Should non-unique gene names be made unique. 
#' @param verbose Logical Should function messages be printed?
#' 
#' @usage Read10X_output(directory, matrix.file, genes.file, barcodes.file, make.feat.unique)
#' 
#' @examples Read10X_output(directory = 'path/to/folder/', matrix.file = 'matrix.mtx', genes.file = 'genes.tsv', barcodes.file = 'barcodes.tsv', make.feat.unique = TRUE)
#'
#' @return Sparse matrix of counts
#'
#' @export

Read10X_output <- function(directory, 
                           matrix.file = 'matrix.mtx', 
                           genes.file = 'genes.tsv', 
                           barcodes.file = 'barcodes.tsv',
                           make.feat.unique = TRUE,
                           verbose = FALSE) {
  
  if(!is.character(directory)) {
    
    stop('directory should be a character string /n')
    
  }
  
  if(!is.character(matrix.file)) {
    
    stop('matrix.file should be a character string /n')
    
  }
  
  if(!is.character(genes.file)) {
    
    stop('genes.file should be a character string /n')
    
  }
  
  if(!is.character(barcodes.file)) {
    
    stop('barcodes.file should be a character string /n')
    
  }
  
  if(!is.logical(make.feat.unique)) {
    
    stop('make.feat.unique should be a logical. TRUE/FALSE /n')
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be a logical. TRUE/FALSE /n')
    
  }
  
  if(isFALSE(dir.exists(directory))) {
    
    stop('directory does not exist /n')
    
  }
  
  dir.files <- list.files(path = directory)
  
  if(!matrix.file %in% dir.files) {
    cat(crayon::cyan(paste0(Sys.time(), ': expected file: ', matrix.file, '\n')))
  }
  
  if(!genes.file %in% dir.files) {
    cat(crayon::cyan(paste0(Sys.time(), ': expected file: ', genes.file, '\n')))
  }
  
  if(!barcodes.file %in% dir.files) {
    cat(crayon::cyan(paste0(Sys.time(), ': expected file: ', barcodes.file, '\n')))
  }
  
  count.file <- paste0(directory, '/', matrix.file)
  
  genes.file <- paste0(directory, '/', genes.file)
  
  barcodes.files <- paste0(directory, '/', barcodes.file)
  
  genes_ensembl <- utils::read.table(file = genes.file, sep = '\t', header = FALSE)
  
  barcodes <- utils::read.table(file = barcodes.files, sep = '\t', header = FALSE)
  
  mm <- Matrix::readMM(count.file)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': files loaded\n')))
    
  }

  true.length <- length(genes_ensembl$V2)
  
  if(true.length != length(unique(genes_ensembl$V2))) {
    
    if(isTRUE(make.feat.unique)) {
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': non-unique features identified\n')))
        
      }
      
      genes_ensembl$V2 <- make.unique(genes_ensembl$V2)
      
    }
    
  }
  
  rownames(mm) <- genes_ensembl$V2
  
  colnames(mm) <- barcodes$V1
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': success: Matrix concatenated\n')))
    
  }
  
  mm <- as_matrix(as(object = mm, Class = 'dgCMatrix'))
  
  return(mm)
}