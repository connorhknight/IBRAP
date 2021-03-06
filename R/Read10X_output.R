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
                           make.feat.unique = TRUE) {
  
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
  cat(crayon::cyan(paste0(Sys.time(), ': files loaded\n')))
  true.length <- length(genes_ensembl$V2)
  if(true.length != length(unique(genes_ensembl$V2))) {
    if(isTRUE(make.feat.unique)) {
      cat(crayon::cyan(paste0(Sys.time(), ': non-unique features identified\n')))
      genes_ensembl$V2 <- make.unique(genes_ensembl$V2)
    }
    
  }
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  cat(crayon::cyan(paste0(Sys.time(), ': success: Matrix concatenated\n')))
  mm <- as.matrix(mm)
  return(mm)
}