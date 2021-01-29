#' Read 10x Cell Ranger output files
#'
#' This function allows you concatenate the three matrix files into an expression matrix.
#'
#' @param directory Paste the ABSOLUTE pathway to the file containing the three output files? Files must be named matrix.mtx, genes.tsv, and barcodes.tsv
#' @keywords
#' @export
#' @examples Read10X_output(directory = 'Users/user/10x_output')
#'

Read10X.output <- function(directory) {
  dir.files <- list.files(path = directory)
  if(!'matrix.mtx' %in% dir.files) {
    cat(crayon::cyan('Expected file: matrix.mtx\n'))
  }
  if(!'genes.tsv' %in% dir.files) {
    cat(crayon::cyan('Expected file: genes.tsv\n'))
  }
  if(!'barcodes.tsv' %in% dir.files) {
    cat(crayon::cyan('Expected file: barcodes.tsv\n'))
  }
  count.file <- paste0(directory, '/matrix.mtx')
  genes.file <- paste0(directory, '/genes.tsv')
  barcodes.files <- paste0(directory, '/barcodes.tsv')
  cat(crayon::cyan('Files read\n'))
  genes_ensembl <- read.table(file = genes.file, sep = '\t', header = FALSE)
  barcodes <- read.table(file = barcodes.files, sep = '\t', header = FALSE)
  mm <- Matrix::readMM(count.file)
  if(isUnique(genes_ensembl$V2) == FALSE) {
    cat(crayon::cyan('Non-unique features identified\n'))
    print(duplicated(genes_ensembl$V2))
    genes_ensembl$V2 <- make.unique(genes_ensembl$V2)
  }
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  cat(crayon::cyan('Success: Matrix concatenated\n'))
  mm <- as.matrix(mm)
  return(mm)
}
