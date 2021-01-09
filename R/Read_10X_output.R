#' Read 10x Cell Ranger output files
#'
#' This function allows you concatenate the three matrix files into an expression matrix.
#'
#' @import utils
#' @import Matrix
#' @param directory Paste the ABSOLUTE pathway to the file containing the three output files? Files must be named matrix.mtx, genes.tsv, and barcodes.tsv
#' @keywords 10x CellRanger
#' @export
#' @examples Read10X_output('place directory here')
#'

Read10X_output <- function(directory) {
  count.file <- paste0(directory, '/matrix.mtx')
  genes.file <- paste0(directory, '/genes.tsv')
  barcodes.files <- paste0(directory, '/barcodes.tsv')
  genes_ensembl <- read.table(file = genes.file, sep = '\t', header = FALSE)
  barcodes <- read.table(file = barcodes.files, sep = '\t', header = FALSE)
  mm <- readMM(count.file)
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  mm <- as.matrix(mm)
  return(mm)
}
