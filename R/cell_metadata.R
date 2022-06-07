#' @name cell_metadata
#' @aliases cell_metadata
#' 
#' @title Generates cell-level metadata
#'
#' @description Calculates total genes and total counts detected per cell
#' 
#' @param assay Counts matrix
#' @param col.prefix Which prefix to add to `'_total.counts' and '_total.feature'` as cell metdata column names
#' 
#' @usage cell_metadata(assay = counts, col.prefix = 'RAW')
#'
#' @return Cell level metadata appended to metadata

cell_metadata <- function(assay, 
                          col.prefix) {
  total.counts <- Matrix::colSums(assay)
  total.features <- Matrix::colSums(assay != 0)
  df <- as.data.frame(as.numeric(total.counts))
  df[['total.features']] <- as.numeric(total.features)
  colnames(df) <- c(paste0(col.prefix, '_total.counts'), 
                    paste0(col.prefix, '_total.features'))
  return(df)
}
