#' @name feature_metadata
#' @aliases feature_metadata
#' 
#' @title Generates feature level metadata
#'
#' @description Calculates total counts and cells detected per gene
#' 
#' @import Matrix
#' 
#' @param assay Counts matrix
#' @param col.prefix Which prefix to add to `'_total.counts' and '_total.cells'` as feature metdata column names
#' 
#' @usage feature_metadata(assay = counts, col.prefix = 'RAW')
#'
#' @return Gene level metadata appended to metadata
#'
#' @export

feature_metadata <- function(assay, 
                             col.prefix) {
  df <- as.data.frame(as.numeric(rowSums(assay)))
  rownames(df) <- rownames(assay)
  df$temp <- as.numeric(rowSums(assay > 0))
  colnames(df) <- c(paste0(col.prefix,'_total.counts'), paste0(col.prefix,'_total.cells'))
  return(df)
}