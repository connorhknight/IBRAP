#' Produces metadata: total features and total counts per cell for different assays
#'
#' @param object A SingleCellExperiment
#' @param assay Provide the name of the assay for this metadata to be generated
#' @param prefix Which prefix should be used for the column names
#' @keywords
#' @export
#' @examples Read10X_output(directory = 'Users/user/10x_output')
#'

metadata.generator <- function(object, 
                               assay, 
                               prefix) {
  cat(crayon::cyan('Annotated with fundamental metadata for assay: ', 'assay', '\n'))
  colData(object)[[paste0(prefix, '_total.counts')]] <- colSums(assay(object, assay))
  colData(object)[[paste0(prefix, '_total.features')]] <- colSums(assay(object, assay) != 0)
  cat(crayon::cyan('Completed\n'))
}