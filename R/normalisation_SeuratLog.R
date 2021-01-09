#' Here performs Seurat log transformation
#'
#' This method logtransforms the data.
#'
#' @import Seurat
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param ... Specify parameters to be passed to the NormalizeData function.
#'
#' @examples normalisation_SeuratLog(object = sce_object)
#' @export


normalisation_SeuratLog <- function(object, ...) {
  y <- object
  tmp <- as.Seurat(x = object, counts = 'counts', data = NULL)
  results <- NormalizeData(object = tmp, ...)
  assay(y, 'data') <- as.matrix(results@assays$RNA@data)
  altExp(object, 'logtransform') <- y
  return(object)
}
