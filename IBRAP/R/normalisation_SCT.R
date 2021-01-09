#' Here performs Seurat::SCTransfrom  normalisation
#'
#' SCTransfrom, Seurats novel normalisation technique
#'
#' @import Seurat
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param ... Specify parameters to be passed to the SCTransform function.
#'
#' @examples normalisation_SCT(object = sce_object)
#' @export

normalisation_SCT <- function(object, ...) {
  y <- object
  tmp <- as.Seurat(x = object, counts = 'counts', data = NULL)
  tmp <- SCTransform(object = tmp, do.scale = FALSE, do.center = FALSE,  ...)
  common_rows <- intersect(rownames(tmp), rownames(y))
  y <- y[common_rows,]
  tmp <- tmp[common_rows,]
  assay(y, 'data') <- as.matrix(tmp@assays$SCT@data)
  altExp(object, 'SCTransform') <- y
  return(object)
}

