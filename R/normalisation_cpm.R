#' Here performs cpm package normalisation
#'
#'
#' @import edgeR
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param ... Here are functions to pass dowen to edgeR::cpm()
#'
#' @examples normalisation_cpm(object = sce_object)
#' @export

normalisation_cpm <- function(object, ...) {
  y <- object
  assay(y, 'data') <- cpm(assay(object, 'counts'), log = TRUE, ...)
  altExp(object, 'cpm') <- y
  return(object)
}
