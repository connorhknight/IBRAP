#' Principal Component Analysis (PCA)
#'
#' @import PCAtools
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param ... Specify parameters to be passed to the pca function from pcaMethods.
#'
#' @examples performPCA(object = sce_object, unwanted.variance = 'percent.mt', scale = TRUE, centre = TRUE)
#' @export

performPCA <- function(object, norm.methods, ...) {

  for(t in norm.methods) {
    t <- altExp(object, as.character(t))
    a <- pca(mat = assay(t, 'scaled'), center = FALSE, scale = FALSE, ...)
    reducedDim(object, 'pca') <- as.matrix(a$rotated[,1:50])
    return(object)
  }

}
