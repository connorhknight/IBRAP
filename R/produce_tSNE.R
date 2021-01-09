#' Uniform Manifold pproximation and Projection (UMAP)
#'
#' performs non-linear dimensionality reduction using UMAP
#'
#' @import ProjectionBasedClustering
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param norm.method Please identify the normalisation techniques.
#' @param init.dims How many dimensions should be used?
#' @param reduction Which reduction method would you like to use?
#' @param ... Further parameters to pass down to the umap function from tSNE
#' @export
#' @examples produce_umap(object = sce_object, multi.norm = TRUE, multi.cluster = TRUE)
#'

produce_tSNE <- function(object, norm.method, init.dims, reduction='pca', ...) {
  for(y in norm.method) {
    print(paste0('Calculating umap for ', y, ' normalisation method.'))
    temp <- altExp(object, y)
    c <- tSNE(DataOrDistances = reducedDim(temp, reduction)[,init.dims],
              Iterations = 1000, ...)$ProjectedPoints
    colnames(c) <- c('tsne_1', 'tsne_2')
    reducedDim(temp, 'tsne') <- c
    altExp(object, y) <- temp
  }
  return(object)
}
