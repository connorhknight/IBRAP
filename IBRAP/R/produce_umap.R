#' Uniform Manifold pproximation and Projection (UMAP)
#'
#' performs non-linear dimensionality reduction using UMAP
#'
#' @import uwot
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param multi.norm TRUE/FALSE, were multiple normalisation techniques used?
#' @param multi.cluster TRUE/FALSE, were multiple clustering techniques used?
#' @param ... Further parameters to pass down to the umap function from uwot.
#' @export
#' @examples produce_umap(object = sce_object, multi.norm = TRUE, multi.cluster = TRUE)
#'

produce_umap <- function(object, norm.methods, reduction='pca', n.dim, ...) {
  for(p in norm.methods) {
    print(paste0('Calculating umap for ', p, ' normalisation method.'))
    temp <- altExp(object, p)
    c <- umap(X = reducedDim(temp, reduction)[,n.dim], verbose = TRUE, ...)
    colnames(c) <- c('umap_1', 'umap_2')
    reducedDim(temp, 'umap') <- c
    altExp(object, p) <- temp
  }
  return(object)
  }

