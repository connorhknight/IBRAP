#' Here performs scran package normalisation
#'
#' Scran normalisation comprises of pre-calculating clustering, computing the sum factors and then log normalising with the factors.
#' Scrans cutting edge is that it functions well wit high zero count data matrices.
#'
#' @import scran
#' @import scater
#' @param object Please specify a SCE object produced using IBRAP functions.
#'
#' @examples normalisation_scran(object = sce_object)
#' @export

normalisation_scran <- function(object, max.cluster.size = 1000, scaling=NULL,
                                log=TRUE, center_size_factors=TRUE) {
  print('initialisaing quickCluster')
  clusters <- quickCluster(object)
  print('quickCluster completed')
  print('initialisaing computeSumFactors')
  sce.scran <- computeSumFactors(object, clusters=clusters, max.cluster.size=max.cluster.size, scaling=scaling)
  print('computeSumFactors completed')
  print('initialisaing LogNormCounts')
  log <- logNormCounts(x = sce.scran, log=log, center_size_factors=center_size_factors, exprs_values='counts')
  print('LogNormCounts completed')
  y <- object
  assay(y, 'data') <- assay(log, 'logcounts')
  altExp(object, 'scran') <- y
  return(object)
}
