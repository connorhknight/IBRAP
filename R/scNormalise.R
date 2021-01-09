#' An integrated normalisation step
#'
#' Includes: normalisation, scale/centre, unwanted variance regression, linear dimensionality reduction (PCA), and benchmarking.
#' This function iterates through all normalisation techniques and benchmarks to indicate the superior method.
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param method Specify 'all' to perform all methods and store them as alternative experiemnts
#' @param scale Should the data be scaled?
#' @param centre Should the data be centre transformed?
#' @param unwanted.variance Specify which metadata columns should be regressed.
#' @param n.features How many highly variable genes should be retained.
#' @param benchmark TRUE/FALSE whether benchmarking should be performed.
#' @param distance.method SpecWhich distance method should be used to calculate distance. Available: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'
#' @param ... If only a single method is applied, you may pass parameters down to the selected function here.
#' @export
#'

scNormalise <- function(object, method = 'scran', scale, centre, unwanted.variance, n.features=2000, benchmark, distance.method,
                        scran.max.cluster.size = 1000, scran.scaling=NULL,
                        scran.log=TRUE, scran.center_size_factors=TRUE) {
  if('scran' %in% method){
    object <- normalisation_scran(object, max.cluster.size = scran.max.cluster.size, scaling=scran.scaling,
                                  log=scran.log, center_size_factors=scran.center_size_factors)
    print('Norm method: Scran, completed.')
  }
  if('logtransform' %in% method){
    object <- normalisation_SeuratLog(object)
    print('Norm method: Seurat LogNorm, completed.')
  }
  if('SCTransform' %in% method){
    object <- normalisation_SCT(object)
    print('Norm method: SCTransform, completed.')
  }
  if('cpm' %in% method){
    object <- normalisation_cpm(object)
    print('Norm method: CPM, completed.')
  }
  object <- find_HVGs(object = object, nfeatures = n.features)
  object <- Scale_data(object = object, unwanted.variance = unwanted.variance, scale = scale, centre = centre)

  ### PCA Analysis ###
  names.of.experiments <- altExpNames(object)
  for(u in names.of.experiments) {
    print(paste0('Calculating PCs for ', u, ' normalisation method.'))
    temp.2 <- altExp(object, u)
    temp.2 <- performPCA(object = temp.2)
    altExp(object, u) <- temp.2
    print(paste0('Finished calculating PCs for ', u, ' normalisation method.'))
  }

  ### benchmarking ###
  if(benchmark==TRUE) {
    names.of.experiments <- altExpNames(object)
    results <- list()
    for(u in names.of.experiments) {
      print(paste0('Calculating benchmarking metrices for ', u, ' normalisation method.'))
      temp.3 <- altExp(object, u)
      bench.results <- benchmark_norm(object = temp.3, DR.assay = 'pca', distance.method = distance.method)
      altExp(object, u) <- bench.results
      print(paste0('Finished calculating benchmarking metrices for ', u, ' normalisation method.'))
      }
  }

  return(object)
}
