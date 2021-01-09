#' An integrated clustering step
#'
#' Includes: normalisation, scale/centre, unwanted variance regression, linear dimensionality reduction (PCA), and benchmarking.
#' This function iterates through all normalisation techniques and benchmarks to indicate the superior method.
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param norm.method Specify which normalisation method you would like to process
#' @param cluster.method Identify the cluster methods used
#' @param num.core How many cores should be used during SC3 processesing.
#' @param reduction Which reduction should be the input.
#' @param dims.use How many dimensions of the reduction should be used?
#' @param nl.reduction.method Which non-linear reduction method should be applied? 'umap', 'tsne'
#' @param benchmark TRUE/FALSE whether benchmarking should be performed.
#' @param dist.method Which distance method should be used to calculate distance. Available: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'
#' @export
#'

sc_Cluster <- function(object, norm.method, cluster.method, res=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5),
                       prune.SNN=0, nn.method='annoy', annoy.metric='euclidean', nn.eps=0.0, num.core, ks, reduction, dims.use, nl.reduction.method,
                       benchmark, dist.method='euclidean', ARI, ground.truth) {
  if('sc3' %in% cluster.method) {
    object <- SC3_cluster(object = object, norm.methods = norm.method, n.core = num.core, ks = ks)
  }
  if('seurat' %in% cluster.method) {
    object <- Seurat_cluster(object = object, norm.method = norm.method, reduction = reduction, dims = dims.use, res = res,
                             prune.SNN = prune.SNN, nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
  }
  if('umap' %in% nl.reduction.method) {
    object <- produce_umap(object = object, norm.methods = norm.method, reduction = reduction)
  }
  if('tsne' %in% nl.reduction.method) {
    object <- produce_tSNE(object = object, norm.method = norm.method, init.dims = dims.use, reduction = reduction)
  }
  if(benchmark == TRUE) {
    object <- benchmarking_clustering(object = object, components = dims.use, reduction = 'pca', dist.method = dist.method,
                                      norm.methods = norm.method, cluster.methods = cluster.method, ARI = ARI, ground.truth = ground.truth)
  }
  return(object)
}
