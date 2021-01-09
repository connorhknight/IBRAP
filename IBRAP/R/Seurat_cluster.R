#' Seurat clustering
#'
#' This function utilises Seurats graph-based clustering alogrithm
#'
#' @import Seurat
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param multi.method TRUE/FALSE if more than one normalisation method is being used.
#' @param reduction Which dimensionality reduction assay would you like to use?
#' @param res This specifies the resolution to be passed down to seurat.
#' @param dims how many dimensionality reduction components should be used.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the strigency of pruning (0 — no pruning, 1 — prune everything).
#' @param nn.method Method for nearest neighbor finding. Options include: rann, annoy
#' @param annoy.metric Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
#' @param nn.eps 	Error bound when performing nearest neighbor seach using RANN; default of 0.0 implies exact nearest neighbor search
#' @param ... Parameters to pass down to the FindClusters function.
#'
#' @examples Seurat_cluster(object = sce_object, multi.method, reduction='pca', res=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5), dims=1:12, prune.SNN=0, nn.method='annoy', annoy.metric='euclidean', nn.eps=0.0)
#' @export

Seurat_cluster <- function(object, norm.method, reduction='pca', res=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5), dims=1:12,
                           prune.SNN=0, nn.method='annoy', annoy.metric='euclidean', nn.eps=0.0, ...) {
  for(o in norm.method) {
    z <- list()
    print(paste0('Calculating Seurat clusters for ', o, ' normalisation method.'))
    temp <- altExp(object, o)
    orig.names <- colnames(colData(object))
    tmp <- as.Seurat(x = temp, counts='counts', data='data')
    tmp <- FindNeighbors(object = tmp, reduction = reduction, verbose = TRUE, dims = dims, compute.SNN = TRUE, prune.SNN = prune.SNN,
                         nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
    tmp <- FindClusters(object = tmp, resolution = res, ...)
    new.names <- colnames(tmp@meta.data)
    sep.names <- new.names[!(new.names %in% orig.names)]
    sep.names <-sep.names[1:length(sep.names)-1]
    new.clusters <- tmp@meta.data[,sep.names]
    for(t in res) {
      z[length(z)+1] <- paste0('Seurat_res_', t)
    }
    colnames(new.clusters) <- unlist(z)
    metadata(temp)$clustering$seurat <- new.clusters
    print(metadata(temp))
    altExp(object, o) <- temp
  }
  return(object)
}
