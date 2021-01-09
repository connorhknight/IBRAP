#' Benchmarking of normalisation techniques
#'
#' @import scran
#' @import edgeR
#' @import Seurat
#' @import PCAtools
#' @import cluster
#' @import clValid
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param DR.assay Which dimensionality reduction assay should be utilised.
#' @param distance.method Which distance method should be used to calculate distance. Available: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'
#' @param n.k Either a single value, a list or a range of k to calculate.
#' @export
#' @examples benchmark_norm(object = sce_object, DR.assay = 'pca', distance.method = 'euclidean', n.k=2:25)
#'

benchmark_norm <- function(object, DR.assay, distance.method, n.k=2:25) {
  PAM_clusters <- data.frame(barcodes=colnames(object))
  DR <- reducedDim(object, as.character(DR.assay))[,1:3]
  dist.matrix <- dist(x = DR, method = distance.method)

  for(v in n.k){
    print(paste0('Clustering, k=', v))
    tmp <- as.data.frame(pam(x = dist.matrix, k = v, diss = TRUE, metric = distance.method, cluster.only = TRUE, do.swap = FALSE))
    tmp <- cbind(rownames(tmp), tmp)
    PAM_clusters <- cbind(PAM_clusters, tmp[,2])
    new.names <- c(colnames(PAM_clusters)[1:length(colnames(PAM_clusters))-1], paste0(v))
    colnames(PAM_clusters) <- new.names
  }
  rownames(PAM_clusters) <- PAM_clusters$barcodes
  PAM_clusters <- PAM_clusters[2:length(colnames(PAM_clusters))]

  sil.results <- data.frame(average_silhoeutte=NA)
  for (v in colnames(PAM_clusters)[1:length(colnames(PAM_clusters))]) {
    print(paste0('Calculating silhouette for k=', v))
    tmp <- silhouette(x = as.numeric(x = as.factor(x = PAM_clusters[,v])), dist = dist.matrix)
    average <- sum(tmp[,3])/length(tmp[,3])
    sil.results[v,] <- average
  }
  sil.results <- sil.results[complete.cases(sil.results),]
  max.AS <- max(sil.results)

  dunn.results <- data.frame(dunn.index=NA)
  for (p in colnames(PAM_clusters)[1:length(colnames(PAM_clusters))]){
    print(paste0('Calculating dunn index for k=', p))
    dunn.results[p,] <- dunn(distance = dist.matrix, clusters = PAM_clusters[,p])
  }
  dunn.results <- dunn.results[complete.cases(dunn.results),]
  max.dunn <- max(dunn.results)

  print(paste0('Maximum average silhouette width: ', max.AS))
  print(paste0('Maximum dunn index: ', max.dunn))
  df <- data.frame(metrics=c(max.AS, max.dunn))
  rownames(df) <- c('max.AS', 'max.dunn')
  results <- cbind(sil.results, dunn.results)
  metadata(object)[['benchmarking_normalisation']][['PAM_clustering']] <- PAM_clusters
  metadata(object)[['benchmarking_normalisation']][['benchmark_results']] <- results
  return(object)
}
