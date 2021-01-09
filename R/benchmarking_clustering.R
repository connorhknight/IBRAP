#' Benchmark clustering results
#'
#' This function benchmarks the c;ustering methods.
#'
#' @import cluster
#' @import clValid
#' @import mclust
#' @param counts An expression MATRIX containing the raw counts of individual samples.
#' @param components How many components to use from dimensionality reduction method
#' @param reduction Which reduction method to use.
#' @param dist.method Which distance method should be used to calculate distance. Available: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'
#' @param multi.norm TRUE/FALSE, whether multiple normalisation methods were used.
#' @param multi.cluster TRUE/FALSE, whether multiple clustering alogirthms were used.
#' @export
#' @examples benchmarking_clustering(counts = counts.matrix, components = 1:19, reduction = 'pca', dist.method = 'euclidean', multi.norm = TRUE, multi.cluster = TRUE)
#'

benchmarking_clustering <- function(object, components, reduction, dist.method='euclidean',
                                    norm.methods, cluster.methods, ARI=FALSE, ground.truth) {
  for(c in norm.methods) {
    temp <- altExp(object, c)
    print(temp)

    for(g in cluster.methods) {
      all.clusters <- as.data.frame(metadata(temp)[['clustering']][[g]])
      print(all.clusters)
      dims <- reducedDim(temp, reduction)[,components]
      dist.matrix <- dist(x = dims, method = dist.method)
      all.clusters <- all.clusters[, colSums(all.clusters != 0) > 0]
      all.clusters <- all.clusters[!is.na(names(all.clusters))]
      print(all.clusters)
      sil.results <- data.frame(average_silhoeutte=NA)
      for (v in colnames(all.clusters)[1:length(colnames(all.clusters))]) {
        print(paste0('Calculating silhouette for ', v))
        tmp <- silhouette(x = as.numeric(x = as.factor(x = all.clusters[,v])), dist = dist.matrix)
        average <- sum(tmp[,3])/length(tmp[,3])
        sil.results[v,] <- average
      }
      sil.results <- sil.results[complete.cases(sil.results),]
      max.AS <- max(sil.results)
      print(max.AS)

      dunn.results <- data.frame(dunn.index=NA)
      for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
        print(paste0('Calculating dunn index for ', p))
        dunn.results[p,] <- dunn(distance = dist.matrix, clusters = as.numeric(x = as.factor(x = all.clusters[,p])))
      }
      dunn.results <- dunn.results[complete.cases(dunn.results),]
      max.dunn <- max(dunn.results)

      conn.results <- data.frame(connectivity=NA)
      for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
        print(paste0('Calculating connectivity for ', p))
        conn.results[p,] <- connectivity(distance = dist.matrix, clusters = all.clusters[,p])
      }
      conn.results <- conn.results[complete.cases(conn.results),]
      max.conn <- max(conn.results)

      if(ARI == TRUE) {
        ARI.results <- data.frame(ARI=NA)
        for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
          print(paste0('Calculating ARI for ', p))
          ARI.results[p,] <- adjustedRandIndex(x = all.clusters[,p], y = ground.truth)
        }
        ARI.results <- ARI.results[complete.cases(ARI.results),]
        max.ARI <- max(ARI.results)
        results <- cbind(sil.results, dunn.results, conn.results, ARI.results)
        rownames(results) <- colnames(all.clusters)
        colnames(results) <- c(paste0(g, '_sil.results'), paste0(g, '_dunn.results'), paste0(g, '_conn.results'), paste0(g, '_ARI.results'))
        metadata(temp)[['benchmarking_clustering']][[as.character(g)]] <- results
      } else {
        results <- cbind(sil.results, dunn.results, conn.results)
        rownames(results) <- colnames(all.clusters)
        colnames(results) <- c(paste0(g, '_sil.results'), paste0(g, '_dunn.results'), paste0(g, '_conn.results'), paste0(g, '_ARI.results'))
        metadata(temp)[['benchmarking_clustering']][[as.character(g)]] <- results
      }
    }
    altExp(object, c) <- temp
  }
  return(object)
}
