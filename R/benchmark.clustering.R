#' @name benchmark.clustering
#' @aliases benchmark.clustering
#' 
#' @title Benchmarks the cluster assignmentws
#'
#' @description Supervised (ARI and NMI) and unsupervised (ASW, Dunn Index, and Connectivity) benchmarking metrics are calculated for cluster assignments. assays, clustering and reduction for distance calculations are iterated through. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param clustering Character. The names of the cluster assignment dataframes to use
#' @param reduction Character. Which reduction(s) within the assay should be supplied for distance calcultions
#' @param n.dims Numerical. How many dimensions of the reduction should be supplied. Default = 1:3
#' @param dist.method Character. Which distance method should be used, options: 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'. Default = 'euclidean'
#' @param ground.truth Vector. If available, supply the vector in the same order as colnames of the ground truth, i.e. true cell type labels. If this is not supplied, only unsupervised methods will be supplied. Default = NULL
#' 
#' @return Benchmarking scores for the supplied cluster assignments 
#'
#' @export

benchmark.clustering <- function(object, 
                                 assay,
                                 clustering,
                                 reduction, 
                                 n.dims = 1:3,
                                 dist.method='euclidean',
                                 ground.truth=NULL) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string(s) \n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('assay: ', x, 'does not exist\n'))
      
    }
    
  }
  
  if(!is.character(clustering)) {
    
    stop('clustering must be character string(s)\n')
    
  }
  
  for(l in assay) {
    
    for(x in clustering) {
      
      if(!x %in% names(object@methods[[l]]@cluster_assignments)) {
        
        stop('clustering: ,', x, 'not present in assay: ', l, '\n')
        
      }
      
    }
    
    for(x in reduction) {
      
      for(i in assay) {
        
        if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                           object@methods[[i]]@visualisation_reductions, 
                           object@methods[[i]]@integration_reductions))) {
          
          stop(paste0('reduction: ', x, ' does not exist\n'))
          
        }
        
      }
      
    }
    
    if(!is.numeric(n.dims)) {
      
      stop('n.dims must be numerical\n')
      
    }
    
    if(!is.character(dist.method)) {
      
      stop('dist.method must be a character string\n')
      
    }
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[l]]@computational_reductions), 
                   names(object@methods[[l]]@integration_reductions),
                   names(object@methods[[l]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[l]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[l]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[l]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[l]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[l]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[l]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        stop('reductions could not be found\n')
        
      }
      
    }
    
    count <- 1
    
    reduction.list <- reduction.list[reduction]
    reduction.list <- reduction.list[match(reduction, names(reduction.list))]
    
    for(k in clustering) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': benchmarking for assay: ', l, ' cluster dataframe: ', k, '\n')))
      
      reduction_sub <- reduction.list[[reduction[count]]][,n.dims]
      
      count <- count + 1
      
      clusters <- object@methods[[l]]@cluster_assignments[[k]]
      
      for(p in colnames(clusters)) {
        
        if(length(unique(clusters[,p])) <= 1) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': cluster column: ', p, ' contains only 1 cluster group, omitting now\n')))
          clusters[,p] <- NULL
          
        }
        
      }
      
      dist.matrix <- dist(x = reduction_sub, method = dist.method)
      
      sil.results <- data.frame(average_silhoeutte=NA)
      
      for (v in colnames(clusters)[1:length(colnames(clusters))]) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': calculating silhouette for ', v, '\n')))
        tmp <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters[,v])), dist = dist.matrix)
        average <- sum(tmp[,3])/length(tmp[,3])
        sil.results[v,] <- average
        
      }
      
      sil.results <- sil.results[complete.cases(sil.results),]
      max.AS <- max(sil.results)
      print(max.AS)
      
      dunn.results <- data.frame(dunn.index=NA)
      for (p in colnames(clusters)[1:length(colnames(clusters))]){
        
        cat(crayon::cyan(paste0(Sys.time(), ': calculating dunn index for ', p, '\n')))
        dunn.results[p,] <- clValid::dunn(distance = dist.matrix, 
                                          clusters = as.numeric(x = as.factor(x = clusters[,p])))
        
      }
      
      dunn.results <- dunn.results[complete.cases(dunn.results),]
      max.dunn <- max(dunn.results)
      print(max.dunn)
      
      conn.results <- data.frame(connectivity=NA)
      for (p in colnames(clusters)[1:length(colnames(clusters))]){
        
        cat(crayon::cyan(paste0(Sys.time(), ': calculating connectivity for ', p, '\n')))
        conn.results[p,] <- clValid::connectivity(distance = dist.matrix, clusters = clusters[,p])
        
      }
      
      conn.results <- conn.results[complete.cases(conn.results),]
      min.conn <- min(conn.results)
      print(min.conn)
      
      if(!is.null(ground.truth)) {
        
        ARI.results <- data.frame(ARI=NA)
        for (p in colnames(clusters)[1:length(colnames(clusters))]){
          
          cat(crayon::cyan(paste0(Sys.time(), ': calculating ARI for ', p, '\n')))
          ARI.results[p,] <- mclust::adjustedRandIndex(x = clusters[,p], y = ground.truth)
          
        }
        
        ARI.results <- ARI.results[complete.cases(ARI.results),]
        max.ARI <- max(ARI.results)
        print(max.ARI)
        
        NMI.results <- data.frame(NMI=NA)
        for (p in colnames(clusters)[1:length(colnames(clusters))]) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': calculating NMI for ', p, '\n')))
          NMI.results[p,] <- aricode::AMI(c1 = clusters[,p], c2 = ground.truth)
          
        }
        
        NMI.results <- NMI.results[complete.cases(NMI.results),]
        max.nmi <- max(NMI.results)
        print(max.nmi)
        
        results <- cbind(sil.results, dunn.results, conn.results, ARI.results, NMI.results)
        rownames(results) <- colnames(clusters)
        colnames(results) <- c(paste0(k, '_sil.results'), 
                               paste0(k, '_dunn.results'), 
                               paste0(k, '_conn.results'), 
                               paste0(k, '_ARI.results'), 
                               paste0(k, '_NMI.results'))
        
        object@methods[[l]]@benchmark_results[[k]] <- as.data.frame(results)
        
      } else {
        
        results <- cbind(sil.results, dunn.results, conn.results)
        rownames(results) <- colnames(clusters)
        colnames(results) <- c(paste0(k, '_sil.results'), paste0(k, '_dunn.results'), paste0(k, '_conn.results'))
        object@methods[[l]]@benchmark_results[[k]] <- as.data.frame(results)
        
      }
    }
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': completed calcualting benchmarking metrices')))
  
  return(object)
  
}