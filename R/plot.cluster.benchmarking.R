#' @name plot.cluster.benchmarking
#' @aliases plot.cluster.benchmarking
#' 
#' @param object An IBRAP S4 class object
#' @param assay Character. Which assay within the object to access
#' @param clustering Character. Which clustering benchmarking results to utilise
#' @param ARI Boolean. Were the ARI-dependent metrices (ARI and NMI) calculated. 
#' 
#' @examples 
#' 
#' # To plot if you supplied ground truth labels during becnhmarking
#' plot.cluster.benchmarking(object = object, assay = 'SCT', clustering = 'pca_nn.v1:louvain', ARI = T)
#' 
#' # To plot if you didn't supplied ground truth labels during becnhmarking
#' plot.cluster.benchmarking(object = object, assay = 'SCT', clustering = 'pca_nn.v1:louvain', ARI = F)
#' 
#' @export plot.cluster.benchmarking

plot.cluster.benchmarking <- function(object, 
                                      assay,
                                      clustering, 
                                      ARI = F){
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string \n')
    
  } else if (is.character(assay)) {
    
    if(!assay %in% names(object@methods)) {
      
      stop('assay is not contained within object@methods \n')
      
    }
    
  }
  
  if(!is.character(clustering)) {
    
    stop('clustering must be a character string \n')
    
  } else {
    
    if(!clustering %in% names(object@methods[[assay]]@benchmark_results$clustering)) {
      
      stop('clustering not contained within benchmark_results \n')
      
    }
    
  }
  
  if(!is.logical(ARI)){
    
    stop('ARI must be boolean. TRUE/FALSE \n')
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  clust.bench <- object@methods[[assay]]@benchmark_results$clustering[[clustering]]
  clust.bench <- as.data.frame(clust.bench)
  
  clust.bench[,'cluster_index'] <- rownames(clust.bench)
  
  if(ARI == TRUE) {
    
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'ARI', 'NMI', 'cluster_index')
    
  } else {
    
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'cluster_index')
    
  }
  
  colnames(clust.bench) <- labels
  
  list.plot <- list()
  
  for(o in 1:sum(length(labels)-2)) {
    
    label <- labels[as.numeric(o)]
    
    fig <- ggplot2::ggplot(clust.bench, ggplot2::aes_string(x = 'cluster_index', y = as.character(label), group = 1)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() + 
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
    
    list.plot[[as.numeric(o)]] <- fig
    
  }
  
  last.label <- labels[as.numeric(sum(length(labels)-1))]
  
  last.fig <- ggplot2::ggplot(clust.bench, ggplot2::aes_string(x = 'cluster_index', y = as.character(last.label), group = 1)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
  
  list.plot[[as.numeric(sum(length(labels)-1))]] <- last.fig
  
  if(ARI == TRUE) {
    
    do.call('ggarrange.tmp', c(plots = list.plot, ncol = 5))
    
  } else {
    
    do.call('ggarrange.tmp', c(plots = list.plot, ncol = 3))
    
  }
  
}
