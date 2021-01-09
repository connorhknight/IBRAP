#' Benchmark metric visualisation tool for clustering
#'
#' Here will display average silhouette width, dunn index, cluster connectivity, and adjusted rand index.
#'
#' @import egg
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param norm.method Identify which normalisation method you would like to access.
#' @param clust.method Which cliuster method would you like to assess?
#' @export
#'

plot_benchmarking <- function(object, clust.method, ARI){
  print('.')
  clust.bench <- metadata(object)[['benchmarking_clustering']][[as.character(clust.method)]]
  print('.')
  clust.bench <- as.data.frame(clust.bench)
  print(clust.bench)
  clust.bench[,'cluster_index'] <- rownames(clust.bench)
  if(ARI == TRUE) {
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'ARI', 'cluster_index')
  } else {
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'cluster_index')
  }

  colnames(clust.bench) <- labels
  print(clust.bench)
  list.plot <- list()

  for(o in 1:sum(length(labels)-2)) {

    label <- labels[as.numeric(o)]
    print(label)
    fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(label), group = 1)) +
      geom_point() +
      geom_line() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    list.plot[[as.numeric(o)]] <- fig

  }

  last.label <- labels[as.numeric(sum(length(labels)-1))]

  last.fig <- fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(last.label), group = 1)) +
    geom_point() +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  list.plot[[as.numeric(sum(length(labels)-1))]] <- last.fig
  print('/')
  do.call('ggarrange', c(plots = list.plot, ncol = 1))
}
