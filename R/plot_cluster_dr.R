#' Dimensionality reduction plot.
#'
#' Here two dimensions will be displayed.
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param reduction Which dimensionality reduction assay should be utilised?
#' @param pt.size What size should be the data points be.
#' @param multi.norm TRUE/FALSE were multiple normalisation techniques performed?
#' @param multi.cluster TRUE/FALSE were multiple clustering techniques performed?
#' @param norm.method Which normalisation method wouyld you like to access?
#' @param cluster.method Which clustering method would you like to access?
#' @param cluster.column Which column within these results would you like to group by?
#' @export
#' @examples plot_cluster_dr(object = sce, reduction='umap', pt.size=0.5, multi.norm=TRUE, multi.cluster=TRUE, norm.method='scran', cluster.method='sc3', cluster.column='SCT_snn_res.0.4')
#'

plot_cluster_dr <- function(object, reduction='umap', pt.size=0.5, norm.method, cluster.method, cluster.column, interactive=FALSE) {
  print('plot_cluster_dr_started')
  iso <- metadata(altExp(object, as.character(norm.method)))[['clustering']][[as.character(cluster.method)]][[as.character(cluster.column)]]
  print('.')
  reduction <- as.data.frame(reducedDim(altExp(object, as.character(norm.method)), as.character(reduction)))
  print('.')
  meta <- colData(altExp(object, as.character(norm.method)))
  print('.')
  meta.colnames <- colnames(meta)
  print('.')
  barcodes <- colnames(altExp(object, as.character(norm.method)))
  print('.')
  x.val <- reduction[,1]
  print('.')
  y.val <- reduction[,2]
  print('.')
  results <- cbind(barcodes, meta, iso, x.val, y.val)
  print('.')
  colnames(results) <- c('barcode', meta.colnames, 'clusters', 'x_value', 'y_value')
  print('.')
  results <- as.data.frame(results)
  print(results)
  plot <- ggplot(data = results, aes(ident = `barcode`, x = `x_value`,
                                 y = `y_value`, color = `clusters`)) +
    geom_point(shape = 16, size = pt.size) +
    scale_color_hue() +
    guides(fill=FALSE, alpha=FALSE, size=FALSE) +
    theme_minimal() +
    xlab(paste0(reduction, '_1')) + ylab(paste0(reduction, '_2')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5))
  if(interactive == TRUE) {
    z <- ggplotly(plot)
    return(z)
  } else {
    return(plot)
  }
  print('plot_cluster_dr_completed')
}
