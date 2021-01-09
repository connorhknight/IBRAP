#' SC3 consensus clustering
#'
#' Functions by using hierachial clustering to determine a consensus across kmeans calculations
#'
#' @import SC3
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param multi.method TRUE/FALSE, were multiple normnalisation techniques used?
#' @param ks Specify k value (number of total clusters)
#' @param n.core SC3 is computationally demanding, and thus may require multiple cores. please specify here. For a local computer use (total cores - 1), in a hpc use all cores
#' @export
#' @examples SC3_cluster(object = sce_object, multi.method = TRUE, ks = 12, n.core = 2)
#'

SC3_cluster <- function(object, norm.methods, ks, n.core=2) {
  for(t in norm.methods){
    print(paste0('Calculating SC3 clusters for ', t, ' normalisation method.'))
    temp <- altExp(object, t)
    temp.2 <- temp
    logcounts(temp.2) <- assay(temp.2, 'data')
    rowData(temp.2)$feature_symbol <- rownames(temp.2)
    temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
    temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
    temp.2 <- sc3_calc_dists(temp.2)
    temp.2 <- sc3_calc_transfs(temp.2)
    temp.2 <- sc3_kmeans(temp.2, ks = ks)
    temp.2 <- sc3_calc_consens(temp.2)
    orig.names <- colnames(colData(temp))
    new.names <- colnames(colData(temp.2))
    sep.names <- new.names[!(new.names %in% orig.names)]
    new.clusters <- colData(temp.2)[,sep.names]
    metadata(temp)$clustering$sc3 <- as.data.frame(new.clusters)
    altExp(object, t) <- temp
    print('SC3 multi-normalisation clustering complete')
  }
  return(object)
}
