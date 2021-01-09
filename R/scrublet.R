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
#'

scrublet <- function(matrix, total_counts = NULL, sim_doublet_ratio = 2.0, n_neighbors = NULL, expected_doublet_rate = 0.075, stdev_doublet_rate = 0.02, random_state = 0L) {
  if(is.matrix(matrix) == FALSE) {
    cat(crayon::cyan('Only an object of class Matrix can be used\n'))
  } else {
    raw_counts <- t(as.data.frame(tmp))
    scrub1 <- scrublet$Scrublet(counts_matrix = raw_counts)
    cat(crayon::cyan('scrublet object created\n'))
    res1 <- scrub1$scrub_doublets(min_counts = 1, min_cells = 1, min_gene_variability_pctl = 85, verbose = TRUE)
    cat(crayon::cyan('doublets detected\n'))
    raw_counts <- t(as.data.frame(raw_counts))
    cat(crayon::cyan('doublet fraction: ', frac, '\n'))
    raw_counts <- raw_counts[,!res1[[2]] == TRUE]
    storage.mode(raw_counts) <- 'integer'
    cat(crayon::cyan('matrix scrubbed\n'))
    return(raw_counts)
    rm(obj, raw_counts, scrub1, res1, scrubbed)
  }

}

mltp <- import('matplotlib.pyplot')

