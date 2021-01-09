#' Benchmark clustering results
#'
#' This function benchmarks the c;ustering methods.
#'
#' @param counts An expression MATRIX containing the raw counts of individual samples.
#' @param norm.method list of normalisation methods used.
#' @param reduction Which reduction method to use.
#' @param n.dim How many dimensions should be used.
#' @param cex.names this value indicates label size
#' @export
#'

plot_reddim_var <- function(object, norm.method, reduction, n.dim, cex.names = 0.6) {
  for(t in norm.method) {
    gg <- apply(X = reducedDim(altExp(object, as.character(t)), as.character(reduction)), MARGIN = 2, FUN = sd)
    barplot(gg, cex.names = cex.names, las = 2, main = paste0(t, '_', reduction, '_var'))
  }
}
