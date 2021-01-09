#' Benchmarking of normalisation techniques
#'
#' This function activates the shiny app for investigating the c;ustering results
#'
#' @import shiny
#' @export
#' @examples benchmark_norm(object = sce_object, DR.assay = 'pca', distance.method = 'euclidean', n.k=2:25)
#'

Run_IBRAP_clustering_browser <- function() {
  appDir <- system.file('shiny-examples', 'myapp', package = 'IBRAP')
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
