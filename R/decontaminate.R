#' Dimensionality reduction plot.
#'
#' Here two dimensions will be displayed.
#'
#' @param batch Numeric or character vector. Batch labels for cells. If batch labels are supplied, DecontX is run on cells from each batch separately. Cells run in different channels or assays should be considered different batches. Default NULL.
#' @param maxIter Integer. Maximum iterations of the EM algorithm. Default 500.
#' @param delta Numeric Vector of length 2. Concentration parameters for the Dirichlet prior for the contamination in each cell. The first element is the prior for the native counts while the second element is the prior for the contamination counts. These essentially act as pseudocounts for the native and contamination in each cell. If estimateDelta = TRUE, this is only used to produce a random sample of proportions for an initial value of contamination in each cell. Then fit_dirichlet is used to update delta in each iteration. If estimateDelta = FALSE, then delta is fixed with these values for the entire inference procedure. Fixing delta and setting a high number in the second element will force decontX to be more aggressive and estimate higher levels of contamination at the expense of potentially removing native expression. Default c(10, 10).
#' @param estimateDelta	Boolean. Whether to update delta at each iteration.
#' @param convergence Numeric. The EM algorithm will be stopped if the maximum difference in the contamination estimates between the previous and current iterations is less than this. Default 0.001.
#' @param iterLogLik Integer. Calculate log likelihood every iterLogLik iteration. Default 10.
#' @param varGenes Integer. The number of variable genes to use in dimensionality reduction before clustering. Variability is calcualted using modelGeneVar function from the 'scran' package. Used only when z is not provided. Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter used in 'dbscan' to estimate broad cell clusters. Used only when z is not provided. Default 1.
#' @param seed Integer. Passed to with_seed. For reproducibility, a default value of 12345 is used. If NULL, no calls to with_seed are made.
#' @param logfile Character. Messages will be redirected to a file named 'logfile'. If NULL, messages will be printed to stdout. Default NULL.
#' @param verbose Logical. Whether to print log messages. Default TRUE.
#' @export
#'

decontaminate <- function(matrix,
                          z = NULL,
                          batch = NULL,
                          maxIter = 500,
                          delta = c(10, 10),
                          estimateDelta = TRUE,
                          convergence = 0.001,
                          iterLogLik = 10,
                          varGenes = 5000,
                          dbscanEps = 1,
                          seed = 12345,
                          logfile = NULL,
                          verbose = TRUE) {
  cat(crayon::cyan('Uninformative features omitted\n'))
  d <- decontX(x = matrix,
               z = z,
               batch = batch,
               maxIter = maxIter,
               delta = delta,
               estimateDelta = estimateDelta,
               convergence = convergence,
               iterLogLik = iterLogLik,
               varGenes = varGenes,
               dbscanEps = dbscanEps,
               seed = seed,
               logfile = logfile,
               verbose = verbose)
  cat(crayon::cyan('Decontamination comlpleted\n'))
  clean.matrix <- d$decontXcounts
  cat(crayon::cyan('Matrix isolated\n'))
  mode(clean.matrix) <- 'integer'
  attr(x = clean.matrix) <- plotDecontXContamination(d)
  cat(crayon::cyan('Converting to integer\n'))
  return(clean.matrix)
  cat(crayon::cyan('Finished\n'))
}
