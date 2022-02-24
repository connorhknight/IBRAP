#' @name perform.decontX
#' @aliases perform.decontX
#' 
#' @title R package: celda, decontX function
#'
#' @description Removes ambient RNA from datasets
#' 
#' @param counts Counts matrix
#' @param z Cluster assignments for cells
#' @param maxIter Maximum number of iterations to be performed
#' @param delta Numeric Vector of length 2. Concentration parameters for the Dirichlet prior for the contamination in each cell. The first element is the prior for the native counts while the second element is the prior for the contamination counts. These essentially act as pseudocounts for the native and contamination in each cell. If estimateDelta = TRUE, this is only used to produce a random sample of proportions for an initial value of contamination in each cell. Then fit_dirichlet is used to update delta in each iteration. If estimateDelta = FALSE, then delta is fixed with these values for the entire inference procedure. Fixing delta and setting a high number in the second element will force decontX to be more aggressive and estimate higher levels of contamination at the expense of potentially removing native expression. Default c(10, 10).
#' @param estimateDelta Boolean. Whether to update delta at each iteration.
#' @param iterLogLik Integer. Calculate log likelihood every iterLogLik iteration. Default 10.
#' @param varGenes Integer. The number of variable genes to use in dimensionality reduction before clustering. Variability is calcualted using modelGeneVar function from the 'scran' package. Used only when z is not provided. Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter used in 'dbscan' to estimate broad cell clusters. Used only when z is not provided. Default 1.
#' @param verbose Logical. Should function information be printed to the console? Default = FALSE
#' @param print.plot Logical. Should the UMAP plot displaying contaimination in each cell be printed? Default = FALSE
#' @param seed Integer. Passed to with_seed. For reproducibility. Default = 1234
#' 
#' @usage perform.decontX(counts = counts)
#' 
#' @return Doublet-omitted sparse matrix
#' 
#' @examples counts <- perform.decontx(counts = counts)
#'
#' @export

perform.decontX <- function(counts,
                            z = NULL,
                            maxIter = 500,
                            delta = c(10, 10),
                            estimateDelta = TRUE,
                            convergence = 0.001,
                            iterLogLik = 10,
                            varGenes = 5000,
                            dbscanEps = 1,
                            print.plot=F,
                            verbose=F,
                            seed = 1234) {
  
  if(!is(object = counts, class2 = 'matrix')) {
    
    if (!is(object = counts, class2 = 'dgCMatrix')) {
      
      stop('counts must be in matrix or dgCMatrix format\n')
      
    }
    
  } 
  
  if(!is.null(z)) {
    
    if(length(z) != ncol(counts)) {
      
      stop('z must have the same length as ncol: counts\n')
      
    }
    
  }
  
  if(!is.numeric(maxIter)) {
    
    stop('maxIter must be numerical\n')
    
  }
  
  if(!is.numeric(delta)) {
    
    stop('delta must be numerical\n')
    
  }
  
  if(!is.logical(estimateDelta)) {
    
    stop('estimateDelta must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.numeric(convergence)) {
    
    stop('convergence must be numerical\n')
    
  }
  
  if(!is.numeric(iterLogLik)) {
    
    stop('iterLogLik must be numerical\n')
    
  }
  
  if(!is.numeric(varGenes)) {
    
    stop('varGenes must be numerical\n')
    
  }
  
  if(!is.numeric(dbscanEps)) {
    
    stop('dbscanEps must be numerical\n')
    
  }
  
  if(!is.logical(print.plot)) {
    
    stop('print.plot should be logical. TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed must be numerical\n')
    
  }
  
  if(isTRUE(verbose)) {
    
    d <- celda::decontX(x = counts,
                        z = z,
                        maxIter = maxIter,
                        delta = delta,
                        estimateDelta = estimateDelta,
                        convergence = convergence,
                        iterLogLik = iterLogLik,
                        varGenes = varGenes,
                        dbscanEps = dbscanEps,
                        seed = seed,
                        verbose = TRUE)
    
  } else if (isFALSE(verbose)) {
    
    d <- celda::decontX(x = counts,
                        z = z,
                        maxIter = maxIter,
                        delta = delta,
                        estimateDelta = estimateDelta,
                        convergence = convergence,
                        iterLogLik = iterLogLik,
                        varGenes = varGenes,
                        dbscanEps = dbscanEps,
                        seed = seed,
                        verbose = FALSE)
    
  }
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': decontamination completed\n')))
    
  }
 
  if(isTRUE(print.plot)) {
    
    print(celda::plotDecontXContamination(x = d))
    
  }
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': ', as.character(formatC(sum(sum(d$contamination)/length(d$contamination)), digits = 2)), '% average contamination\n')))
    
  }

  clean.matrix <- d$decontXcounts
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': matrix isolated\n')))
    
  }

  clean.matrix <- round(clean.matrix)
  
  zero.samples <- Matrix::colSums(as.matrix(clean.matrix)) > 0
  
  clean.matrix <- clean.matrix[,zero.samples]
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': converted to integer\n')))
    
  }

  return(clean.matrix)
  
}