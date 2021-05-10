#' @name perform.decontX
#' @aliases perform.decontX
#' 
#' @title R package: celda, decontX function
#'
#' @description Removes ambient RNA from datasets
#' 
#' @import crayon
#' @import celda
#' 
#' @param counts Counts matrix
#' @param z Cluster assignments for cells
#' @param batch Batches for each cell if multiple samples are present
#' @param maxIter Maximum number of iterations to be performed
#' @param delta Numeric Vector of length 2. Concentration parameters for the Dirichlet prior for the contamination in each cell. The first element is the prior for the native counts while the second element is the prior for the contamination counts. These essentially act as pseudocounts for the native and contamination in each cell. If estimateDelta = TRUE, this is only used to produce a random sample of proportions for an initial value of contamination in each cell. Then fit_dirichlet is used to update delta in each iteration. If estimateDelta = FALSE, then delta is fixed with these values for the entire inference procedure. Fixing delta and setting a high number in the second element will force decontX to be more aggressive and estimate higher levels of contamination at the expense of potentially removing native expression. Default c(10, 10).
#' @param estimateDelta Boolean. Whether to update delta at each iteration.
#' @param iterLogLik Integer. Calculate log likelihood every iterLogLik iteration. Default 10.
#' @param varGenes Integer. The number of variable genes to use in dimensionality reduction before clustering. Variability is calcualted using modelGeneVar function from the 'scran' package. Used only when z is not provided. Default 5000.
#' @param dbscanEps Numeric. The clustering resolution parameter used in 'dbscan' to estimate broad cell clusters. Used only when z is not provided. Default 1.
#' @param seed Integer. Passed to with_seed. For reproducibility, a default value of 12345 is used. If NULL, no calls to with_seed are made.
#' 
#' @usage perform.decontX(counts = counts)
#' 
#' @return Doublet-omitted sparse matrix
#'
#' @export

perform.decontX <- function(counts,
                            z = NULL,
                            batch = NULL,
                            maxIter = 500,
                            delta = c(10, 10),
                            estimateDelta = TRUE,
                            convergence = 0.001,
                            iterLogLik = 10,
                            varGenes = 5000,
                            dbscanEps = 1,
                            seed = 12345) {
  
  if(!is(object = counts, class2 = 'matrix')) {
    
    if (!is(object = counts, class2 = 'dgCMatrix')) {
      
      cat(cyan('counts must be in matrix or dgCMatrix format\n'))
      return(counts)
      
    }
    
  } 
  
  if(!is.null(z)) {
    
    if(length(z) != ncol(counts)) {
      
      cat(cyan('z must have the same length as ncol: counts\n'))
      return(counts)
      
    }
    
  }
  
  if(!is.numeric(maxIter)) {
    
    cat(cyan('maxIter must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(delta)) {
    
    cat(cyan('delta must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.logical(estimateDelta)) {
    
    cat(cyan('estimateDelta must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.numeric(convergence)) {
    
    cat(cyan('convergence must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(iterLogLik)) {
    
    cat(cyan('iterLogLik must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(varGenes)) {
    
    cat(cyan('varGenes must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(dbscanEps)) {
    
    cat(cyan('dbscanEps must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(seed)) {
    
    cat(cyan('seed must be numerical\n'))
    return(counts)
    
  }
  
  if(is.null(batch)) {
    
    d <- decontX(x = counts,
                        z = z,
                        batch = NULL,
                        maxIter = maxIter,
                        delta = delta,
                        estimateDelta = estimateDelta,
                        convergence = convergence,
                        iterLogLik = iterLogLik,
                        varGenes = varGenes,
                        dbscanEps = dbscanEps,
                        seed = seed,
                        verbose = TRUE)
    
  } else {
    
    d <- decontX(x = counts,
                        z = z,
                        batch = object@sample_metadata[,batch],
                        maxIter = maxIter,
                        delta = delta,
                        estimateDelta = estimateDelta,
                        convergence = convergence,
                        iterLogLik = iterLogLik,
                        varGenes = varGenes,
                        dbscanEps = dbscanEps,
                        seed = seed,
                        verbose = TRUE)
    
  }
  
  cat(cyan('Decontamination completed\n'))
  
  print(plotDecontXContamination(x = d))
  
  cat(cyan(paste0(formatC(sum(d$contamination)/length(d$contamination), 
                                  digits = 2), 
                          '% average contamination\n')))
  
  clean.matrix <- d$decontXcounts
  cat(cyan('Matrix isolated\n'))
  clean.matrix <- round(clean.matrix)
  zero.samples <- colSums(as.matrix(clean.matrix)) > 0
  clean.matrix <- clean.matrix[,zero.samples]
  cat(cyan('converted to integer\n'))
  return(clean.matrix)
  
}