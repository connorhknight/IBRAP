#' @name perform.pca
#' @aliases perform.pca
#' 
#' @title Performs PCA reduction
#'
#' @description Performs PCA reduction on defined method-assays. Data should be HVG subset, normalised and scaled (in the norm.scaled assay)
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param n.pcs Numerical. How many principal components should be produced. Default = 50
#' @param reduction.save Character. What should this reduction be saved as in computation_reduction. Default = 'pca'
#' @param print.variance Logical. Should the plot be printed to the console
#' @param ... Arguments to be passed to PCAtools::pca
#' 
#' @return PCA reductions contained within the computational_reduction list in the defined assays
#' 
#' @examples 
#' 
#' object <- perform.pca(object = object, 
#'                       assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                       n.pcs = 50, 
#'                       reduction.save = 'pca')
#'
#' @export

perform.pca <- function(object, 
                        assay,
                        slot='norm.scaled',
                        n.pcs=50,
                        reduction.save='pca', 
                        print.variance = FALSE, 
                        ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('assay: ', x, ' does not exist\n'))
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be character string')
    
  }
  
  if(!is.numeric(n.pcs)) {
    
    stop('n.pcs must be numerical')
    
  }
  
  if(!is.character(reduction.save)) {
    
    stop('reduction.save must be numerical')
    
  }
  
  if(!is.logical(print.variance)) {
    
    stop('print.plot must be logical, TRUE/FALSE \n')
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  list.of.figs <- list()
  
  for(t in assay) {
    
    if(is.null(object@methods[[t]]@highly.variable.genes)) {
      
      stop(paste0('no variable features have been identified for assay: ', t))
      
    }
    
    mat <- as.matrix(object@methods[[t]][[slot]][rownames(object@methods[[t]][[slot]]) %in% object@methods[[t]]@highly.variable.genes,])
    
    cat(crayon::cyan(paste0(Sys.time(), ': initialising PCA for assay:', t, '\n')))
    
    ass <- strsplit(x = names(object@methods)[which(names(object@methods)==t)], split = '_')[[1]][1]
    
    if(ass %in% c('SCT','SCRAN','TPM')) {
      
      a <- suppressWarnings(PCAtools::pca(mat = mat, center = F, scale = F, ...))

      cat(crayon::cyan(paste0(Sys.time(), ': PCA completed\n')))
      
      object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,1:n.pcs])
      
    } else if (ass == 'SCANPY') {
      
      sc <- reticulate::import('scanpy')
      
      scobj <- sc$AnnData(X = t(as.matrix(object@methods[[t]][['norm.scaled']])))
      scobj$obs_names <- as.factor(colnames(object@methods[[t]][['norm.scaled']]))
      scobj$var_names <- as.factor(rownames(object@methods[[t]][['norm.scaled']]))
      
      sc$tl$pca(data = scobj, n_comps = as.integer(n.pcs), use_highly_variable = as.logical(F))
      
      tmp <- scobj$obsm[['X_pca']]
      rownames(tmp) <- colnames(object)
      
      pc.names <- list()
      count <- 1
      
      for(x in 1:n.pcs) {
        
        pc.names[[count]] <- paste0('PC',count)
        
        count <- count + 1
        
      }
      
    colnames(tmp) <- unlist(pc.names)
    rownames(tmp) <- colnames(object@methods[[t]][['norm.scaled']])
    
    cat(crayon::cyan(paste0(Sys.time(), ': PCA completed\n')))
    
    object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(tmp)
      
    } else {
      
      a <- suppressWarnings(PCAtools::pca(mat = mat, center = F, scale = F, ...))

      cat(crayon::cyan(paste0(Sys.time(), ': PCA completed\n')))
      
      object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,1:n.pcs])
      
    }

  }
  
  if(isTRUE(print.variance)) {
    
    print(IBRAP::plot.variance(object = object, assay = assay, reduction = reduction.save))
    
  }
  
  return(object)
  
}