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
#' @param ... Arguments to be passed to PCAtools::pca
#' 
#' @return PCA reductions contained within the computational_reduction list in the defined assays
#'
#' @export

perform.pca <- function(object, 
                        assay,
                        slot='norm.scaled',
                        n.pcs=50,
                        reduction.save='pca', 
                        ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('assay: ', x, ' does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be character string'))
    return(object)
    
  }
  
  if(!is.numeric(n.pcs)) {
    
    cat(crayon::cyan('n.pcs must be numerical'))
    return(object)
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be numerical'))
    return(object)
    
  }
  
  for(t in assay) {
    
    mat <- object@methods[[t]][[slot]]
    cat(crayon::cyan('Initialising PCA for assay:', t, '\n'))
    a <- PCAtools::pca(mat = mat, center = F, scale = F, ...)
    b <- PCAtools::findElbowPoint(a$variance)
    
    p <- PCAtools::screeplot(pcaobj = a, components = 1:sum(as.numeric(b)+10), 
                             title = paste0(assay,'_PCA_variance'), vline = b) +
      ggplot2::geom_label(ggplot2::aes(x = b, y = 50,
                                       label = 'Elbow point', vjust = -1, size = 8)) +
      ggplot2::ggtitle(paste0(t,'_screeplot'))
    
    print(p)
    
    cat(crayon::cyan('PCA completed\n'))
    
    object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,n.pcs])
    
  }
  
  return(object)
  
}