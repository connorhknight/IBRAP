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
#' @param plot.var Boolean. Should the explained variance of PCs be plotted
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
                        plot.var = TRUE,
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
    eig <- a$sdev^2/sum(a$sdev^2)
    eig <- as.data.frame(eig*100)
    temp <- as.character(colnames(a$rotated))
    eig[,2] <- factor(x = temp, levels = unique(temp))
    colnames(eig) <- c ('Variance', 'PCs')

    p <- ggplot2::ggplot(data = eig[n.pcs,], mapping = ggplot2::aes(x = PCs, y = Variance)) + ggplot2::geom_point() + egg::theme_article() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) + ggplot2::ylab('Explained Variance (%)')
    
    print(p)
    
    cat(crayon::cyan('PCA completed\n'))
    
    object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,n.pcs])
    
  }
  
  return(object)
  
}