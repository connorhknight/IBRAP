#' @name plot.variance
#' @aliases plot.variance
#' 
#' @title Plot reduction explained variance
#'
#' @description Shows the explained variance for a reduction for selecting how many dimensions to use downstream
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. Which reductions to access within the supplied assays
#' @param ... Arguments to be passed to ggpubr::ggarrange
#' 
#' @return Plots showing the supplied reductions explained variance
#' 
#' @examples 
#' 
#' object <- plot.variance(object = object, 
#'                         assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                         reduction = 'pca')
#'
#' @export

plot.variance <- function(object, assay, reduction, ...) {
  
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
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        stop(paste0('reduction: ', x, '  does not exist\n'))
        
      }
      
    }
    
  }
  
  list.of.plots <- list()
  
  count <- 1
  
  for(x in assay) {
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[x]]@computational_reductions), 
                   names(object@methods[[x]]@integration_reductions),
                   names(object@methods[[x]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[x]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[x]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[x]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[x]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[x]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[x]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        stop('reductions could not be found\n')
        
      }
      
    }
    
    for(p in reduction) {
      
      red <- reduction.list[[reduction]]
      
      sdev <- apply(red, MARGIN = 2, FUN = sd)
      
      eig <- sdev^2/sum(sdev^2)
      eig <- as.data.frame(eig*100)
      temp <- as.character(colnames(red))
      eig[,2] <- factor(x = temp, levels = unique(temp))
      colnames(eig) <- c ('Variance', 'PCs')
      
      list.of.plots[[count]] <- ggplot2::ggplot(data = eig, mapping = ggplot2::aes(x = PCs, y = Variance)) + 
        ggplot2::geom_point() + 
        egg::theme_article() + 
        ggplot2::ylab('Explained Variance (%)') +
        ggplot2::xlab('Components')
        ggplot2::ggtitle(paste0(x,'_',p)) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1), 
                       plot.title = ggplot2::element_text(hjust = 0.5))
      
      count <- count + 1
      
    }
    
  }
  
  ggpubr::ggarrange(plotlist = list.of.plots, ...)
  
}