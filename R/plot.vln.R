#' @name plot.vln
#' @aliases plot.vln
#' 
#' @title Plot of violin plot of defined features
#'
#' @param object An IBRAP S4 class object
#' @param assay Character. Which assay within the object to access
#' @param slot Character. Which expression matrix would you like to access. Default = 'normalised'
#' @param features Character. Which features should be plotted, it is not recommended to exceed 3 but can be more if your screen size is larger. 
#' @param group.by Vector. What should the cell types be split by. 
#' @param title Character. What should be the overall title of the plot. Default = NULL
#' @param xlab Character. What should by the x axis title. Default = 'group'
#' @param ylab Character. What shoul be the y axis title. Default = 'expression'
#' 
#' @return A violin plot of defined gene expression and groups
#'
#' @export plot.features

plot.vln <- function(object, 
                     assay,
                     slot,
                     features, 
                     group.by, 
                     title = NULL, 
                     xlab = 'group', 
                     ylab = 'expression') {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string \n')
    
  } else if (is.character(assay)) {
    
    if(!assay %in% names(object@methods)) {
      
      stop('assay is not contained within object@methods \n')
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be character string \n')
    
  } else if(is.character(slot)) {
    
    if(!slot %in% c('counts','normalised','norm.scaled')) {
      
      stop('slot must be either counts, normalised or norm.scaled \n')
      
    }
    
  }
  
  for(x in features) {
    
    if(!x %in% rownames(object@methods[[assay]][[slot]])) {
      
      stop(paste0(x, ' not contained within expression matrix \n'))
      
    }
    
  }
  
  if(length(group.by) != ncol(object)) {
    
    stop('the group.by vector length does not \n')
    
  }
  
  if(!is.character(title)) {
    
    stop('title must be character string \n')
    
  }
  
  if(!is.character(xlab)) {
    
    stop('xlab must be character string \n')
    
  }
  
  if(!is.character(ylab)) {
    
    stop('ylab must be character string \n')
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  plot.list <- list()
  
  count <- 1
  
  for(x in features) {
    ass <- t(object@methods[[assay]][[slot]][x,])
    df <- data.frame(barcodes = rownames(object@sample_metadata))
    rownames(df) <- df$barcodes
    df[,'feature'] <- t(ass)
    df[,'group'] <- group.by
    
    if(count < length(features)) {
      
      p <- ggplot2::ggplot(data = df, ggplot2::aes(x = group, y = feature, color = group)) + 
        ggplot2::geom_violin() + ggplot2::geom_boxplot() +
        ggplot2::theme_light() + ggplot2::labs(title = x, y = as.character(ylab), x = '') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold', size = 20),
                       legend.position = 'none', axis.text.x = ggplot2::element_blank())
      
    } else {
      
      p <- ggplot2::ggplot(data = df, ggplot2::aes(x = group, y = feature, color = group)) + 
        ggplot2::geom_violin() + ggplot2::geom_boxplot() +
        ggplot2::theme_light() + ggplot2::labs(title = x, y = as.character(ylab), x = as.character(xlab)) + 
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold', size = 20), 
                       axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1),
                       legend.position = 'none')
      
    }
    
    plot.list[[x]] <- p
    
    count <- count + 1
    
  }
  
  if(length(plot.list) > 1) {
    do.call('ggarrange.tmp', c(plots = plot.list, ncol = 1))
  } else {
    plot.list[1]
  }
}
