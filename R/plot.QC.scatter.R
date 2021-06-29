#' @name plot.QC.scatter
#' @aliases plot.QC.scatter
#' 
#' @title Plots two QC metrices in scatter format
#'
#' @description Plots two QC metric in scatter format
#' 
#' @param object IBRAP S4 class object
#' @param x Character string of which metadata column to subset for the x axis 
#' @param y Character string of which metadata column to subset for the y axis 
#' @param split.by Character string of which metadata column to separate cells by
#' 
#' @usage plot.QC.scatter(object = obj, x = 'RAW_total.features', y = 'RAW_total.counts', split.by = 'original.project')
#' 
#' @return A Scatter plot
#'
#' @export plot.QC.scatter

plot.QC.scatter <- function(object, 
                            x, 
                            y, 
                            split.by) {
  
  metadata <- object@sample_metadata
  
  if(!x %in% colnames(metadata)) {
    
    cat(crayon::cyan('X variable does not exist\n'))
    return(NULL)
    
  }
  
  if(!y %in% colnames(metadata)) {
    
    cat(crayon::cyan('Y variable does not exist\n'))
    return(NULL)
    
  }
  
  if(!is.null(split.by)) {
    
    if(!split.by %in% colnames(metadata)){
      
      cat(crayon::cyan('split.by variable does not exist\n'))
      return(NULL)
      
    }
    
  }
  
  new.df <- data.frame(as.factor(rownames(metadata)))
  new.df$x <- metadata[,x]
  new.df$y <- metadata[,y]
  new.df$project <- metadata[,split.by]
  
  proj.length <- length(unique(new.df$project))
  
  if(proj.length < 3) {
    
    proj.length.new <- 3
    cols.proj <- RColorBrewer::brewer.pal(n = proj.length.new, name = 'Pastel1')
    cols.proj <- cols.proj[1:proj.length]
    
  } else {
    
    cols.proj <- RColorBrewer::brewer.pal(n = proj.length, name = 'Pastel1')
    
  }
  
  p <- ggplot2::ggplot(data = new.df, mapping = ggplot2::aes(x = x, y = y, col = project)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + ggplot2::ggtitle(paste0(x,'_vs_',y)) + 
    ggplot2::ylab(y) + ggplot2::xlab(x) + ggplot2::labs(color='identifier') + ggplot2::scale_color_manual(values=cols.proj)
  
  print(p)
  
}