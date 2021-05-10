#' @name plot.QC.scatter
#' @aliases plot.QC.scatter
#' 
#' @title Plots two QC metrices in scatter format
#'
#' @description Plots two QC metric in scatter format
#' 
#' @import crayon
#' @import egg
#' @import ggplots
#' @import RColorBrewer
#' 
#' @param object IBRAP S4 class object
#' @param x Character string of which metadata column to subset for the x axis 
#' @param y Character string of which metadata column to subset for the y axis 
#' @param split.by Character string of which metadata column to separate cells by
#' 
#' @usage plot.QC.vln(object = obj, metadata.columns = `c('RAW_total.features', 'RAW_total.counts')`, split.by = 'original.project')
#' 
#' @return A Violin plot
#'
#' @export

plot.QC.scatter <- function(object, 
                            x, 
                            y, 
                            split.by) {
  
  metadata <- object@sample_metadata
  
  if(!x %in% colnames(metadata)) {
    
    cat(cyan('X variable does not exist\n'))
    
  }
  
  if(!y %in% colnames(metadata)) {
    
    cat(cyan('Y variable does not exist\n'))
    
  }
  
  if(!is.null(split.by)) {
    
    if(!split.by %in% colnames(metadata)){
      
      cat(cyan('split.by variable does not exist\n'))
      
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
    
  }
  
  p <- ggplot(data = new.df, mapping = aes(x = x, y = y, col = project)) + 
    geom_point() + theme_classic() + ggtitle(paste0(x,'_vs_',y)) + 
    ylab(y) + xlab(x) + labs(color='identifier') + scale_color_manual(values=cols.proj)
  
  print(p)
  
}