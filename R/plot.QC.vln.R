#' @name plot.QC.vln
#' @aliases plot.QC.vln
#' 
#' @title Plots a given QC metric
#'
#' @description Plots a defined QC metric in violin format
#' 
#' @param object IBRAP S4 class object
#' @param metadata.columns A character string or vector of character strings of which metadata columns to subset and plot
#' @param split.by A character string indicating which metadata column should be used to divide samples
#' 
#' @usage plot.QC.vln(object = obj, metadata.columns = `c('RAW_total.features', 'RAW_total.counts')`, split.by = 'original.project')
#' 
#' @return A Violin plot
#'
#' @export

plot.QC.vln <- function(object, 
                        metadata.columns=c('RAW_total.features', 
                                           'RAW_total.counts'), 
                        split.by='original.project') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  plots.list <- list()
  metadata <- object@sample_metadata
  
  for(m in metadata.columns) {
    
    if(!m %in% colnames(metadata)) {
      
      cat(crayon::cyan('Provided column names do not exist\n'))
      return(NULL)
      
    }
    
  }
  
  if(!split.by %in% colnames(object@sample_metadata)) {
    
    cat(crayon::cyan(paste0(split.by, ' does not exist\n')))
    return(object)
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  cols <- RColorBrewer::brewer.pal(n = length(metadata.columns), name = 'Pastel2')
  
  count <- 1
  
  for(o in metadata.columns) {
    
    new.metadata <- data.frame(project=as.factor(object[[split.by]]))
    new.metadata$sample <- as.factor(colnames(object))
    new.metadata$variable <- object[[o]]
    proj.length <- length(unique(new.metadata$project))
    
    if(proj.length < 3) {
      
      proj.length.new <- 3
      cols.proj <- RColorBrewer::brewer.pal(n = proj.length.new, name = 'Pastel1')
      cols.proj <- cols.proj[1:proj.length]
      
    }
    
    plots.list[[o]] <- ggplot2::ggplot(data = new.metadata, 
                                       mapping = ggplot2::aes(x=variable, y=project, fill=project)) + 
      ggplot2::geom_violin() + ggplot2::coord_flip() + ggplot2::ggtitle(o) + 
      ggplot2::xlab('') + ggplot2::ylab('project') + ggplot2::theme_classic() + 
      ggplot2::geom_boxplot(lwd = 0.6, width = 0.09, fill = cols[[count]]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(face = 'bold', angle = 45, vjust = 1, hjust=1), 
                     legend.position="none", plot.title = ggplot2::element_text(hjust=0.5)) + 
      ggplot2::scale_fill_manual(values=cols.proj)
    count <- count + 1
  }
  
  do.call(what = 'ggarrange.tmp', args = list(plots = plots.list, nrow=1, ncol=length(plots.list)))
  
}
