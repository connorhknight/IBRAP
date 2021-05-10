#' @name plot.QC.vln
#' @aliases plot.QC.vln
#' 
#' @title Plots a given QC metric
#'
#' @description Plots a defined QC metric in violin format
#' 
#' @import crayon
#' @import egg
#' @import ggplots
#' @import RColorBrewer
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
    
    cat(cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  plots.list <- list()
  metadata <- object@sample_metadata
  
  for(m in metadata.columns) {
    
    if(!m %in% colnames(metadata)) {
      
      cat(cyan('Provided column names do not exist\n'))
      return(NULL)
      
    }
    
  }
  
  if(!split.by %in% colnames(object@sample_metadata)) {
    
    cat(cyan(paste0(split.by, ' does not exist\n')))
    return(object)
    
  }
  
  cols <- brewer.pal(n = length(metadata.columns), name = 'Pastel2')
  
  count <- 1
  
  for(o in metadata.columns) {
    
    new.metadata <- data.frame(project=as.factor(object[[split.by]]))
    new.metadata$sample <- as.factor(colnames(object))
    new.metadata$variable <- object[[o]]
    proj.length <- length(unique(new.metadata$project))
    
    if(proj.length < 3) {
      
      proj.length.new <- 3
      cols.proj <- brewer.pal(n = proj.length.new, name = 'Pastel1')
      cols.proj <- cols.proj[1:proj.length]
      
    }
    
    plots.list[[o]] <- ggplot(data = new.metadata, 
                                       mapping = aes(x=variable, y=project, fill=project)) + 
      geom_violin() + coord_flip() + ggtitle(o) + 
      xlab('') + ylab('project') + theme_classic() + 
      geom_boxplot(lwd = 0.6, width = 0.09, fill = cols[[count]]) +
      theme(axis.text.x = element_text(face = 'bold', angle = 45, vjust = 1, hjust=1), 
                     legend.position="none", plot.title = element_text(hjust=0.5)) + 
      scale_fill_manual(values=cols.proj)
    count <- count + 1
  }
  
  do.call(what = 'ggarrange', args = list(plots = plots.list, nrow=1, ncol=length(plots.list)))
  
}