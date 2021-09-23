#' @name plot.features
#' @aliases plot.features
#' 
#' @title Plot of reduced dimensions and features
#'
#' @param object An IBRAP S4 class object
#' @param assay Character. Which assay within the object to access
#' @param slot Character. Which expression matrix would you like to access. Default = 'normalised'
#' @param reduction Character. Which reduction should be used for the plot
#' @param features Character Which features should be plotted
#' @param order Boolean. Should datapoints be plotted in order of expression intensity. Default = TRUE
#' @param percentile. Numerical. What percentile of datapoint expression should be plotted. Default = c(0.1,0.9)
#' @param pt_size Numeric. what size should the inidividual plot sizes be
#' @param cells Numeric. Which cells should be subset for plotting, Default = NULL
#' 
#' @return A plot of reduced dimensions annotated with assignments
#'
#' @export plot.features

plot.features <- function(object, 
                          assay, 
                          slot,
                          reduction,
                          features,
                          order = TRUE,
                          percentile = c(0.1,0.9),
                          pt_size = 0.5,
                          cells = NULL) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(reduction)) {
    
    stop('reduction must be character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('assay does not exist\n')
    
  }
  
  if(!reduction %in% names(c(object@methods[[assay]]@computational_reductions, 
                             object@methods[[assay]]@visualisation_reductions, 
                             object@methods[[assay]]@integration_reductions))) {
    
    stop('reduction does not exist\n')
    
  }
  
  if(!is.character(features)) {
    
    stop('features must be character string \n')
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  plot.list <- list()
  
  for(x in features) {
    
    results <- as.data.frame(object@methods[[assay]]@visualisation_reductions[[reduction]])[,1:2]
    
    orig.colnames <- colnames(object@methods[[assay]]@visualisation_reductions[[reduction]][,1:2])

    if(is.null(object@methods[[assay]][[slot]][x,])) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': feature ', x, ' was not present in the defined assay, resorting to counts matrix', '\n')))

      iso <- object@methods[[1]]@counts[x,]
      
    } else {
      
      iso <- object@methods[[assay]][[slot]][x,]
      
    }
    
    colnames(results) <- c('red_1', 'red_2')
    
    results[,x] <- iso
    
    colnames(results)[3] <- 'feature'
    
    if(!is.null(cells)) {
      
      results <- results[cells,]
      
    }
    
    if(isTRUE(order)) {
      
      results <- results[order(as.numeric(results$feature)),]
      
    }
    
    lower <- as.numeric(quantile(results$feature, percentile)[1])
    upper <- as.numeric(quantile(results$feature, percentile)[2])
    
    results$feature[which(results$feature <= upper)] <- 0
    results$feature[which(results$feature >= lower)] <- 0

    plot.list[[x]] <- ggplot2::ggplot(data = results[order(results$feature),], 
                                      ggplot2::aes(x = red_1, y = red_2)) + 
      ggplot2::geom_point(ggplot2::aes(color=feature), size = pt_size)+ 
      ggplot2::scale_color_gradient(low = '#C0C0C0', high = '#4169E1') + 
      ggplot2::theme_classic() + ggplot2::labs(title=x, x=orig.colnames[1], y=orig.colnames[2]) + 
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = 'bold', size = 20))

  }
  
  is.even <- function(x) {
    
    if(as.integer(x %% 2) == 0) {
      
      return(TRUE)
      
    } else {
      
      return(FALSE)
      
    }
    
  }
  
  if(!is.even(length(plot.list))) {
    
    plot.list[length(plot.list)+1] <- plot.list[1] + ggplot2::geom_blank()
    
  }
  
  if(length(plot.list) > 1){

    if(length(plot.list) <3) {

      do.call('ggarrange.tmp', c(plots = plot.list, ncol = 1, nrow = 2))
      
    } else if(length(plot.list) <5) {

      do.call('ggarrange.tmp', c(plots = plot.list, ncol = 2, nrow = 2))
      
    } else if(length(plot.list) <7) {

      do.call('ggarrange.tmp', c(plots = plot.list, ncol = 3, nrow = 2))
      
    } else {
      
      do.call('ggarrange.tmp', c(plots = plot.list))
      
    }
    
  } else {

    plot.list[1]
    
  }
  
}
