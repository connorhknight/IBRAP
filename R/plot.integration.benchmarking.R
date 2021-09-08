#' @name plot.integration.benchmarking
#' @aliases plot.integration.benchmarking
#' 
#' @param object An IBRAP S4 class object
#' @param assay Character. Which assay within the object to access
#' 
#' @export plot.integration.benchmarking

plot.integration.benchmarking <- function(object, assay) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }

  for(x in assay) {
    
    if(is.null(object@methods[[x]]@benchmark_results$integration)) {
      
      stop(paste0('integration benchmarking has not been performed on this assay: ', x))
      
    }
    
    if(!is.list(object@methods[[x]]@benchmark_results$integration)) {
      
      stop('integration benchmarking results must be in list format')
      
    }
    
  }
  
  results <- data.frame(state = character(), ASW = numeric(), assay = character())
  
  count <- 1
  
  for(x in assay) {
    
    res <- object@methods[[x]]@benchmark_results$integration
    
    for(i in 1:length(res)) {
      
      results[count,'state'] <- strsplit(x = names(res)[[i]], split = '_')[[1]][2]
      results[count,'ASW'] <- res[[i]]
      results[count,'assay'] <- x
      
      count <- count + 1
      
    }
    
  }
  
  p <- ggplot2::ggplot(results, ggplot2::aes(fill=state, y=ASW, x=assay)) + 
    ggplot2::geom_bar(position="dodge", stat="identity") +
    ggplot2::theme_classic()
  
  print(p)
  
}
