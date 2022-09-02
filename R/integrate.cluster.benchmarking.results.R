#' @name integrate.cluster.benchmarking.results
#' @aliases integrate.cluster.benchmarking.results
#' 
#' @param object IBRAP S4 class object
#' 
#' @return data frame of all clustering results separated by pipeline
#' 
#'
#' @export integrate.cluster.benchmarking.results

integrate.cluster.benchmarking.results <- function(object) {
  
  new.list <- list()
  
  for(x in names(object@methods)) {
    
    if(length(object[[x]]@benchmark_results$clustering)!=0) {
      
      for(y in names(object[[x]]@benchmark_results$clustering)) {
        
        new.list[[paste0(x,'_',y)]] <- object[[x]]@benchmark_results$clustering[[y]]
        new.list[[paste0(x,'_',y)]]$normalisation_method <- x
        new.list[[paste0(x,'_',y)]]$pipeline <- y
        
      }
      
    }
    
  }
  
  return(do.call(what = rbind, args = new.list))
  
}
