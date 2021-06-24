#' @name double_bracket-methods
#' 
#' @title Method override for `'[[<-'` subset function regarding methods S4 object
#' 
#'  @exportMethod  

setMethod(f = '[[<-', signature = 'methods',
          function(x, 
                   i, 
                   j, 
                   value){
            
            y <- as.list(x)
            y[[i]] <- value
            
            x@counts <- y$counts
            if(length(as.matrix(y$normalised)) != 0) {
              
              x@normalised <- y$normalised
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@norm.scaled <- y$norm.scaled
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@highly.variable.genes <- y$highly.varaible.genes
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@feature_metadata <- y$feature_metadata
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@neighbours <- y$neighbours
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@computational_reductions <- y$computational_reductions
              
            }
            if(length(as.matrix(y$norm.scaled))) {
              
              x@integration_reductions <- y$integration_reductions
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@visualisation_reductions <- y$visualisation_reductions
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@cluster_assignments <- y$cluster_assignments
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@benchmark_results <- y$benchmark_results
              
            }
            
            return(x)
            
          })
