#' @name doublet_bracket
#' 
#' @exportMethod `[[`

setMethod(f = '[[', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   ...) {
            
            if(i %in% names(x@sample_metadata)) {
              
              x@sample_metadata[[i, exact = T]]
              
            } else if(i %in% names(x@methods)) {
              
              x@methods[[i, exact = T]]
              
            } else if(i == 'pipelines'){
              
              x@pipelines
              
            }
            
          })

setMethod(f = '[[', signature = 'methods',
          function(x, 
                   i, 
                   j, 
                   ...) {
            
            as.list.methods <- function(x) {
              
              new.list <- list(counts = x@counts, 
                               normalised = x@normalised, 
                               norm.scaled = x@norm.scaled,
                               highly.variable.genes = x@highly.variable.genes,
                               feature_metadata = x@feature_metadata,
                               neighbours = x@neighbours,
                               computational_reductions = x@computational_reductions,
                               integration_reductions = x@integration_reductions,
                               visualisation_reductions = x@visualisation_reductions,
                               cluster_assignments = x@cluster_assignments,
                               benchmark_results = x@benchmark_results)
              return(new.list)
              
            }
            
            y <- as.list.methods(x)
            
            if(i == 'counts') {
              
              return(y$counts)
              
            }
            
            if(i == 'normalised') {
              
              return(y$normalised)
              
            }
            
            if(i == 'norm.scaled') {
              
              return(y$norm.scaled)
              
            }
            
            if(i == 'highly.variable.genes') {
              
              return(y$highly.variable.genes)
              
            }
            
            if(i %in% names(y$neighbours)) {
              
              return(y$neighbours[[i]])
              
            }
            
            if(i %in% names(y$computational_reductions)) {
              
              return(y$computational_reductions[[i]])
              
            }
            
            if(i %in% names(y$integration_reductions)) {
              
              return(y$integration_reductions[[i]])
              
            }
            
            if(i %in% names(y$visualisation_reductions)) {
              
              return(y$visualisation_reductions[[i]])
              
            }
            
            if(i %in% names(y$cluster_assignments)) {
              
              return(y$cluster_assignments[[i]])
              
            }
            
            if(i %in% names(y$benchmark_results$clustering)) {
              
              return(y$benchmark_results$clustering[[i]])
              
            }
            
            if(i %in% names(y$benchmark_results$integration)) {
              
              return(y$benchmark_results$integration[[i]])
              
            }
            
          })
