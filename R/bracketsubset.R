#' @name doublet_bracket
#' 
#' @exportMethod `[[`

setMethod(f = '[[', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   ...) {
            x@sample_metadata[[i, exact = TRUE]]
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
            y[[i]]
            
          })