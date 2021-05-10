#' @title Method override for as.list function regarding methods S4 object
#'
#' @export

setMethod(f = 'as.list', signature = 'methods',
          function(x) {
            
            new.list <- list(counts = x@counts, 
                             normalised = x@normalised, 
                             norm.scaled = x@norm.scaled,
                             highly.variable.genes = x@highly.variable.genes,
                             feature_metadata = x@feature_metadata,
                             graphs = x@graphs,
                             computational_reductions = x@computational_reductions,
                             integration_reductions = x@integration_reductions,
                             visualisation_reductions = x@visualisation_reductions,
                             cluster_assignments = x@cluster_assignments,
                             benchmark_results = x@benchmark_results,
                             alt_objects = x@alt_objects)
            return(new.list)
            
          })