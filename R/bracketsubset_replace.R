#' @name doublet_bracket_replace
#' 
#' @exportMethod `[[<-`

setMethod(f = '[[<-', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   value) {
            
            if(i %in% names(x@sample_metadata)) {
              
              x@sample_metadata[[i, exact = T]] <- value
              
            } else if(i %in% names(x@methods)) {
              
              x@methods[[i]] <- value
              
            } else if(i == 'pipelines'){
              
              x@pipelines <- value
              
            }
            
            x@sample_metadata[[i]] <- value
            
            return(x)
          })

setMethod(f = '[[<-', signature = 'methods',
          function(x, 
                   i, 
                   j, 
                   value){
            
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
            
            if(i %in% 'counts') {
              
              y@counts <- value
              
            }
            
            if(i %in% 'normalised') {
              
              y@normalised <- value
              
            }
            
            if(i %in% 'norm.scaled') {
              
              y@norm.scaled <- value
              
            }
            
            if(i %in% 'highly.variable.genes') {
              
              y@highly.variable.genes <- value
              
            }
            
            if(i %in% names(y@neighbours)) {
              
              y@neighbours[[i]] <- value
              
            }
            
            if(i %in% names(y@computational_reductions)) {
              
              y@computational_reductions[[i]] <- value
              
            }
            
            if(i %in% names(y@integration_reductions)) {
              
              y@integration_reductions[[i]] <- value
              
            }
            
            if(i %in% names(y@visualisation_reductions)) {
              
              y@visualisation_reductions[[i]] <- value
              
            }
            
            if(i %in% names(y@cluster_assignments)) {
              
              y@cluster_assignments[[i]] <- value
              
            }
            
            if(i %in% names(y@benchmark_results$clustering)) {
              
              y@benchmark_results$clustering[[i]] <- value
              
            }
            
            if(i %in% names(y@benchmark_results$integration)) {
              
              y@benchmark_results$integration[[i]] <- value
              
            }
            
          })
