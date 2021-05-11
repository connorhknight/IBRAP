#' @title Method override for `'['` subset function regarding IBRAP S4 object

setMethod(f = '[', signature = 'IBRAP',
          function(x, 
                   i, 
                   j, 
                   ..., 
                   drop = FALSE) {
            
            if(!missing(i) & missing(j)) {
              
              if(!is.character(i)) {
                
                ii <- .convert_subset_index(x = i, rownames(x))
                
              } else {
                
                ii <- i
                
              }
              
              .sample_metadata <- x@sample_metadata
              
              list.methods <- list()
              
              for(p in names(x@methods)) {
                
                if(length(as.matrix(x@methods[[p]]@counts)) != 0) {
                  
                  .counts <- x@methods[[p]]@counts[ii, , drop = FALSE]
                  
                } else {
                  
                  .counts <- x@methods[[p]]@counts
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@normalised)) != 0) {
                  
                  .normalised <- x@methods[[p]]@normalised[ii, , drop = FALSE]
                  
                } else {
                  
                  .normalised <- x@methods[[p]]@normalised
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@norm.scaled)) != 0) {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled[ii, , drop = FALSE]
                  
                } else {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@feature_metadata)) != 0) {
                  
                  .feature_metadata <- x@methods[[p]]@feature_metadata[ii, , drop = FALSE]
                  
                } else {
                  
                  .feature_metadata <- x@methods[[p]]@feature_metadata
                  
                }
                
                .highly.variable.genes <- x@methods[[p]]@highly.variable.genes
                
                .graphs <- x@methods[[p]]@graphs
                
                .computational_reductions <- x@methods[[p]]@computational_reductions
                
                .integration_reductions <- x@methods[[p]]@integration_reductions
                
                .visualisation_reductions <- x@methods[[p]]@visualisation_reductions
                
                .cluster_assignments <- x@methods[[p]]@cluster_assignments
                
                .benchmark_results <- x@methods[[p]]@benchmark_results
                
                .alt_objects <- x@methods[[p]]@alt_objects
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         graphs = .graphs,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         benchmark_results = .benchmark_results,
                                         alt_objects = .alt_objects)
                
              }
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata))
              
            } 
            
            if(missing(i) & !missing(j)) {
              
              if(!is.character(j)) {
                
                jj <- .convert_subset_index(x = j, colnames(x))
                
              } else {
                
                jj <- j
                
              }
              
              .sample_metadata <- x@sample_metadata[jj, , drop = FALSE]
              
              list.methods <- list()
              
              for(p in names(x@methods)) {
                
                if(length(as.matrix(x@methods[[p]]@counts)) != 0) {
                  
                  .counts <- x@methods[[p]]@counts[ , jj, drop = FALSE]
                  
                } else {
                  
                  .counts <- x@methods[[p]]@counts
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@normalised)) != 0) {
                  
                  .normalised <- x@methods[[p]]@normalised[ , jj, drop = FALSE]
                  
                } else {
                  
                  .normalised <- x@methods[[p]]@normalised
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@norm.scaled)) != 0) {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled[ , jj, drop = FALSE]
                  
                } else {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled
                  
                }
                
                if(length(x@methods[[p]]@computational_reductions) != 0) {
                  
                  .computational_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@computational_reductions)) {
                    
                    .computational_reductions[[g]] <- x@methods[[p]]@computational_reductions[jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .computational_reductions <-x@methods[[p]]@computational_reductions
                  
                }
                
                if(length(x@methods[[p]]@integration_reductions) != 0) {
                  
                  .integration_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@integration_reductions)) {
                    
                    .integration_reductions[[g]] <- x@methods[[p]]@integration_reductions[jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .integration_reductions <-x@methods[[p]]@integration_reductions
                  
                }
                
                if(length(x@methods[[p]]@visualisation_reductions) != 0) {
                  
                  .visualisation_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@visualisation_reductions)) {
                    
                    .visualisation_reductions[[g]] <- x@methods[[p]]@visualisation_reductions[jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .visualisation_reductions <-x@methods[[p]]@visualisation_reductions
                  
                }
                
                .highly.variable.genes <- x@methods[[p]]@highly.variable.genes
                
                .feature_metadata <- x@methods[[p]]@feature_metadata
                
                .graphs <- x@methods[[p]]@graphs
                
                .benchmark_results <- x@methods[[p]]@benchmark_results
                
                .alt_objects <- x@methods[[p]]@alt_objects
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         graphs = .graphs,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         benchmark_results = .benchmark_results,
                                         alt_objects = .alt_objects)
                
              }
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata))
              
            }
            
            if(!missing(i) & !missing(j)) {
              
              if(!is.character(i)) {
                
                ii <- .convert_subset_index(x = i, rownames(x))
                
              } else {
                
                ii <- i
                
              }
              
              if(!is.character(j)) {
                
                jj <- .convert_subset_index(x = j, colnames(x))
                
              } else {
                
                jj <- j
                
              }
              
              .sample_metadata <- x@sample_metadata[jj, , drop = FALSE]
              
              list.methods <- list()
              
              for(p in names(x@methods)) {
                
                if(length(as.matrix(x@methods[[p]]@counts)) != 0) {
                  
                  .counts <- x@methods[[p]]@counts[ii , jj, drop = FALSE]
                  
                } else {
                  
                  .counts <- x@methods[[p]]@counts
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@normalised)) != 0) {
                  
                  .normalised <- x@methods[[p]]@normalised[ii , jj, drop = FALSE]
                  
                } else {
                  
                  .normalised <- x@methods[[p]]@normalised
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@norm.scaled)) != 0) {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled[ii , jj, drop = FALSE]
                  
                } else {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled
                  
                }
                
                if(length(as.matrix(x@methods[[p]]@feature_metadata)) != 0) {
                  
                  .feature_metadata <- x@methods[[p]]@feature_metadata[ii, , drop = FALSE]
                  
                } else {
                  
                  .feature_metadata <- x@methods[[p]]@feature_metadata
                  
                }
                
                if(length(x@methods[[p]]@computational_reductions) != 0) {
                  
                  .computational_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@computational_reductions)) {
                    
                    .computational_reductions[[g]] <- x@methods[[p]]@computational_reductions[jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .computational_reductions <-x@methods[[p]]@computational_reductions
                  
                }
                
                if(length(x@methods[[p]]@integration_reductions) != 0) {
                  
                  .integration_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@integration_reductions)) {
                    
                    .integration_reductions[[g]] <- x@methods[[p]]@integration_reductions[jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .integration_reductions <-x@methods[[p]]@integration_reductions
                  
                }
                
                if(length(x@methods[[p]]@visualisation_reductions) != 0) {
                  
                  .visualisation_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@visualisation_reductions)) {
                    
                    .visualisation_reductions[[g]] <- x@methods[[p]]@visualisation_reductions[jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .visualisation_reductions <-x@methods[[p]]@visualisation_reductions
                  
                }
                
                .highly.variable.genes <- x@methods[[p]]@highly.variable.genes
                
                .graphs <- x@methods[[p]]@graphs
                
                .benchmark_results <- x@methods[[p]]@benchmark_results
                
                .alt_objects <- x@methods[[p]]@alt_objects
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         graphs = .graphs,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         benchmark_results = .benchmark_results,
                                         alt_objects = .alt_objects)
                
              }
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata))
            }
          })