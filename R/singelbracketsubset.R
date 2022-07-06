#' @title Method override for `'['` subset function regarding IBRAP S4 object
#' 
#' @exportMethod `[`

setMethod(f = '[', signature = 'IBRAP',
          function(x, 
                   i, 
                   j, 
                   ..., 
                   drop = FALSE) {
            
            if(!missing(i) & missing(j)) {
              
              # Just features
              
              if(!is.character(i)) {
                
                ii <- .convert_subset_index(x = i, rownames(x))
                
              } else {
                
                ii <- i
                
              }
              
              .sample_metadata <- x@sample_metadata
              
              list.methods <- list()
              
              for(p in names(x@methods)) {
                
                if(is.na(x@methods[[p]]@counts[1])) {

                  genes <- ii[ii %in% rownames(x@methods[[p]]@counts)]
                  .counts <- x@methods[[p]]@counts[genes , , drop = FALSE]
                  
                } else {
                  
                  .counts <- x@methods[[p]]@counts
                  
                }
                
                if(is.na(x@methods[[p]]@normalised[1])) {

                  genes <- ii[ii %in% rownames(x@methods[[p]]@normalised)]
                  .normalised <- x@methods[[p]]@normalised[genes, ,drop = FALSE]
                  
                } else {
                  
                  .normalised <- x@methods[[p]]@normalised
                  
                }
                
                if(is.na(x@methods[[p]]@norm.scaled[1])) {
                  
                  genes <- ii[ii %in% rownames(x@methods[[p]]@norm.scaled)]
                  .norm.scaled <- x@methods[[p]]@norm.scaled[genes, ,drop = FALSE]
                  
                } else {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled
                  
                }
                
                if(is.na(x@methods[[p]]@feature_metadata[1])) {
                  
                  genes <- ii[ii %in% rownames(x@methods[[p]]@feature_metadata)]
                  .feature_metadata <- x@methods[[p]]@feature_metadata[genes, , drop = FALSE]
                  
                } else {
                  
                  .feature_metadata <- x@methods[[p]]@feature_metadata
                  
                }
                
                if(!is.null(x@methods[[p]]@highly.variable.genes)) {
                  
                  .highly.variable.genes <- x@methods[[p]]@highly.variable.genes[x@methods[[p]]@highly.variable.genes %in% rownames(.counts)]
                  
                } else {
                  
                  .highly.variable.genes <- x@methods[[p]]@highly.variable.genes
                  
                }

                .neighbours <- x@methods[[p]]@neighbours
                
                .computational_reductions <- x@methods[[p]]@computational_reductions
                
                .integration_reductions <- x@methods[[p]]@integration_reductions
                
                .visualisation_reductions <- x@methods[[p]]@visualisation_reductions
                
                .cluster_assignments <- x@methods[[p]]@cluster_assignments
                
                .benchmark_results <- x@methods[[p]]@benchmark_results
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         neighbours = .neighbours,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         cluster_assignments = .cluster_assignments,
                                         benchmark_results = .benchmark_results)
                
              }
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata))
              
            } 
            
            if(missing(i) & !missing(j)) {
              
              # Just cells
              
              if(!is.character(j)) {
                
                jj <- .convert_subset_index(x = j, colnames(x))
                
              } else {
                
                jj <- j
                
              }
              
              .sample_metadata <- x@sample_metadata[jj, , drop = FALSE]
              
              list.methods <- list()
              
              for(p in names(x@methods)) {
                
                if(is.na(x@methods[[p]]@counts[1])) {
                  
                  .counts <- x@methods[[p]]@counts[ , jj, drop = FALSE]
                  
                } else {
                  
                  .counts <- x@methods[[p]]@counts
                  
                }
                
                if(is.na(x@methods[[p]]@normalised[1])) {
                  
                  .normalised <- x@methods[[p]]@normalised[ , jj, drop = FALSE]
                  
                } else {
                  
                  .normalised <- x@methods[[p]]@normalised
                  
                }
                
                if(is.na(x@methods[[p]]@norm.scaled[1])) {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled[ , jj, drop = FALSE]
                  
                } else {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled
                  
                }
                
                if(length(x@methods[[p]]@computational_reductions) != 0) {
                  
                  .computational_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@computational_reductions)) {
                    
                    .computational_reductions[[g]] <- x@methods[[p]]@computational_reductions[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .computational_reductions <-x@methods[[p]]@computational_reductions
                  
                }
                
                if(length(x@methods[[p]]@integration_reductions) != 0) {
                  
                  .integration_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@integration_reductions)) {
                    
                    .integration_reductions[[g]] <- x@methods[[p]]@integration_reductions[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .integration_reductions <-x@methods[[p]]@integration_reductions
                  
                }
                
                if(length(x@methods[[p]]@visualisation_reductions) != 0) {
                  
                  .visualisation_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@visualisation_reductions)) {
                    
                    .visualisation_reductions[[g]] <- x@methods[[p]]@visualisation_reductions[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .visualisation_reductions <- x@methods[[p]]@visualisation_reductions
                  
                }
                
                if(length(x@methods[[p]]@cluster_assignments) != 0) {
                  
                  .cluster_assignments <- list()
                  
                  for(g in names(x@methods[[p]]@cluster_assignments)) {
                    
                    .cluster_assignments[[g]] <- x@methods[[p]]@cluster_assignments[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .cluster_assignments <- x@methods[[p]]@cluster_assignments
                  
                }
                
                .highly.variable.genes <- x@methods[[p]]@highly.variable.genes
                
                .feature_metadata <- x@methods[[p]]@feature_metadata
                
                if(length(x@methods[[p]]@neighbours) != 0) {
                  
                  list.neighbours <- list()
                  
                  for(l in names(x@methods[[p]]@neighbours)) {
                    
                    .neighbours <- x@methods[[p]]@neighbours
                    
                    list.neighbours.sub <- list()

                    for(t in names(x@methods[[p]]@neighbours[[l]])) {

                      list.neighbours.sub[[t]] <- x@methods[[p]]@neighbours[[l]][[t]][jj, jj, drop = FALSE]
                      
                    }
                    
                    list.neighbours[[l]] <- list.neighbours.sub
                    
                  }
                  
                  .neighbours <- list.neighbours
                  
                } else {
                  
                  .neighbours <- x@methods[[p]]@neighbours
                  
                }
                
                .benchmark_results <- x@methods[[p]]@benchmark_results
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         neighbours = .neighbours,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         cluster_assignments = .cluster_assignments,
                                         benchmark_results = .benchmark_results)
                
              }
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata))
              
            }
            
            if(!missing(i) & !missing(j)) {
              
              # features and cells
              
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
                
                if(is.na(x@methods[[p]]@counts[1])) {

                  cells <- jj[jj %in% colnames(x@methods[[p]]@counts)]
                  genes <- ii[ii %in% rownames(x@methods[[p]]@counts)]
                  .counts <- x@methods[[p]]@counts[genes , cells, drop = FALSE]

                } else {
                  
                  .counts <- x@methods[[p]]@counts
                  
                }
                
                if(is.na(x@methods[[p]]@normalised[1])) {
                  
                  cells <- jj[jj %in% colnames(x@methods[[p]]@normalised)]
                  genes <- ii[ii %in% rownames(x@methods[[p]]@normalised)]
                  .normalised <- x@methods[[p]]@normalised[genes , cells, drop = FALSE]
                  
                } else {
                  
                  .normalised <- x@methods[[p]]@normalised
                  
                }
                
                
                if(is.na(x@methods[[p]]@norm.scaled[1])) {
                  
                  cells <- jj[jj %in% colnames(x@methods[[p]]@norm.scaled)]
                  genes <- ii[ii %in% rownames(x@methods[[p]]@norm.scaled)]
                  .norm.scaled <- x@methods[[p]]@norm.scaled[genes , cells, drop = FALSE]
                  
                } else {
                  
                  .norm.scaled <- x@methods[[p]]@norm.scaled
                  
                }

                if(is.na(x@methods[[p]]@feature_metadata[1])) {
                  
                  genes <- ii[ii %in% rownames(x@methods[[p]]@feature_metadata)]
                  .feature_metadata <- x@methods[[p]]@feature_metadata[genes, , drop = FALSE]
                  
                } else {
                  
                  .feature_metadata <- x@methods[[p]]@feature_metadata
                  
                }
                
                if(length(x@methods[[p]]@computational_reductions) != 0) {
                  
                  .computational_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@computational_reductions)) {
                    
                    .computational_reductions[[g]] <- x@methods[[p]]@computational_reductions[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .computational_reductions <-x@methods[[p]]@computational_reductions
                  
                }
                
                if(length(x@methods[[p]]@integration_reductions) != 0) {
                  
                  .integration_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@integration_reductions)) {
                    
                    .integration_reductions[[g]] <- x@methods[[p]]@integration_reductions[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .integration_reductions <-x@methods[[p]]@integration_reductions
                  
                }
                
                if(length(x@methods[[p]]@visualisation_reductions) != 0) {
                  
                  .visualisation_reductions <- list()
                  
                  for(g in names(x@methods[[p]]@visualisation_reductions)) {
                    
                    .visualisation_reductions[[g]] <- x@methods[[p]]@visualisation_reductions[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .visualisation_reductions <-x@methods[[p]]@visualisation_reductions
                  
                }
                
                if(length(x@methods[[p]]@cluster_assignments) != 0) {
                  
                  .cluster_assignments <- list()
                  
                  for(g in names(x@methods[[p]]@cluster_assignments)) {
                    
                    .cluster_assignments[[g]] <- x@methods[[p]]@cluster_assignments[[g]][jj, , drop = FALSE]
                    
                  }
                  
                } else {
                  
                  .cluster_assignments <- x@methods[[p]]@cluster_assignments
                  
                }
                
                if(!is.null(x@methods[[p]]@highly.variable.genes)) {
                  
                  .highly.variable.genes <- x@methods[[p]]@highly.variable.genes[x@methods[[p]]@highly.variable.genes %in% rownames(.counts)]
                  
                } else {
                  
                  .highly.variable.genes <- x@methods[[p]]@highly.variable.genes
                  
                }
                
                if(length(x@methods[[p]]@neighbours) != 0) {
                  
                  list.neighbours <- list()
                  
                  for(l in names(x@methods[[p]]@neighbours)) {
                    
                    .neighbours <- x@methods[[p]]@neighbours
                    
                    list.neighbours.sub <- list()
                    
                    for(t in names(x@methods[[p]]@neighbours[[l]])) {
                      
                      list.neighbours.sub[[t]] <- x@methods[[p]]@neighbours[[l]][[t]][jj, jj, drop = FALSE]
                      
                    }
                    
                    list.neighbours[[l]] <- list.neighbours.sub
                    
                  }
                  
                  .neighbours <- list.neighbours
                  
                } else {
                  
                  .neighbours <- x@methods[[p]]@neighbours
                  
                }
                
                .benchmark_results <- x@methods[[p]]@benchmark_results
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         neighbours = .neighbours,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         cluster_assignments = .cluster_assignments,
                                         benchmark_results = .benchmark_results)
                
              }
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata))
            }
          })
