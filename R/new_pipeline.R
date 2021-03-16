data.obj <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/sce_full/sce_full_Koh.rds')
kumar.obj <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/sce_full/sce_full_Kumar.rds')
mat <- assay(data.obj,'counts')
mat <- round(mat)

reticulate::use_python('/Users/knight05/Library/r-miniconda/envs/r-reticulate/bin/python')
reticulate::import('scrublet', convert = FALSE)

# Read 10x

require(Biobase)

Read10X_output <- function(directory, 
                           matrix.file = 'matrix.mtx', 
                           genes.file = 'genes.tsv', 
                           barcodes.file = 'barcodes.tsv',
                           make.feat.unique = TRUE) {
  
  dir.files <- list.files(path = directory)
  if(!matrix.file %in% dir.files) {
    cat(crayon::cyan(paste0('Expected file: ', matrix.file, '\n')))
  }
  if(!genes.file %in% dir.files) {
    cat(crayon::cyan('Expected file: ', genes.file, '\n'))
  }
  if(!barcodes.file %in% dir.files) {
    cat(crayon::cyan('Expected file: ', barcodes.file, '\n'))
  }
  count.file <- paste0(directory, '/', matrix.file)
  genes.file <- paste0(directory, '/', genes.file)
  barcodes.files <- paste0(directory, '/', barcodes.file)
  genes_ensembl <- utils::read.table(file = genes.file, sep = '\t', header = FALSE)
  barcodes <- utils::read.table(file = barcodes.files, sep = '\t', header = FALSE)
  mm <- Matrix::readMM(count.file)
  cat(crayon::cyan('Files loaded\n'))
  true.length <- length(genes_ensembl$V2)
  if(true.length != length(unique(genes_ensembl$V2))) {
    if(isTRUE(make.feat.unique)) {
      cat(crayon::cyan('Non-unique features identified\n'))
      genes_ensembl$V2 <- make.unique(genes_ensembl$V2)
    }
    
  }
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  cat(crayon::cyan('Success: Matrix concatenated\n'))
  mm <- as.matrix(mm)
  return(mm)
}

marrow_E <- Read10X_output(directory = '/Users/knight05/Raw_Data/Database_samples/healthy_references/BMMC_atlas/marrow_Ck')

# Metadata generator

library(Matrix)

cell_metadata <- function(assay, 
                          col.prefix) {
  total.counts <- colSums(assay)
  total.features <- colSums(assay != 0)
  df <- as.data.frame(as.numeric(total.counts))
  df[['total.features']] <- as.numeric(total.features)
  colnames(df) <- c(paste0(col.prefix, '_total.counts'), 
                    paste0(col.prefix, '_total.features'))
  return(df)
}

feature_metadata <- function(assay, col.prefix) {
  df <- as.data.frame(as.numeric(rowSums(assay)))
  rownames(df) <- rownames(assay)
  df$temp <- rowSums(assay > 0)
  colnames(df) <- c(paste0(col.prefix,'_total.counts'), paste0(col.prefix,'_total.cells'))
  return(df)
}

subset_IBRAP <- function(object, 
                         features = NULL, 
                         cells = NULL) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be class IBRAP'))
    return(NULL)
    
  }
  
  if(is.null(features)) {
    
    if(is.numeric(cells)) {
      
      features <- length(rownames(object@methods[[object@active.method]]@counts))
      
    } else if(is.character(cells)) {
      
      features <- rownames(object@methods[[object@active.method]]@counts)
      
    } else (
      
      cat(crayon::cyan('cells must be a character string or numeric'))
      
    )
    
  }
  
  if(is.null(cells)) {
    
    if(is.numeric(features)) {
      
      cells <- length(colnames(object@methods[[object@active.method]]@counts))
      
    } else if(is.character(features)) {
      
      cells <- colnames(object@methods[[object@active.method]]@counts)
      
    } else (
      
      cat(crayon::cyan('features must be a character string or numeric'))
      
    )
    
  }
  
  .sample_metadata <- object@sample_metadata[cells,]
  
  list.methods <- list()
  
  for(x in names(object@methods)) {
    
    sub <- object@methods[[x]]
    
    if(length(as.matrix(sub@counts)) != 0){
      
      .counts <- sub@counts[features,cells]
      
    } else {
      
      .counts <- sub@counts
      
    }
    
    if(length(as.matrix(sub@normalised)) != 0) {
      
      .normalised <- sub@normalised[features,]
      
    } else {
      
      .normalised <- sub@normalised
      
    }
    
    if(length(as.matrix(sub@norm.scaled)) != 0) {
      
      .norm.scaled <- sub@norm.scaled[features,cells]
      
    } else {
      
      .norm.scaled <- sub@norm.scaled
      
    }
    
    if(length(as.matrix(sub@feature_metadata)) != 0) {
      
      .feature_metadata <- sub@feature_metadata[features,]
      
    } else {
      
      .feature_metadata <- sub@feature_metadata
      
    }
    
    if(length(names(sub@computational_reductions)) != 0) {
      
      .computational_reductions <- list()
      
      for(comp in names(sub@computational_reductions))
        
        .computational_reductions[[comp]] <- sub@computational_reductions[[comp]][cells,]
      
    } else {
      
      .computational_reductions <- sub@computational_reductions
      
    }
    
    if(length(names(sub@visualisation_reductions)) != 0) {
      
      .visualisation_reductions <- list()
      
      for(vis in names(sub@visualisation_reductions)) {
        
        .visualisation_reductions[[vis]] <- sub@visualisation_reductions[[vis]][cells,] 
        
      }
      
    } else {
      
      .visualisation_reductions <- sub@visualisation_reductions
      
    }
    
    if(length(names(sub@cluster_assignments)) != 0) {
      
      .cluster_assignments <- sub@cluster_assignments[cells,]
      
    } else {
      
      .cluster_assignments <- sub@cluster_assignments
      
    }
    
    list.methods[[x]] <- new(Class = 'method', 
                             counts = .counts, 
                             normalised = .normalised, 
                             norm.scaled = .norm.scaled,
                             highly.variable.genes = sub@highly.variable.genes,
                             feature_metadata = .feature_metadata,
                             graphs = sub@graphs,
                             computational_reductions = .computational_reductions,
                             visualisation_reductions = .visualisation_reductions,
                             cluster_assignments = .cluster_assignments,
                             benchmark_results = sub@benchmark_results,
                             alt_objects = sub@alt_objects,
                             method.name = sub@method.name)
    
  }
  
  ibrap <- new(Class = 'IBRAP', 
               methods = list.methods, 
               sample_metadata = .sample_metadata,
               pipelines = object@pipelines, 
               active.variable = object@active.variable,
               doublet.barcodes = object@doublet.barcodes,
               active.method = object@active.method
  )
  
  print(ibrap)
  return(ibrap)
}

filter_cell_IBRAP <- function(object, assay, slot, min.features) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be class IBRAP'))
    return(NULL)
    
  }
  
  counts <- object@methods[[assay]][[slot]]
  
  nfeatures <- Matrix::colSums(x = counts > 0)
  counts <- counts[, which(x = nfeatures >= min.features)]
  
  object <- object[,colnames(counts)]
  
  return(object)
  
}

filter_feature_IBRAP <- function(object, assay, slot, min.cells,
                                 ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be class IBRAP'))
    return(NULL)
    
  }
  
  counts <- object@methods[[assay]][[slot]]
  
  num.cells <- Matrix::rowSums(x = counts > 0)
  counts <- counts[which(x = num.cells >= min.cells), ]
  
  object <- object[rownames(counts),]
  
  return(object)
  
}


isUnique <- function(vector){
  return(!any(duplicated(vector)))
}

# define classes 

setClass(Class = 'IBRAP', 
         representation = representation(
           methods = 'list', 
           sample_metadata = 'data.frame',
           pipelines = 'data.frame', 
           active.variable = 'character',
           doublet.barcodes = 'character',
           active.method = 'character'
         ))

setMethod(f = 'show', signature = 'IBRAP', definition = function(object) {
  cat(crayon::white(paste0('An object of class ', 
                           class(object), 
                           '\n')))
  
  cat(crayon::white(paste0('  ', 
                           nrow(object@methods[[object@active.method]]@counts), 
                           ' features by ', 
                           ncol(object@methods[[object@active.method]]@counts), 
                           ' samples\n')))
  
  cat(crayon::white(paste0('  ', 'Active method:', 
                           object@active.method, ' (features:', nrow(object@methods[[object@active.method]]@counts), 
                           ', samples:', paste0(ncol(object@methods[[object@active.method]]@counts),')'), '\n')))
  
  lol <- names(object@methods)[1]
  
  if(length(names(object@methods)) > 1) {
    
    for(x in names(object@methods)[2:length(names(object@methods))]) {
      
      lol <- paste0(lol, ', ', x)
      
    }
    
  }
  
  cat(crayon::white(paste0('  Available methods: ', lol, '\n')))
  
})

setClass(Class = 'methods',
         representation = representation(
           counts = 'dgCMatrix', 
           decontaminated = 'dgCMatrix',
           normalised = 'dgCMatrix', 
           norm.scaled = 'matrix',
           highly.variable.genes = 'character',
           feature_metadata = 'data.frame',
           graphs = 'list',
           computational_reductions = 'list',
           integration_reductions = 'list',
           visualisation_reductions = 'list',
           cluster_assignments = 'list',
           benchmark_results = 'list',
           alt_objects = 'list',
           method.name = 'character'
         ))

### add decontaminated slot

setMethod(f = 'show', signature = 'methods', definition = function(object) {
  cat(crayon::white(paste0(object@method.name, ', an IBRAP method\n')))
  cat(crayon::white(paste0('  Contains: ', 
                           nrow(object@counts), 
                           ' features by ', 
                           ncol(object@counts), 
                           ' samples\n')))
})

setMethod(f = 'rownames', signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            rownames(x@methods[[x@active.method]]@counts)
          })

setMethod(f = 'colnames', signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@methods[[x@active.method]]@counts)
          })

setMethod(f = 'dim', signature = 'IBRAP',
          function(x) {
            dim(x@methods[[x@active.method]]@counts)
          })

setMethod(f = '[[', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   ...) {
            x@sample_metadata[[i, exact = TRUE]]
          })

setMethod(f = '[[<-', signature = 'IBRAP', 
          function(x, 
                   i, 
                   j, 
                   value) {
            x@sample_metadata[[i]] <- value
            return(x)
          })

setMethod(f = 'as.list', signature = 'methods',
          function(x) {
            
            new.list <- list(counts = x@counts, 
                             decontaminated = x@decontaminated,
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
                             alt_objects = x@alt_objects,
                             method.name = x@method.name)
            return(new.list)
            
          })

setMethod(f = '[[', signature = 'methods',
          function(x, 
                   i, 
                   j, 
                   ...) {
            
            y <- as.list(x)
            y[[i]]
            
          })

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
            if(length(as.matrix(y$decontaminated))) {
              
              x@decontaminated <- y@decontaminated
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@highly.variable.genes <- y$highly.varaible.genes
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@feature_metadata <- y$feature_metadata
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@graphs <- y$graphs
              
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
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@alt_objects <- y$alt_objects
              
            }
            if(length(as.matrix(y$norm.scaled)) != 0) {
              
              x@method.name <- y$method.name
              
            }
            
            return(x)
            
          })

setMethod(f = '$', signature = 'IBRAP',
          function(x, 
                   name){
            
            x@sample_metadata[[name]]
            
          })

.convert_subset_index <- function(x, 
                                  y) {
  
  return(y[x])
  
}

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
                
                if(length(as.matrix(x@methods[[p]]@decontaminated)) != 0) {
                  
                  .decontaminated <- x@methods[[p]]@decontaminated[ii, , drop = FALSE]
                  
                } else {
                  
                  .decontaminated <- x@methods[[p]]@decontaminated
                  
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
                
                .method.name <- x@methods[[p]]@method.name
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         decontaminated = .decontaminated,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         graphs = .graphs,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         benchmark_results = .benchmark_results,
                                         alt_objects = .alt_objects,
                                         method.name = .method.name)
                
              }
              
              .pipelines <- x@pipelines
              
              .active.variable <- x@active.variable
              
              .doublet.barcodes <- x@doublet.barcodes
              
              .active.method <- x@active.method
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata,
                         pipelines = .pipelines, 
                         active.variable = .active.variable,
                         doublet.barcodes = .doublet.barcodes,
                         active.method = .active.method))
              
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
                
                if(length(as.matrix(x@methods[[p]]@decontaminated))) {
                  
                  .decontaminated <- x@methods[[p]]@decontaminated[ , jj, drop = FALSE]
                  
                } else {
                  
                  .decontaminated <- x@methods[[p]]@decontaminated
                  
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
                
                .method.name <- x@methods[[p]]@method.name
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         decontaminated = .decontaminated,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         graphs = .graphs,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         benchmark_results = .benchmark_results,
                                         alt_objects = .alt_objects,
                                         method.name = .method.name)
                
              }
              
              .pipelines <- x@pipelines
              
              .active.variable <- x@active.variable
              
              .doublet.barcodes <- x@doublet.barcodes
              
              .active.method <- x@active.method
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata,
                         pipelines = .pipelines, 
                         active.variable = .active.variable,
                         doublet.barcodes = .doublet.barcodes,
                         active.method = .active.method))
              
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
                
                if(length(as.matrix(x@methods[[p]]@decontaminated))) {
                  
                  .decontaminated <- x@methods[[p]]@decontaminated[ii, jj, drop = FALSE]
                  
                } else {
                  
                  .decontaminated <- x@methods[[p]]@decontaminated
                  
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
                
                .method.name <- x@methods[[p]]@method.name
                
                list.methods[[p]] <- new(Class = 'methods', 
                                         counts = .counts,
                                         decontaminated = .decontaminated,
                                         normalised = .normalised,
                                         norm.scaled = .norm.scaled,
                                         highly.variable.genes = .highly.variable.genes,
                                         feature_metadata = .feature_metadata,
                                         graphs = .graphs,
                                         computational_reductions = .computational_reductions,
                                         integration_reductions = .integration_reductions,
                                         visualisation_reductions = .visualisation_reductions,
                                         benchmark_results = .benchmark_results,
                                         alt_objects = .alt_objects,
                                         method.name = .method.name)
                
              }
              
              .pipelines <- x@pipelines
              
              .active.variable <- x@active.variable
              
              .doublet.barcodes <- x@doublet.barcodes
              
              .active.method <- x@active.method
              
              return(new(Class = 'IBRAP', 
                         methods = list.methods, 
                         sample_metadata = .sample_metadata,
                         pipelines = .pipelines, 
                         active.variable = .active.variable,
                         doublet.barcodes = .doublet.barcodes,
                         active.method = .active.method))
            }
          })

setMethod(f = 'merge', signature = 'IBRAP',
          function(x, 
                   y){
            
            items <- c(x,y)
            
            for(i in items) {
              
              if(length(i@methods[[1]]) > 1) {
                
                cat(crayon::cyan('No analysis can be performed prior to merging'))
                return(NULL)
                
              }
              
            }
            
            column.names <- list()
            counts.list <- list()
            sample.list <- list()
            
            for(i in 1:length(items)) {
              
              column.names[[i]] <- colnames(items[[i]]@methods[[1]]@counts)
              counts.list[[i]] <- as.matrix(items[[i]]@methods[[1]]@counts)
              sample.list[[i]] <- as.data.frame(items[[i]]@sample_metadata)
              
            }
            
            column.names <- unlist(column.names)
            column.names[!isUnique(column.names)] <- make.unique(column.names[!isUnique(column.names)])
            
            for(q in 1:length(items)) {
              
              colnames(counts.list[[q]]) <- column.names[1:length(colnames(counts.list[[q]]))]
              rownames(sample.list[[q]]) <- column.names[1:length(rownames(sample.list[[q]]))]
              column.names <- column.names[sum(length(colnames(counts.list[[q]]))+1):length(column.names)]
              
            }
            
            .counts <- counts.list[[1]]
            .sample_metadata <- sample.list[[1]]
            .feature_metadata <- as.data.frame(items[[1]]@methods[[1]]@feature_metadata)
            
            for(t in 2:length(items)) {
              
              .counts <- merge(x = .counts, y = counts.list[[t]], by = 'row.names', all = T)
              rownames(.counts) <- .counts$Row.names
              .counts$Row.names <- NULL
              
              .sample_metadata <- rbind(.sample_metadata, sample.list[[t]])
              
              .feature_metadata <- merge(.feature_metadata, 
                                         as.data.frame(items[[t]]@methods[[1]]@feature_metadata), 
                                         by = 'row.names', all = T)
              .feature_metadata[is.na(.feature_metadata)] <- 0
              rownames(.feature_metadata) <- .feature_metadata$Row.names
              .feature_metadata$Row.names <- NULL
              .feature_metadata$total.cells.x <- .feature_metadata$total.cells.x + .feature_metadata$total.cells.y
              .feature_metadata$total.counts.x <- .feature_metadata$total.counts.x + .feature_metadata$total.counts.y
              .feature_metadata <- .feature_metadata[,1:2]
              colnames(.feature_metadata) <- c('total.cells', 'total.counts')
              
            }
            
            .counts[is.na(.counts)] <- 0
            
            .counts <- Matrix::Matrix(data = as.matrix(.counts), sparse = T)
            .sample_metadata[match(colnames(.counts), rownames(.sample_metadata)),]
            .feature_metadata[match(rownames(.counts), rownames(.feature_metadata)),]
            
            new.method <- list()
            
            new.method[[names(x@methods)[1]]] <- new(Class = 'methods',
                                                     counts = .counts,
                                                     feature_metadata = .feature_metadata,
                                                     method.name = names(x@methods)[1]
            )
            
            ibrap <- new(Class = 'IBRAP',
                         methods = new.method, 
                         sample_metadata = .sample_metadata,
                         active.method = names(x@methods)[1])
            
            return(ibrap)
            
          })

pancreas.data <- readRDS(file = "~/Raw_Data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS('~/Raw_Data/pancreas_v3_files/pancreas_metadata.rds')
pancreas.data <- as.matrix(pancreas.data)
pancreas.data <- round(pancreas.data)

### removing scrublets

perform.scrublet <- function(counts,
                             total_counts = NULL, 
                             sim_doublet_ratio = 2.0, 
                             n_neighbors = NULL, 
                             expected_doublet_rate = 0.075, 
                             stdev_doublet_rate = 0.02, 
                             random_state = 0L,
                             synthetic_doublet_umi_subsampling = 1.0,
                             use_approx_neighbors = TRUE,
                             distance_metric = 'euclidean',
                             get_doublet_neighbor_parents = FALSE,
                             min_counts = 3L,
                             min_cells = 3L,
                             min_gene_variability_pctl = 85L,
                             log_transform = FALSE,
                             mean_center = TRUE, 
                             normalize_variance = TRUE,
                             n_prin_comps = 30L,
                             svd_solver = 'arpack') {
  
  cat(crayon::cyan('Initialising scrublet\n'))
  scrublet <- reticulate::import('scrublet', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  
  scrub1 <- scrublet$Scrublet(counts_matrix = t(as.data.frame(as.matrix(counts))), 
                              total_counts = total_counts, 
                              sim_doublet_ratio = sim_doublet_ratio, 
                              n_neighbors = n_neighbors, 
                              expected_doublet_rate = expected_doublet_rate, 
                              stdev_doublet_rate = stdev_doublet_rate, 
                              random_state = random_state)
  
  cat(crayon::cyan('scrublet object created\n'))
  
  res1 <- reticulate::py_to_r(scrub1$scrub_doublets(synthetic_doublet_umi_subsampling = synthetic_doublet_umi_subsampling,
                                                    use_approx_neighbors = use_approx_neighbors, 
                                                    distance_metric = distance_metric, 
                                                    get_doublet_neighbor_parents = get_doublet_neighbor_parents, 
                                                    min_counts = min_counts,
                                                    min_cells = min_cells, 
                                                    min_gene_variability_pctl = min_gene_variability_pctl,
                                                    log_transform = log_transform,
                                                    mean_center = mean_center,
                                                    normalize_variance = normalize_variance,
                                                    n_prin_comps = n_prin_comps,
                                                    svd_solver = svd_solver,
                                                    verbose = TRUE))
  
  sim.plot <- ggplot2::qplot(as.vector(reticulate::py_to_r(scrub1$doublet_scores_sim_)), 
                             geom = 'histogram') + 
    ggplot2::stat_bin(bins = 100) + 
    ggplot2::xlab('doublet scores') + 
    ggplot2::ylab('frequency') + 
    ggplot2::ggtitle(paste0('simulated_doublets')) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  obs.plot <- ggplot2::qplot(as.vector(res1)[[1]], 
                             geom = 'histogram') + 
    ggplot2::stat_bin(bins = 80) + 
    ggplot2::xlab('doublet scores') + 
    ggplot2::ylab('frequency') + 
    ggplot2::ggtitle(paste0('observed doublets')) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  comb.plot <- cowplot::plot_grid(sim.plot, obs.plot, ncol = 2, nrow = 1)
  print(comb.plot)
  
  cat(crayon::cyan('doublets detected\n'))
  counts <- counts[,!res1[[2]] == TRUE]
  cat(crayon::cyan('matrix scrubbed\n'))
  
  return(counts)
}

pancreas.data <- perform.scrublet(counts = pancreas.data)

# decontamination

perform.decontX <- function(counts,
                            z = NULL,
                            batch = NULL,
                            maxIter = 500,
                            delta = c(10, 10),
                            estimateDelta = TRUE,
                            convergence = 0.001,
                            iterLogLik = 10,
                            varGenes = 5000,
                            dbscanEps = 1,
                            seed = 12345,
                            logfile = NULL,
                            verbose = TRUE) {
  
  if(is.null(batch)) {
    
    d <- celda::decontX(x = counts,
                        z = z,
                        batch = NULL,
                        maxIter = maxIter,
                        delta = delta,
                        estimateDelta = estimateDelta,
                        convergence = convergence,
                        iterLogLik = iterLogLik,
                        varGenes = varGenes,
                        dbscanEps = dbscanEps,
                        seed = seed,
                        logfile = logfile,
                        verbose = verbose)
    
  } else {
    
    d <- celda::decontX(x = counts,
                        z = z,
                        batch = object@sample_metadata[,batch],
                        maxIter = maxIter,
                        delta = delta,
                        estimateDelta = estimateDelta,
                        convergence = convergence,
                        iterLogLik = iterLogLik,
                        varGenes = varGenes,
                        dbscanEps = dbscanEps,
                        seed = seed,
                        logfile = logfile,
                        verbose = verbose)
    
  }
  
  cat(crayon::cyan('Decontamination comlpleted\n'))
  
  cat(crayon::cyan(paste0(formatC(sum(d$contamination)/length(d$contamination), 
                                  digits = 2), 
                          '% average contamination\n')))
  
  clean.matrix <- d$decontXcounts
  cat(crayon::cyan('Matrix isolated\n'))
  clean.matrix <- round(clean.matrix)
  zero.samples <- colSums(clean.matrix) > 0
  clean.matrix <- clean.matrix[,zero.samples]
  cat(crayon::cyan('converted to integer\n'))
  return(clean.matrix)
  
}

pancreas.data <- perform.decontX(counts = pancreas.data)

# IBRAP object

createIBRAPobject <- function(counts, 
                              original.project, 
                              method.name = 'RAW', 
                              meta.data = NULL,
                              min.cells=NULL,
                              min.features=NULL) {
  
  if(!is.character(original.project)) {
    
    cat(crayon::cyan('original.project must be a character string\n'))
    return(NULL)
    
  }
  
  if(!is.character(method.name)) {
    
    cat(crayon::cyan('method.name must be a character string\n'))
    return(NULL)
    
  }
  
  if(class(counts)[1] != 'dgCMatrix') {
    
    cat(crayon::cyan('Converting counts input into sparse matrix\n'))
    counts <- Matrix::Matrix(counts, sparse = T)

  }
  
  cat(crayon::cyan(paste0('Adding ', original.project, ' as barcode prefix\n')))
  
  if(!is.null(meta.data)) {
    
    if(!rownames(meta.data) %in% colnames(counts)) {
      
      cat(crayon::cyan('meta.data rownames must be the same as counts colnames'))
      
      return(NULL)
      
    } else {
      
      rownames(meta.data) <- paste0(original.project, '_', rownames(meta.data))
      
    }
    
  }
  
  colnames(counts) <- paste0(original.project, '_', colnames(counts))
  
  if(!is.null(min.features)) {
    
    nfeatures <- Matrix::colSums(x = counts > 0)
    
    counts <- counts[, which(x = nfeatures >= min.features)]
    
  }
  
  if(!is.null(min.cells)) { 
    
    num.cells <- Matrix::rowSums(x = counts > 0)
    
    counts <- counts[which(x = num.cells >= min.cells), ]
    
  }
  
  meta <- as.data.frame(replicate(n = length(colnames(counts)), expr = original.project))
  
  colnames(meta) <- 'original.project'
  
  meta.2 <- cell_metadata(assay = counts, col.prefix = method.name)
  
  meta <- cbind(meta, meta.2)
  
  rownames(meta) <- colnames(counts)
  
  f.metadata <- feature_metadata(assay = counts, col.prefix = method.name)
  
  if(!is.null(meta.data)) {
    
    cat(crayon::cyan('Concatenating metadata\n'))
    
    l1 <- colnames(meta)
    
    l2 <- colnames(meta.data)
    
    if(isFALSE(isUnique(c(l1,l2)))) {
      
      cat(crayon::cyan('Column names from meta.data cannot be named:', 'original.project, counts_total.counts or counts_total.features\n'))
      
      return(NULL)
      
    }
    
    meta <- meta[match(rownames(meta.data), rownames(meta)),]
    
    meta <- cbind(meta, meta.data)
    
    meta <- meta[match(colnames(counts), rownames(meta)),]
    
  }
  
  f.metadata <- f.metadata[match(rownames(counts), rownames(f.metadata)),]
  
  ##########################################################
  
  first.method <- new('methods', 
                      counts = counts,
                      feature_metadata = f.metadata,
                      method.name = method.name)
  
  methods <- list()
  
  methods[[as.character(method.name)]] <- first.method
  
  IBRAP.obj <- new(Class = 'IBRAP', 
                   methods = methods, 
                   sample_metadata = meta,
                   active.method = as.character(method.name))
  
  return(IBRAP.obj)
  
}

# Load practice data

# data.kumar <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/sce_full/sce_full_Kumar.rds')
# mat.kumar <- assay(data.kumar,'counts')
# mat.kumar <- round(mat.kumar)
# extmeta.kumar <- as.data.frame(data.kumar$phenoid)
# colnames(extmeta.kumar) <- 'celltype'

# kumar <- createIBRAPobject(counts = mat.kumar, 
#                            original.project = 'kumar', 
#                            method.name = 'RAW', 
#                            meta.data = extmeta.kumar, 
#                            min.cells = 3)
# 
# data.obj <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/sce_full/sce_full_Koh.rds')
# mat <- assay(data.obj,'counts')
# mat <- round(mat)
# extmeta <- as.data.frame(data.obj$phenoid)
# colnames(extmeta) <- 'celltype'
# 
# koh <- createIBRAPobject(counts = mat, 
#                            original.project = 'koh', 
#                            method.name = 'RAW', 
#                            meta.data = extmeta, 
#                            min.cells = 3)

panc <- createIBRAPobject(counts = FL1,
                          original.project = 'marrow_Ck', 
                          method.name = 'RAW')

panc <- panc[,panc$tech == 'celseq2']

find_percentage_genes <- function(object, 
                                  pattern='^MT-', 
                                  assay='RAW', 
                                  slot='counts',
                                  column.name = 'RAW_percent.mt') {
  cat(crayon::cyan('Calculating percentage\n'))
  mat <- as.matrix(object@methods[[assay]][[slot]])
  temp <- colSums(mat[grep(pattern = pattern, x = rownames(mat)),]) / colSums(mat) * 100
  
  cat(crayon::cyan('Percentage calculated\n'))
  temp <- as.data.frame(temp)
  colnames(temp) <- column.name
  if(column.name %in% colnames(object@sample_metadata)) {
    cat(crayon::cyan('Removing old metadata column\n'))
    object@sample_metadata <- object@sample_metadata[,colnames(object@sample_metadata) != column.name]
  }
  cat(crayon::cyan('Appending new column\n'))
  colnames(temp) <- column.name
  object@sample_metadata <- cbind(object@sample_metadata, temp)
  return(object)
}

panc <- find_percentage_genes(object = panc, pattern = '^MT-', 
                              assay = 'RAW', slot = 'counts',
                              column.name = 'RAW_percent.mt')
panc <- find_percentage_genes(object = panc, pattern = 'RPL', 
                              assay = 'RAW', slot = 'counts',
                              column.name = 'RAW_percent.rp')

require(egg)

plot.QC.vln <- function(object, 
                        metadata.columns=c('RAW_total.features', 'RAW_total.counts'), 
                        split.by='original.project') {
  
  plots.list <- list()
  metadata <- object@sample_metadata
  
  for(m in metadata.columns) {
    
    if(!m %in% colnames(metadata)) {
      
      cat(crayon::cyan('Provided column names do not exist'))
      return(NULL)
      
    }
    
  }
  
  for(o in metadata.columns) {
    
    new.metadata <- data.frame(project=as.factor(object[[split.by]]))
    new.metadata$sample <- as.factor(colnames(object))
    new.metadata$variable <- object[[o]]
    plots.list[[o]] <- ggplot2::ggplot(data = new.metadata, 
                                       mapping = ggplot2::aes(x=variable, y=project, fill=project)) + 
      ggplot2::geom_violin() + ggplot2::coord_flip() + ggplot2::ggtitle(o) + 
      ggplot2::xlab('') + ggplot2::ylab('project') + ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(face = 'bold', angle = 45, vjust = 1, hjust=1), 
                     legend.position="none", plot.title = ggplot2::element_text(hjust=0.5))
  }
  
  do.call(what = 'ggarrange', args = list(plots = plots.list, nrow=1, ncol=length(plots.list)))
  
}

plot.QC.vln(object = panc)

plot.QC.scatter <- function(object, 
                            x, 
                            y, 
                            split.by) {
  
  metadata <- object@sample_metadata
  
  if(!x %in% colnames(metadata)) {
    
    cat(crayon::cyan('X variable does not exist'))
    
  }
  
  if(!y %in% colnames(metadata)) {
    
    cat(crayon::cyan('Y variable does not exist'))
    
  }
  
  if(!is.null(split.by)) {
    
    if(!split.by %in% colnames(metadata)){
      
      if(!split.by %in% colnames(metadata)) {
        
        cat(crayon::cyan('split.by variable does not exist'))
        
      }
      
    }
    
  }
  
  new.df <- data.frame(as.factor(rownames(metadata)))
  new.df$x <- metadata[,x]
  new.df$y <- metadata[,y]
  new.df$project <- metadata[,split.by]
  
  p <- ggplot2::ggplot(data = new.df, mapping = ggplot2::aes(x = x, y = y, col = project)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + ggplot2::ggtitle(paste0(x,'_vs_',y)) + 
    ggplot2::ylab(y) + ggplot2::xlab(x) + labs(color='identifier') 
  
  print(p)
  
}

plot.QC.scatter(object = panc, x = 'RAW_total.counts', y = 'RAW_total.features', split.by = 'original.project')

panc <- perform.decontX(object = panc)


panc <- filter_cell_IBRAP(object = panc, assay = 'RAW', slot = 'counts', min.features = 200)
panc <- filter_feature_IBRAP(object = panc, assay = 'RAW', slot = 'counts', min.cells = 3)

cell_metadata()

panc <- find_percentage_genes(object = panc, pattern = '^MT-', 
                              assay = 'RAW', slot = 'decontaminated',
                              column.name = 'decontaminated_percent.mt')
panc <- find_percentage_genes(object = panc, pattern = 'RPL', 
                              assay = 'RAW', slot = 'decontaminated',
                              column.name = 'decontaminated_percent.rp')

sd.value <- sd(panc$decontaminated_total.features)
med.value <- median(panc$decontaminated_total.features)
max.features <- (sd.value*3)+med.value

panc <- filter_cell_IBRAP(object = panc, decontaminated_total.features < max.features)

add.cell.cycle <- function(object, 
                           assay,
                           slot,
                           transform, ...) {
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/data/Homo_sapiens.csv', header = TRUE, sep = ',')
  cat(crayon::cyan('Cell cycle genes loaded\n'))
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]])
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Cell cycle scores identified\n'))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    object@sample_metadata <- cbind(object@sample_metadata, seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][['counts']])
    seuobj@assays$RNA@data <- object@methods[[assay]][[slot]]
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Data transformed\n'))
    
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    
    object@sample_metadata <- cbind(object@sample_metadata, seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
    
  }
  return(object)
}

panc <- add.cell.cycle(object = panc, assay = 'RAW', slot = 'counts', transform = TRUE)

add.feature.score <- function(object, 
                              assay, 
                              slot,
                              transform, 
                              features, 
                              column.name,
                              ...) {
  genes <- rownames(object)
  genes <- list(genes[genes %in% features])
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]])
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = genes, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    object@sample_metadata[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][['counts']])
    seuobj@assays$RNA@data <- object@methods[[assay]][[slot]]
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    for(o in names(seuobj@meta.data)) {
      
      if(o %in% names(object@sample_metadata)) {
        
        cat(crayon::cyan(paste0('found duplicated column name: ',o, 'removing old column names.\n')))
        object@sample_metadata[,o] <- NULL
        
      }
      
    }
    object@sample_metadata[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}

panc <- add.feature.score(object = panc, 
                          assay = 'RAW', 
                          slot = 'decontaminated',
                          transform = TRUE, 
                          features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
                                       'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
                                       'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
                                       'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
                                       'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
                          column.name = 'StressScore')

perform.sct.normalisation <- function(object, 
                                      assay,
                                      slot,
                                      new.assay.name = 'SCT',
                                      prefix = 'sctransform',
                                      do.scale = TRUE,
                                      do.center = TRUE,
                                      min_cells = 3,
                                      save.seuratobject = TRUE,
                                      ...) {
  
  cat(crayon::cyan('Converting to Seurat object\n'))
  seuratobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]], project = 'NA')
  cat(crayon::cyan('Initiating SCTransform\n'))
  seuratobj <- Seurat::SCTransform(object = seuratobj, do.scale = do.scale, do.center = do.center, min_cells = min_cells, ...)
  cat(crayon::cyan('SCTransform completed!\n'))
  .highly.variable.genes <- as.character(seuratobj@assays$SCT@var.features)
  .counts <- as(object = as.matrix(seuratobj@assays$SCT@counts), Class = 'dgCMatrix')
  .normalised <- as(as.matrix(seuratobj@assays$SCT@data), Class = 'dgCMatrix')
  .norm.scaled <- as.matrix(seuratobj@assays$SCT@scale.data)
  feat.meta <- feature_metadata(assay = .counts, col.prefix = new.assay.name)
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = .counts, 
                                          normalised = .normalised, 
                                          norm.scaled = .norm.scaled,
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feat.meta,
                                          method.name = as.character(new.assay.name))
  if(isTRUE(save.seuratobject)) {
    
    object@methods[[new.assay.name]]@alt_objects[['seurat']] <- seuratobj
    
  }
  cat(crayon::cyan('Populated IBRAP object\n'))
  return(object)
}

panc <- perform.sct.normalisation(object = panc, assay = 'RAW', slot = 'counts')

perform.scran.normalisation <- function(object, 
                                        assay = 'RAW',
                                        slot = 'decontaminated',
                                        batch=NULL,
                                        vars.to.regress=NULL,
                                        do.scale=TRUE,
                                        do.center=TRUE,
                                        new.assay.name = 'SCRAN',
                                        n.genes=1500,
                                        max.cluster.size = 1000,
                                        center_size_factors=TRUE) {
  
  assay.list <- list()
  mat <- object@methods[[assay]][[slot]]
  assay.list[[slot]] <- mat
  sce <- SingleCellExperiment::SingleCellExperiment(assay.list)
  
  clusters <- scran::quickCluster(mat)
  
  cat(crayon::cyan('quickCluster completed\n'))
  sce <- scran::computeSumFactors(sce, clusters=clusters, max.cluster.size=max.cluster.size, assay.type=slot)
  sce <- scuttle::logNormCounts(x = sce, log = F, center.size.factors=center_size_factors, exprs_values=slot)
  assay(sce, 'non-logged') <- assay(sce, 'normcounts')
  sce <- scuttle::logNormCounts(x = sce, log = T, center.size.factors=center_size_factors, exprs_values=slot)
  .counts <- assay(sce, 'non-logged')
  .normalised <- assay(sce, 'normcounts')
  cat(crayon::cyan('normalisation completed\n'))
  feat.meta <- feature_metadata(assay = .counts, col.prefix = new.assay.name)
  if(!is.null(batch)) {
    
    dec <- scran::modelGeneVar(sce, assay.type='normcounts', block=object@sample_metadata[[batch]])
    
  } else {
    
    dec <- scran::modelGeneVar(sce, assay.type='normcounts')
    
  }
  
  top.hvgs <- scran::getTopHVGs(stats = dec, n = n.genes)
  
  cat(crayon::cyan('HVGs identified\n'))
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  
  if(!is.null(vars.to.regress)) {
    
    vars.to.regress.df <- as.data.frame(object@sample_metadata[,vars.to.regress])
    colnames(vars.to.regress.df) <- vars.to.regress
    rownames(vars.to.regress.df) <- colnames(object)
    
    vars.to.regress.df <- vars.to.regress.df[match(rownames(seuobj@meta.data), rownames(vars.to.regress.df)),]
    seuobj@meta.data <- cbind(seuobj@meta.data,vars.to.regress.df)
    
    colnames(seuobj@meta.data) <- c(names(seuobj@meta.data)[1:sum(length(names(seuobj@meta.data))-length(vars.to.regress))], vars.to.regress)
    
    seuobj@assays$RNA@data <- .normalised[top.hvgs,]
    seuobj <- Seurat::ScaleData(object = seuobj, vars.to.regress=vars.to.regress, do.scale=do.scale, do.center=do.center)
    
  } else {
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center)
    
  }
  
  .norm.scaled <- seuobj@assays$RNA@scale.data
  
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = as(.counts, 'dgCMatrix'), 
                                          normalised = as(.normalised, 'dgCMatrix'), 
                                          norm.scaled = as.matrix(.norm.scaled),
                                          highly.variable.genes = top.hvgs,
                                          feature_metadata = feat.meta,
                                          method.name = as.character(new.assay.name))
  cat(crayon::cyan('Done\n'))
  return(object)
}

library(SingleCellExperiment)
panc <- perform.scran.normalisation(object = panc, assay = 'RAW', slot = 'counts', vars.to.regress = 'RAW_total.counts')

perform.tpm.normalisation <- function(object, 
                                      assay, 
                                      slot,
                                      n.genes,
                                      do.scale,
                                      do.center,
                                      vars.to.regress,
                                      new.assay.name = 'TPM') {
  
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/data/mart_export.csv', header = TRUE, sep = ',')
  r$Gene.length <- r$Gene.end..bp. - r$Gene.start..bp.
  
  subset <- r[r$Gene.name %in% rownames(object),]
  
  cat(crayon::cyan('Matrix subsetted\n'))
  
  rownames(subset) <- make.unique(names = as.character(subset$Gene.name), '.')
  
  cat(crayon::cyan('Rownames added\n'))
  
  meta <- object@methods[[assay]]@feature_metadata[intersect(rownames((object@methods[[assay]]@feature_metadata)), rownames(subset)),]
  
  cat(crayon::cyan('Gene names interesected\n'))
  
  mat <- as.matrix(object@methods[[assay]][[slot]])
  
  mat <- mat[intersect(rownames(mat), rownames(subset)),]
  
  ordered <- subset[match(rownames(mat), rownames(subset)),]
  
  cat(crayon::cyan('Matrices ordered\n'))
  
  cat(crayon::cyan('Calculated counts/feature length\n'))
  
  calc <- sweep(mat, 1, as.numeric(ordered$Gene.length), `/`)
  
  scale.factor <- colSums(calc)/1000000
  
  .counts <- sweep(calc, 2, as.numeric(scale.factor), `/`)
  
  cat(crayon::cyan('Calculations completed\n'))
  
  cat(crayon::cyan('log2(x+1) transforming\n'))
  .normalised <- log2(.counts+1)
  cat(crayon::cyan('Transformation completed\n'))
  
  dec <- scran::modelGeneVar(x = .normalised)
  .highly.variable.genes <- scran::getTopHVGs(stats = dec, n=n.genes)
  seuobj <- Seurat::CreateSeuratObject(counts = mat)
  seuobj@assays$RNA@data <- .normalised[.highly.variable.genes,]
  
  if(!is.null(vars.to.regress)) {
    
    vars.to.regress.df <- as.data.frame(object@sample_metadata[,vars.to.regress])
    colnames(vars.to.regress.df) <- vars.to.regress
    rownames(vars.to.regress.df) <- colnames(object)
    
    vars.to.regress.df <- vars.to.regress.df[match(rownames(seuobj@meta.data), 
                                                   rownames(vars.to.regress.df)),]
    seuobj@meta.data <- cbind(seuobj@meta.data,vars.to.regress.df)
    
    colnames(seuobj@meta.data) <- c(names(seuobj@meta.data)[1:sum(length(names(seuobj@meta.data))-length(vars.to.regress))], vars.to.regress)
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center,vars.to.regress=vars.to.regress)
    
  } else {
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center)
    
  }
  
  .norm.scaled <- seuobj@assays$RNA@scale.data
  
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = as(.counts, 'dgCMatrix'), 
                                          normalised = as(.normalised, 'dgCMatrix'), 
                                          norm.scaled = as.matrix(.norm.scaled),
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feature_metadata(assay = .counts, col.prefix = new.assay.name),
                                          method.name = as.character(new.assay.name))
  
  cat(crayon::cyan('Completed!\n'))
  
  return(object)
}

panc <- perform.tpm.normalisation(object = panc, assay = 'RAW', slot = 'counts', n.genes = 1500, 
                                  do.scale = T, do.center = T, vars.to.regress = 'RAW_total.counts')

perform.scanpy.normalisation <- function(object, 
                                         assay='RAW', 
                                         slot='counts', 
                                         new.assay.name='SCANPY', 
                                         target_sum = 1e6, 
                                         exclude_highly_expressed = FALSE,  
                                         max_fraction = 0.05, 
                                         key_added = 'scanpy_norm_factor',
                                         
                                         n_top_genes = 1500, 
                                         max_mean = 6, 
                                         min_mean = 0.0125, 
                                         min_disp = 0.5, 
                                         span = 0.3, 
                                         n_bins = 20, 
                                         flavor = 'seurat', 
                                         batch_key = NULL,
                                         
                                         do.scale=TRUE,
                                         vars.to.regress=NULL, 
                                         n_jobs = NULL, 
                                         zero_center = TRUE, 
                                         max_value = NULL, 
                                         obsm = NULL,
                                         
                                         save.anndata = TRUE
                                         ) {
  sc <- reticulate::import('scanpy')
  scobj <- sc$AnnData(X = t(as.matrix(object@methods[[assay]][[slot]])))
  scobj$obs_names <- as.factor(colnames(object))
  scobj$var_names <- as.factor(rownames(object))
  
  if(length(names(object@sample_metadata)) >= 1) {
    scobj$obs <- object@sample_metadata
  }
  
  cat(crayon::cyan('Normalising counts\n'))
  
  if(!is.null(target_sum) & !is.null(key_added)) {
    sc$pp$normalize_total(adata = scobj, target_sum = as.integer(target_sum), 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction), key_added = as.character(key_added))
  } else if (!is.null(target_sum)) {
    sc$pp$normalize_total(adata = scobj, target_sum = as.integer(target_sum), 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction))
  } else if (!is.null(key_added)) {
    sc$pp$normalize_total(adata = scobj, 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction))
  } else {
    sc$pp$normalize_total(adata = scobj)
  }
  
  .counts <- t(scobj$X)
  rownames(.counts) <- rownames(object)
  colnames(.counts) <- colnames(object)
  
  feat.metadata <- feature_metadata(assay = .counts, col.prefix = 'SCANPY')
  
  cat(crayon::cyan('Logging data\n'))
  
  sc$pp$log1p(scobj)
  
  .normalised <- t(scobj$X)
  rownames(.normalised) <- rownames(object)
  colnames(.normalised) <- colnames(object)
  
  cat(crayon::cyan('Computing highly variable genes\n'))
  
  if (!is.null(n_top_genes) & !is.null(batch_key)) {
    
    sc$pp$highly_variable_genes(adata = scobj, 
                                n_top_genes = as.integer(n_top_genes), 
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor), 
                                batch_key = as.character(batch_key))
    
  } else if (!is.null(n_top_genes)) {
    
    sc$pp$highly_variable_genes(adata = scobj, 
                                n_top_genes = as.integer(n_top_genes), 
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor))
    
  } else if (!is.null(batch_key)) {
    
    sc$pp$highly_variable_genes(adata = scobj,
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor))
    
  } else {
    
    sc$pp$highly_variable_genes(adata = scobj, 
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor))
    
  }
  
  .highly.variable.genes <- rownames(object)[scobj$var[['highly_variable']]]
  
  scobj2 <- sc$AnnData(X = t(.normalised[.highly.variable.genes,]))
  
    if(length(names(object@sample_metadata)) >= 1) {
      scobj2$obs <- object@sample_metadata
    }
    
  if(!is.null(do.scale)) {
    
    if(!is.null(vars.to.regress)) {
      if(!is.null(n_jobs)) {
        sc$pp$regress_out(adata = scobj2, keys = vars.to.regress, n_jobs = as.integer(n_jobs))
      } else {
        sc$pp$regress_out(adata = scobj2, keys = vars.to.regress)
      }
    }

    if(!is.null(max_value) & !is.null(obsm)) {
      sc$pp$scale(scobj2, zero_center = as.logical(zero_center), max_value = as.integer(max_value), obsm = as.character(obsm))
    } else if(!is.null(max_value)) {
      sc$pp$scale(scobj2, zero_center = as.logical(zero_center), max_value = as.integer(max_value))
    } else if(!is.null(obsm)) {
      sc$pp$scale(scobj2, zero_center = as.logical(zero_center), obsm = as.character(obsm))
    } else {
      sc$pp$scale(scobj2)
    }
    
    .norm.scaled <- t(scobj2$X)
    colnames(.norm.scaled) <- colnames(object)
    rownames(.norm.scaled) <- rownames(object)[rownames(object) %in% .highly.variable.genes]
    
  } else {
    
    .norm.scaled <- .normalised[.highly.variable.genes,]
    colnames(.norm.scaled) <- colnames(object)
    rownames(.norm.scaled) <- rownames(object)[rownames(object) %in% .highly.variable.genes]
    
  }
  
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = as(.counts, 'dgCMatrix'), 
                                          normalised = as(.normalised, 'dgCMatrix'), 
                                          norm.scaled = as.matrix(.norm.scaled),
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feat.metadata,
                                          method.name = as.character(new.assay.name))
  
  if(isTRUE(save.anndata)) {
    
    object@methods[[new.assay.name]]@alt_objects$anndata
    
  }
  
  
  
  return(object)
}

panc <- perform.scanpy.normalisation(object = panc, vars.to.regress = 'RAW_total.counts')

perform.bbknn <- function(object, 
                          reduction,
                          graph.name = 'bbknn',
                          column.correct,
                          n_pcs = NULL, 
                          trim = NULL, 
                          n_trees = 10,
                          use_faiss = TRUE,
                          set_op_mix_ratio = 1.0,
                          local_connectivity= 1) {
  
  sc <- reticulate::import('scanpy')
  print('.')
  scobj <- sc$AnnData(X = reducedDim(object, reduction))
  print('.')
  scobj$obs_names <- as.factor(colnames(object))
  print('.')
  scobj$var_names <- as.factor(colnames(reducedDim(object, reduction)))
  print('.')
  scobj$obsm$update(X_pca = reducedDim(object, reduction))
  print(scobj)
  
  if(as.data.frame(colData(object)) >= 1) {
    print('obs')
    print(as.data.frame(colData(object)))
    scobj$obs <- pd$DataFrame(data = as.data.frame(colData(object)))
  }
  
  if(is.null(n_pcs)) {
    print('npcs calculated')
    n_pcs <- as.integer(length(colnames(reducedDim(object, reduction))))
  }
  
  if(is.null(trim)) {
    print('initialising bbknn')
    sc$external$pp$bbknn(scobj, 
                         batch_key = as.character(column.correct), 
                         approx = as.logical(FALSE), 
                         n_pcs = n_pcs,
                         n_trees = as.integer(n_trees),
                         use_faiss = as.logical(use_faiss),
                         set_op_mix_ratio = set_op_mix_ratio,
                         local_connectivity = local_connectivity)
  } else if (!is.null(trim)) {
    sc$external$pp$bbknn(scobj, 
                         batch_key= as.character(column.correct), 
                         approx = as.logical(FALSE), 
                         n_pcs = n_pcs,
                         trim = as.integer(trim),
                         n_trees = as.integer(n_trees),
                         use_faiss = as.logical(use_faiss),
                         set_op_mix_ratio = set_op_mix_ratio,
                         local_connectivity = local_connectivity)
  }
  
  graph.list <- list()
  connectivities <- scobj$obsp[['connectivities']]
  colnames(connectivities) <- colnames(object)
  rownames(connectivities) <- colnames(object)
  distances <- scobj$obsp[['distances']]
  colnames(distances) <- colnames(object)
  rownames(distances) <- colnames(object)
  
  graph.list[['connectivities']] <- connectivities
  graph.list[['distances']] <- distances
  print(graph.list)
  metadata(object)[['graphs']][[graph.name]] <- graph.list
  
  return(object)
  
}

test <- perform.bbknn(object = pancreas.test, 
                      reduction = 'uncorrected_dbmap', 
                      graph.name = 'dbmap_bbknn', 
                      column.correct = 'tech')

perform.pca <- function(object, 
                        assay,
                        slot='norm.scaled',
                        n.pcs,
                        reduction.save='pca', 
                        ...) {
  for(t in assay) {
    
    mat <- object@methods[[t]][[slot]]
    cat(crayon::cyan('Initialising PCA for assay:', t, '\n'))
    a <- PCAtools::pca(mat = mat, center = F, scale = F)
    b <- PCAtools::findElbowPoint(a$variance)
    
    p <- PCAtools::screeplot(pcaobj = a, components = 1:sum(as.numeric(b)+10), 
                        title = paste0(assay,'_PCA_variance'), vline = b) +
      geom_label(aes(x = b, y = 50,
                     label = 'Elbow point', vjust = -1, size = 8)) +
      ggtitle(paste0(t,'_screeplot'))
    
    print(p)
    
    cat(crayon::cyan('PCA completed\n'))
    
    object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,n.pcs])
    
  }
  
  return(object)
  
}

panc <- perform.pca(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), n.pcs = 1:50, reduction.save = 'pca')

perform.dbmap <- function(object, 
                          assay, 
                          slot='normalised',
                          n_components = 100, 
                          n_neighbors = 15, 
                          reduction.save='dbmap',
                          save.object = TRUE) {
  
  scipy.sparse <- reticulate::import('scipy.sparse', convert = FALSE)
  
  dbmap <- reticulate::import('dbmap', convert = FALSE)
  
  for(o in assay) {
    
    cat(crayon::cyan(paste0('calculating dbmap for assay: ', o,'\n')))
    
    cellnames <- colnames(object@methods[[o]][[slot]])
    
    data <- scipy.sparse$csr_matrix(reticulate::r_to_py(t(object@methods[[o]][[slot]][object@methods[[o]]@highly.variable.genes,])))
    
    if(!is.null(n_components)) {
      
      diff <- dbmap$diffusion$Diffusor(n_components = as.integer(n_components), n_neighbors = as.integer(n_neighbors),
                                       transitions = as.logical(F),
                                       norm = as.logical(F), ann_dist = as.character('cosine'),
                                       n_jobs = as.integer(10), kernel_use = as.character('simple'))$fit(data)
      
    } else {
      
      diff <- dbmap$diffusion$Diffusor(n_neighbors = as.integer(n_neighbors),
                                       transitions = as.logical(F),
                                       norm = as.logical(F), ann_dist = as.character('cosine'),
                                       n_jobs = as.integer(10), kernel_use = as.character('simple'))$fit(data)
      
    }
    
    dbmap_components <- reticulate::py_to_r(diff$transform(data))
    
    res <- diff$return_dict()

    rownames(dbmap_components) <- cellnames
    
    dim.names <- list()
    for(t in 1:length(colnames(dbmap_components))) {
      dim.names[[t]] <- paste0('dbmap_', t)
    }
    colnames(dbmap_components) <- unlist(dim.names)
    
    object@methods[[o]]@computational_reductions[[reduction.save]] <- as.matrix(dbmap_components)
    
    if(isTRUE(save.object)) {
      
      object@methods[[o]]@alt_objects[[reduction.save]] <- diff
      
    }
    
  }
  
  return(object)
  
}

panc <- perform.dbmap(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'),  reduction.save = 'dbmap')

perform.umap <- function(object, 
                         assay,
                         reduction='pca',
                         reduction.save='umap',
                         n.dim=NULL, 
                         n_components = 3, 
                         ...) {
  
  for(u in assay) {
    
    cat(crayon::cyan('Processing', reduction, 'for assay:', u,'\n'))
    if(!is.null(n.dim)) {
      
      c <- uwot::umap(X = object@methods[[u]]@computational_reductions[[reduction]][,n.dim], n_components = n_components, verbose = TRUE, ...)
      
    } else {
      
      c <- uwot::umap(X = object@methods[[u]]@computational_reductions[[reduction]], n_components = n_components, verbose = TRUE, ...)
      
    }
    dim.names <- list()
    for(l in 1:n_components) {
      dim.names[[l]] <- paste0('umap_', l)
    }
    colnames(c) <- unlist(dim.names)
    object@methods[[u]]@visualisation_reductions[[reduction.save]] <- c
    
  }
  
  return(object)
}

panc <- perform.umap(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = 'pca', reduction.save = 'pca_umap', n.dim = 1:42, n_components = 3)
panc <- perform.umap(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = 'dbmap', reduction.save = 'dbmap_umap', n_components = 3)

perform.tsne <- function(object, 
                         assay,
                         reduction, 
                         reduction.save, 
                         n.dim=NULL,
                         n_components = 3, 
                         ...) {
  if(is(object = object, class2 = 'IBRAP') == FALSE) {
    cat(crayon::cyan('Must be IBRAP class object\n'))
  }
  if(is.null(assay)) {
    cat(crayon::cyan(paste0('Please provide assay\n')))
  }
  
  for(g in assay) {
    
    cat(crayon::cyan('t-SNE reduction initialising for assay:', g, '\n'))
    
    if(!is.null(n.dim)) {
      
      c <- ProjectionBasedClustering::tSNE(DataOrDistances = object@methods[[g]]@computational_reductions[[reduction]][,n.dim], 
                                           OutputDimension = n_components, Iterations = 1000, ...)$ProjectedPoints
      
    } else {
      
      c <- ProjectionBasedClustering::tSNE(DataOrDistances = object@methods[[g]]@computational_reductions[[reduction]], 
                                           OutputDimension = n_components, Iterations = 1000, ...)$ProjectedPoints
      
    }
    
    cat(crayon::cyan('t-SNE reduction completed\n'))
    dim.names <- list()
    for(t in 1:n_components) {
      dim.names[[t]] <- paste0('tsne_', t)
    }
    colnames(c) <- unlist(dim.names)
    rownames(c) <- colnames(object)
    object@methods[[g]]@visualisation_reductions[[reduction.save]] <- c
    
    cat(crayon::cyan('t-SNE data added\n'))
    
  }

  return(object)
}

panc <- perform.tsne(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = 'pca', reduction.save = 'pca_tsne', n.dim = 1:42, n_components = 3)
panc <- perform.tsne(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = 'dbmap', reduction.save = 'dbmap_tsne', n_components = 3)

library(harmony)

perform.harmony <- function(object, assay, group.by.vars, reduction = 'pca', dims.use = NULL,
                            theta = NULL, lambda = NULL, sigma = 0.1, nclust = NULL,
                            tau = 0, block.size = 0.05, max.iter.harmony = 10,
                            max.iter.cluster = 20, epsilon.cluster = 1e-05,
                            epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                            reference_values = NULL, reduction.save = "harmony", ...) {
  
  for(p in assay) {
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[p]]@computational_reductions), 
                   names(object@methods[[p]]@integration_reductions),
                   names(object@methods[[p]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[p]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        cat(crayon::cyan('reductions could not be found\n'))
        return(object)
        
      }
      
    }
    
    print('.')
    
    count <- 1
    
    for(g in reduction) {
      
      reduction.list <- reduction.list[[g]]
      
      red.save <- reduction.save[[count]]
      
      dims <- dims.use[[count]]
      
      if(is.null(dims)) {
        
        dims <- 1:length(colnames(reduction.list))
        
      }
      
      cat(crayon::cyan('Initialising harmony\n'))
      
      harm <- harmony::HarmonyMatrix(data_mat = reduction.list[,dims], meta_data = object@sample_metadata, vars_use = group.by.vars, do_pca = FALSE, 
                                     theta = theta, lambda = lambda, sigma = sigma, nclust = nclust, tau = tau, 
                                     block.size = block.size, max.iter.harmony = max.iter.harmony, max.iter.cluster = max.iter.cluster, 
                                     epsilon.cluster = epsilon.cluster, epsilon.harmony = epsilon.harmony, plot_convergence = plot_convergence, 
                                     return_object = FALSE, verbose = verbose, reference_values = reference_values)
      
      object@methods[[t]]@integration_reductions[[red.save]] <- harm
      
      count <- count + 1
      
    }
    
    
    
  }

  cat(crayon::cyan('Harmony completed\n'))
  return(object)
}

perform.harmony(object = panc, assay = c('SCRAN'), group.by.vars = 'celltype', reduction = c('pca', 'dbmap'))

perform.scanorama <- function(object, 
                              assay, 
                              split.by, 
                              reduced.dim, 
                              n.dims = 50, 
                              reduction.save='scanorama', 
                              assay.save='scanorama', 
                              batch_size = as.integer(5000), 
                              approx = TRUE, 
                              sigma = as.integer(15), 
                              alpha = as.integer(0.1), 
                              knn = as.integer(20), 
                              hvg = NULL) {
  cat(crayon::cyan('Initialising scanorama\n'))
  scanorama <- reticulate::import('scanorama', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  list.matrix <- list()
  column.names <- list()
  sep <- unique(colData(object)[,split.by])
  print('.')
  mat <- assay
  print('.')
  counter <- 1
  print('.')
  for(x in sep) {
    print('.')
    column.names[[counter]] <- colnames(mat[,colData(object)[,split.by] == x])
    print('.')
    list.matrix[[counter]] <- t(mat[,colData(object)[,split.by] == x])
    print('.')
    counter <- counter + 1
  }
  cat(crayon::cyan('Matrices isolated\n'))
  gene.list <- list()
  for(x in 1:length(sep)) {
    gene.list[[x]] <- rownames(mat[,colData(object)[,split.by] == x])
  }
  cat(crayon::cyan('Genes identified\n'))
  
  cat(crayon::cyan('Corrections starting\n'))
  integrated.corrected.data <- scanorama$correct(datasets_full = reticulate::r_to_py(list.matrix), 
                                                 genes_list = reticulate::r_to_py(gene.list), 
                                                 dimred = as.integer(n.dims), 
                                                 return_dimred=TRUE, 
                                                 return_dense=TRUE, 
                                                 verbose = TRUE, 
                                                 batch_size = as.integer(batch_size), 
                                                 approx = approx, 
                                                 sigma = as.integer(sigma), 
                                                 alpha = as.integer(alpha), 
                                                 knn = as.integer(knn), 
                                                 hvg = hvg)
  
  if(reduced.dim == TRUE) {
    dims <- list()
    cat(crayon::cyan('Isolating scanorama reduced dimensions\n'))
    dim.names <- list()
    for(c in 1:n.dims) {
      dim.names[[c]] <- paste0('scanorama_', c)
    }
    
    dim.names <- unlist(dim.names)
    
    for(x in 1:length(sep)) {
      transposed <- t(reticulate::py_to_r(integrated.corrected.data)[[1]][[x]])
      colnames(transposed) <- column.names[[x]]
      rownames(transposed) <- dim.names
      dims[[x]] <- transposed
    }
    cat(crayon::cyan('Combining samples\n'))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan('Samples concatenated\n'))
    reducedDim(object, reduction.save) <- t(combined)
    
    return(object)
  } else if (reduced.dim == FALSE) {
    dims <- list()
    cat(crayon::cyan('Isolating scanorama corrected gene matrix\n'))
    for(x in 1:length(sep)) {
      transposed <- t(reticulate::py_to_r(integrated.corrected.data)[[2]][[x]])
      print('.')
      colnames(transposed) <- column.names[[x]]
      print('.')
      rownames(transposed) <- reticulate::py_to_r(integrated.corrected.data)[[3]]
      print('.')
      dims[[x]] <- transposed
    }
    cat(crayon::cyan('Combining samples\n'))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan('Samples concatenated\n'))
    reducedDim(object, reduction.save) <- as.matrix(t(combined))
    cat(crayon::cyan('Scanorama completed\n'))
    return(object)
  }
  
}

perform.bbknn <- function(object, 
                          reduced.dim, 
                          dims, 
                          split.by, 
                          graph.save = 'bbknn') {
  anndata <- reticulate::import('anndata', convert = FALSE)
  sc <- reticulate::import('scanpy', convert = FALSE)
  bbknn <- reticulate::import('bbknn', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  pca <- reticulate::r_to_py(reducedDim(object, reduced.dim)[,dims])
  batch <- reticulate::r_to_py(object[[split.by]])
  adata <- anndata$AnnData(X = pca, obs = batch)
  cat(crayon::cyan('Anndata object created\n'))
  sc$tl$pca(adata)
  adata$obsm$X_pca <- pca
  cat(crayon::cyan('Correcting data\n'))
  bbknn$bbknn(adata = adata, batch_key = 0)
  cat(crayon::cyan('Data corrected\n'))
  
  mat  <- reticulate::py_to_r(adata$obsp[['connectivities']])
  rownames(mat) <- colnames(object)
  colnames(mat) <- colnames(object)
  
  metadata(object)[['graphs']][[graph.save]] <- as.Graph(as.matrix(mat))
  cat(crayon::cyan('BBKNN completed\n'))
  return(object)
}

perform.seurat.cluster <- function(object, 
                                   assay,
                                   reduction='pca', 
                                   assignment.df.name='seurat',
                                   res=c(0.1,0.2,0.3,0.4,0.5,
                                         0.6,0.7,0.8,0.9,1,
                                         1.1,1.2,1.3,1.4,1.5), 
                                   dims=NULL,
                                   prune.SNN=0, 
                                   nn.method='annoy', 
                                   annoy.metric='euclidean', 
                                   nn.eps=0.0, 
                                   ...) {
  
  for(p in assay) {
    
    cat(crayon::cyan(paste0('Calculating Seurat clusters for assay: ', p, '\n')))
    tmp <- Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts)
    tmp@reductions[[reduction]] <- Seurat::CreateDimReducObject(embeddings = object@methods[[p]]@computational_reductions[[reduction]], key = reduction)
    orig.name <- names(tmp@meta.data)
    if(!is.null(dims)) {
      
      tmp <- Seurat::FindNeighbors(object = tmp, reduction = reduction, verbose = TRUE, dims = dims, compute.SNN = TRUE, prune.SNN = prune.SNN,
                                   nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
      
    } else {
      
      tmp <- Seurat::FindNeighbors(object = tmp, reduction = reduction, verbose = TRUE, dims = 1:length(colnames(object@methods[[p]]@computational_reductions[[reduction]])), compute.SNN = TRUE, prune.SNN = prune.SNN,
                                   nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
      
    }
    
    cat(crayon::cyan('Neighbours identified\n'))
    tmp <- Seurat::FindClusters(object = tmp, resolution = res, ...)
    cat(crayon::cyan('Clusters identified\n'))
    post.name <- names(tmp@meta.data)
    
    post.name <- post.name[!post.name %in% orig.name]
    post.name <- post.name[1:length(post.name)-1]
    new.clusters <- tmp@meta.data[,post.name]
    new.clusters <- new.clusters[match(rownames(object@sample_metadata), rownames(new.clusters)),]
    object@methods[[p]]@cluster_assignments[[assignment.df.name]] <- new.clusters
    cat(crayon::cyan('Seurat clusters added\n'))
    
  }
  
  return(object)
  
}

panc <- perform.seurat.cluster(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = 'pca', 
                               assignment.df.name = 'pca_seurat', dims = 1:42)

panc <- perform.seurat.cluster(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = 'dbmap', 
                               assignment.df.name = 'dbmap_seurat')

library(SC3)

perform.sc3.reduction.cluster <- function(object, 
                                assay,
                                reduction,
                                dims,
                                assignment.df.name,
                                ks, 
                                n.core = 3) {
  
  cat(crayon::cyan('Initialising SC3 clustering\n'))
  
  for(p in assay) {
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[p]]@computational_reductions), 
                   names(object@methods[[p]]@integration_reductions),
                   names(object@methods[[p]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[p]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        cat(crayon::cyan('reductions could not be found\n'))
        return(object)
        
      }
      
    }
  
    reduction.list <- reduction.list[[reduction]]
    
    if(is.null(dims)) {
      
      dims <- 1:length(colnames(reduction.list))
      
    }
    
    print(t(reduction.list))
    
    temp.2 <- SingleCellExperiment(list('logcounts' = t(reduction.list)[dims,]))
    rowData(temp.2)$feature_symbol <- rownames(temp.2)
    temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
    temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
    temp.2 <- sc3_calc_dists(temp.2)
    temp.2 <- sc3_calc_transfs(temp.2)
    temp.2 <- sc3_kmeans(temp.2, ks = ks)
    temp.2 <- sc3_calc_consens(temp.2)
    object@methods[[p]]@cluster_assignments[[assignment.df.name]] <- as.data.frame(colData(temp.2))
    cat(crayon::cyan('SC3 clustering completed\n'))
    
  }

  return(object)
  
}

test <- perform.sc3.reduction.cluster(object = panc, 
                              assay = c('SCRAN'), 
                              reduction = 'pca', 
                              dims = 1:42, 
                              assignment.df.name = 'pca_SC3', 
                              ks = 10:16, 
                              n.core = 3)

perform.sc3.slot.cluster <- function(object, 
                                          assay,
                                          slot,
                                          HVGs,
                                          assignment.df.name,
                                          ks, 
                                          n.core = 3) {
  
  cat(crayon::cyan('Initialising SC3 clustering\n'))
  
  for(p in assay) {
    
    mat <- object@methods[[p]][[slot]]
    
    if(!is.null(HVGs)) {
      
      mat <- mat[object@methods[[p]]@highly.variable.genes,]
      
    }
    
    temp.2 <- SingleCellExperiment(list('logcounts' = mat))
    rowData(temp.2)$feature_symbol <- rownames(temp.2)
    temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
    temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
    temp.2 <- sc3_calc_dists(temp.2)
    temp.2 <- sc3_calc_transfs(temp.2)
    temp.2 <- sc3_kmeans(temp.2, ks = ks)
    temp.2 <- sc3_calc_consens(temp.2)
    object@methods[[p]]@cluster_assignments[[assignment.df.name]] <- as.data.frame(colData(temp.2))
    cat(crayon::cyan('SC3 clustering completed\n'))
    
  }
  
  return(object)
  
}

perform.tsne.kmeans <- function(object, 
                                reduction=NULL, 
                                k=NULL,
                                data_frame_name,
                                method='kmeans',
                                ...) {
  
  clusters <- data.frame(kmeans=numeric(length(colnames(object))))
  rownames(clusters) <- colnames(object)
  
  if(is.null(k)) {
    cat(crayon::cyan('Specify number of clusters\n'))
  }
  if(is.null(reduction)) {
    cat(crayon::cyan('Provide assay\n'))
  }
  if(method == 'pam') {
    for(i in k) {
      clusters[,paste0('pam_clustering_K_', i)] <- cluster::pam(x = reducedDim(object, reduction), k = i, ...)$clustering
    }
  }
  if(method == 'kmeans') {
    for(i in k) {
      clusters[[paste0('kmeans_clustering_K_', i)]] <- kmeans(x = reducedDim(object, reduction), centers = i, ...)$cluster
    }
  } else {
    cat(crayon::cyan('Please specify method: pam or kmeans\n'))
  }
  print(clusters)
  print(data_frame_name)
  clusters <- clusters[,2:length(colnames(clusters))]
  metadata(object)[['clustering']][[data_frame_name]] <- clusters
  return(object)
}

benchmark.clustering <- function(object, 
                                 assay,
                                 clustering,
                                 reduction, 
                                 components, 
                                 dist.method='euclidean',
                                 ground.truth=NULL) {
  
  for(l in assay) {
    
    print(l)
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[l]]@computational_reductions), 
                   names(object@methods[[l]]@integration_reductions),
                   names(object@methods[[l]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[l]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[l]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[l]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[l]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[l]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[l]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        cat(crayon::cyan('reductions could not be found\n'))
        return(object)
        
      }
      
    }
    
    count <- 1
    
    reduction.list <- reduction.list[reduction]
    reduction.list <- reduction.list[match(reduction, names(reduction.list))]
    
    for(k in clustering) {
      
      reduction_sub <- reduction.list[[reduction[count]]][,1:3]
      
      count <- count + 1
      
      clusters <- object@methods[[l]]@cluster_assignments[[k]]
      
      for(p in colnames(clusters)) {
        
        if(length(unique(clusters[,p])) <= 1) {
          
          cat(crayon::cyan(paste0('cluster column: ', p, ' contains only 1 cluster group, omitting now\n')))
          clusters[,p] <- NULL
          
        }
        
      }
      
      dist.matrix <- dist(x = reduction_sub, method = dist.method)
      
      sil.results <- data.frame(average_silhoeutte=NA)
      
      for (v in colnames(clusters)[1:length(colnames(clusters))]) {
       
        cat(crayon::cyan(paste0('Calculating silhouette for ', v, '\n')))
        tmp <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters[,v])), dist = dist.matrix)
        average <- sum(tmp[,3])/length(tmp[,3])
        sil.results[v,] <- average
        
      }
      
      sil.results <- sil.results[complete.cases(sil.results),]
      max.AS <- max(sil.results)
      print(max.AS)
      
      dunn.results <- data.frame(dunn.index=NA)
      for (p in colnames(clusters)[1:length(colnames(clusters))]){
        
        cat(crayon::cyan(paste0('Calculating dunn index for ', p, '\n')))
        dunn.results[p,] <- clValid::dunn(distance = dist.matrix, 
                                          clusters = as.numeric(x = as.factor(x = clusters[,p])))
        
      }
      
      dunn.results <- dunn.results[complete.cases(dunn.results),]
      max.dunn <- max(dunn.results)
      print(max.dunn)
      
      conn.results <- data.frame(connectivity=NA)
      for (p in colnames(clusters)[1:length(colnames(clusters))]){
        
        cat(crayon::cyan(paste0('Calculating connectivity for ', p, '\n')))
        conn.results[p,] <- clValid::connectivity(distance = dist.matrix, clusters = clusters[,p])
        
      }
      
      conn.results <- conn.results[complete.cases(conn.results),]
      min.conn <- min(conn.results)
      print(min.conn)
      
      if(!is.null(ground.truth)) {
        
        ARI.results <- data.frame(ARI=NA)
        for (p in colnames(clusters)[1:length(colnames(clusters))]){
          
          cat(crayon::cyan(paste0('Calculating ARI for ', p, '\n')))
          ARI.results[p,] <- mclust::adjustedRandIndex(x = clusters[,p], y = ground.truth)
          
        }
        
        ARI.results <- ARI.results[complete.cases(ARI.results),]
        max.ARI <- max(ARI.results)
        print(max.ARI)
        
        NMI.results <- data.frame(NMI=NA)
        for (p in colnames(clusters)[1:length(colnames(clusters))]) {
          
          cat(crayon::cyan(paste0('Calculating NMI for ', p, '\n')))
          NMI.results[p,] <- aricode::AMI(c1 = clusters[,p], c2 = ground.truth)
          
        }
        
        NMI.results <- NMI.results[complete.cases(NMI.results),]
        max.nmi <- max(NMI.results)
        print(max.nmi)
        
        results <- cbind(sil.results, dunn.results, conn.results, ARI.results, NMI.results)
        rownames(results) <- colnames(clusters)
        colnames(results) <- c(paste0(k, '_sil.results'), 
                               paste0(k, '_dunn.results'), 
                               paste0(k, '_conn.results'), 
                               paste0(k, '_ARI.results'), 
                               paste0(k, '_NMI.results'))
        
        object@methods[[l]]@benchmark_results[[k]] <- results
        
      } else {
        
        results <- cbind(sil.results, dunn.results, conn.results)
        rownames(results) <- colnames(all.clusters)
        colnames(results) <- c(paste0(k, '_sil.results'), paste0(k, '_dunn.results'), paste0(k, '_conn.results'))
        object@methods[[l]]@benchmark_results[[k]] <- results
        
      }
    }
  }
  
  return(object)
  
}

panc <- benchmark.clustering(object = panc, 
                     assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                     clustering = c('pca_seurat', 'dbmap_seurat'), 
                     reduction = c('pca_umap', 'dbmap_umap'), 
                     components = 1:3, 
                     dist.method = 'euclidean', 
                     ground.truth = panc$celltype)

library(plotly)
library(ggplot2)
library(egg)

panc <- readRDS("~/Results/scRNA-seq/IBRAP_development/Example_processed_results/RShiny_test.RData")

plot.reduced.dim <- function(object,
                             reduction,
                             assay,
                             clust.method,
                             column,
                             pt.size=5, 
                             metadata.access='clustering',
                             sub.access='metadata',
                             group.by, 
                             dimensions) {
  
  project.met <- object@sample_metadata
  assay.met <- object@methods[[assay]]@cluster_assignments[[clust.method]]
  assay.met <- assay.met[match(rownames(project.met), rownames(assay.met)),]
  meta <- cbind(project.met, assay.met)
  
  
  
  print('plot_cluster_dr_started')
  results <- as.data.frame(object@methods[[assay]]@visualisation_reductions[[reduction]])
  print('.')
  results[,'variable'] <- meta[,column]
  rownames(results) <- colnames(object)
  print('.')
  if(dimensions == 3) {
    print('3')
    print(plotly::plot_ly(data = results,
                          x = as.formula(paste0('~', colnames(results)[1])), 
                          y = as.formula(paste0('~', colnames(results)[2])),
                          z = as.formula(paste0('~', colnames(results)[3])), 
                          color = as.formula(paste0('~',colnames(results)[4])), 
                          colors = colorspace::qualitative_hcl(n = length(unique(results[,4])), palette = 'Dark 3'),
                          mode = "markers", 
                          marker = list(size = pt.size, width=0.5), 
                          text=as.formula(paste0('~',colnames(results)[4])), 
                          hoverinfo="text", plot_bgcolor = 'black'))
    
  } else if (dimensions == 2) {
    print('2')
    print(plotly::plot_ly(data = as.data.frame(results), 
                          x = as.formula(paste0('~', colnames(results)[1])), 
                          y = as.formula(paste0('~', colnames(results)[2])), 
                          color = as.formula(paste0('~',colnames(results)[4])),
                          colors = colorspace::qualitative_hcl(n = length(unique(results[,4])), palette = 'Dark 3'), 
                          mode = "markers", 
                          marker = list(size = pt.size, width=0.5), 
                          text=as.formula(paste0('~',colnames(results)[4])), 
                          hoverinfo="text", plot_bgcolor = 'black'))
  }
  
}

plot.reduced.dim(object = panc, reduction = 'dbmap_umap', assay = 'SCRAN', 
                 clust.method = 'dbmap_seurat', column = 'RNA_snn_res.1.5', dimensions = 2, pt.size = 2)

plot.features <- function(object, 
                          reduction='', 
                          pt.size=10, 
                          assay, 
                          feature,
                          dimensions) {
  results <- as.data.frame(reducedDim(object, reduction))[,1:3]
  print(dimensions)
  iso <- assay(object, assay)[feature,]
  print('.')
  results[,feature] <- iso
  print('.')
  colnames(results)[4] <- 'feature'
  print('.')
  if(dimensions == 3){
    print('3')
    p <- plotly::plot_ly(data = results[order(results$feature),], 
                         x = as.formula(paste0('~', colnames(results)[1])), 
                         y = as.formula(paste0('~', colnames(results)[2])),
                         z = as.formula(paste0('~', colnames(results)[3])), 
                         color = as.formula(paste0('~',colnames(results)[4])), 
                         mode = "markers", colors = RColorBrewer::brewer.pal(n = 9, name = 'Blues')[3:6],
                         marker = list(size = 5, width=5), 
                         text=as.formula(paste0('~',colnames(results)[4])), 
                         hoverinfo="text", plot_bgcolor = 'black')
    p <- p %>% layout(title=feature)
    print(p)
  } 
  if (dimensions == 2) {
    print('2')
    p <- plotly::plot_ly(data = results[order(results$feature),], 
                         x = as.formula(paste0('~', colnames(results)[1])), 
                         y = as.formula(paste0('~', colnames(results)[2])),
                         color = as.formula(paste0('~',colnames(results)[4])),  
                         mode = "markers",
                         colors = RColorBrewer::brewer.pal(n = 9, name = 'Blues')[3:6],
                         marker = list(size = 5, width=5), 
                         text=as.formula(paste0('~',colnames(results)[4])), 
                         hoverinfo="text", plot_bgcolor = 'black')
    p <- p %>% layout(title=feature)
    print(p)
  } 
}

#plot.features(object = panc, reduction = 'scanorama_reduced_umap', 
#              pt.size = 1, assay = 'decontXcounts', feature = "A1BG-AS1", 
#              dimensions = 2)

plot.features.multiple <- function(object, 
                                   assay, 
                                   reduction, 
                                   lab.key, 
                                   features) {
  
  plot.list <- list()
  
  for(x in features) {
    print('1')
    results <- as.data.frame(reducedDim(object, reduction))[,1:2]
    print('2')
    print(assay)
    iso <- assay(object, assay)[x,]
    print('3')
    colnames(results) <- c('red_1', 'red_2')
    print('4')
    results[,x] <- iso
    print('5')
    colnames(results)[3] <- 'feature'
    print('6')
    plot.list[[x]] <- ggplot(data = results[order(results$feature),], aes(x = red_1, y = red_2)) + 
      geom_point(aes(color=feature)) + 
      scale_color_gradient(low = '#FFFF00', high = '#FF0000') + 
      theme_bw() + labs(title=x, x=paste0(lab.key,'_1'), y=paste0(lab.key,'_2')) + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10))
    print('7')
  }
  
  is.even <- function(x) {
    if(as.integer(x %% 2) == 0) {
      print(TRUE)
    } else {
      print(FALSE)
    }
  }
  
  if(!is.even(length(plot.list))) {
    plot.list[length(plot.list)+1] <- plot.list[1] + geom_blank()
  }
  
  if(length(plot.list) > 1){
    print('multi')
    if(length(plot.list) <3) {
      print('1x2')
      do.call('ggarrange', c(plots = plot.list, ncol = 1, nrow = 2))
    } else if(length(plot.list) <5) {
      print('2x2')
      do.call('ggarrange', c(plots = plot.list, ncol = 2, nrow = 2))
    } else if(length(plot.list) <7) {
      print('3x2')
      do.call('ggarrange', c(plots = plot.list, ncol = 3, nrow = 2))
    }
  } else {
    print('single')
    plot.list[1]
  }
}

# plot.features.multiple(object = panc, assay = 'sctransform', 
#                        reduction = 'scanorama_reduced_umap', 
#                        lab.key = 'umap', features = c('GCG', 'MAFA', 'MAFB'))

plot.barplot <- function(object, 
                         x.value, 
                         y.value) {
  object[['var']] <- y.value 
  object[['group']] <- x.value
  print(length(unique(object[['var']])))
  p <- dittoSeq::dittoBarPlot(object = object, 
                              var = 'var', 
                              group.by = 'group') + 
    labs(title='') + 
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
  
  print(p)
}

plot.barplot(object = panc,
             y.value = metadata(panc)$clustering[['seurat:pca_harmony']][['Seurat_res_0.2']],
             x.value = panc$tech)

plot.benchmarking <- function(object, 
                              clust.method, 
                              ARI){
  print('.')
  clust.bench <- metadata(object)[['benchmarking_clustering']][[as.character(clust.method)]]
  print('.')
  clust.bench <- as.data.frame(clust.bench)
  print(clust.bench)
  clust.bench[,'cluster_index'] <- rownames(clust.bench)
  if(ARI == TRUE) {
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'ARI', 'NMI', 'cluster_index')
  } else {
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'cluster_index')
  }
  
  colnames(clust.bench) <- labels
  print(clust.bench)
  list.plot <- list()
  
  for(o in 1:sum(length(labels)-2)) {
    
    label <- labels[as.numeric(o)]
    print(label)
    fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(label), group = 1)) +
      geom_point() +
      geom_line() + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    list.plot[[as.numeric(o)]] <- fig
    
  }
  
  last.label <- labels[as.numeric(sum(length(labels)-1))]
  
  last.fig <- fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(last.label), group = 1)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  list.plot[[as.numeric(sum(length(labels)-1))]] <- last.fig
  print('/')
  do.call('ggarrange', c(plots = list.plot, ncol = 5))
}

plot.heatmap <- function(object, 
                         assay, 
                         features, 
                         ...) {
  
  ass <- object[features,]
  ass <- assay(ass, assay)
  cat(crayon::cyan('Isolated assay\n'))
  
  logged <- log(ass+1)
  means <- apply(X = logged, MARGIN = 1, FUN = mean)
  standev <- apply(X = logged, MARGIN = 1, FUN = sd)
  z_scores <- logged - means / standev 
  cat(crayon::cyan('z-scores calculated\n'))
  
  p <- gplots::heatmap.2(x = z_scores, 
                         key.title = 'colour key', 
                         key.xlab = 'z-score', 
                         key.ylab = '', xlab = 'cells', 
                         ylab = 'genes', trace = 'none', margins = c(7,7))
  p
}

# plot.heatmap(object = pancreas, assay = 'counts', features = features)

plot.vln <- function(object, 
                     assay,
                     features, 
                     group.by, 
                     title = NULL, 
                     xlab = 'group', 
                     ylab = 'expression') {
  
  if(!length(group.by) == length(colnames(object))) {
    cat(crayon::cyan('group.by variable must match the length of column names in SCE.'))
    return(NULL)
  }
  
  plot.list <- list()
  
  for(x in features) {
    ass <- t(assay(object[x,], assay))
    df <- data.frame(groups = colnames(object))
    rownames(df) <- df$barcodes
    df <- cbind(df, ass)
    df[,'groups'] <- group.by
    colnames(df) <- c('group', 'feature')
    
    p <- ggplot(data = df, aes(x = group, y = feature, color = group)) + 
      geom_violin() + 
      theme_light() + labs(title = x, y = 'expression', x = 'group') + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))
    
    plot.list[[x]] <- p
  }
  
  if(length(plot.list) > 1) {
    do.call('ggarrange', c(plots = plot.list, ncol = 1))
  } else {
    plot.list[1]
  }
}

# plot.vln(object = panc, 
#          assay = 'sctransform', 
#          features = c("A1BG-AS1", "A1BG", "A1CF"),
#          group.by = metadata(panc)[['clustering']][['seurat:scanorama_reduced']][['Seurat_res_0.2']])

perform.seurat.diffexp <- function(object, counts.assay, norm.assay, scaled.assay, test, cluster.method, cluster.column, ...) {
  seuobj <- Seurat::as.Seurat(x = object, counts = counts.assay, data = norm.assay)
  seuratobj$clusters <- metadata(pancreas)[['clustering']][['seurat:scanorama_reduced']][['Seurat_res_0.2']]
  Seurat::Idents(seuratobj) <- 'clusters'
  results <- Seurat::FindAllMarkers(object = seuobj, test.use = test)
  return(results)
}

tester <- perform.seurat.diffexp(object = panc, 
                                 counts.assay = 'sctransform_counts', 
                                 norm.assay = 'sctransform_data', 
                                 test = 'MAST', 
                                 cluster.method = 'seurat:scanorama_reduced', 
                                 cluster.column = 'Seurat_res_0.2')





