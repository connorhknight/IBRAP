data.obj <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/sce_full/sce_full_Koh.rds')
kumar.obj <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/sce_full/sce_full_Kumar.rds')
mat <- assay(data.obj,'counts')
mat <- round(mat)

reticulate::use_python('/Users/knight05/Library/r-miniconda/envs/r-reticulate/bin/python')
reticulate::import('scrublet', convert = FALSE)

# Read 10x

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

marrow_E <- Read10X_output(directory = '/Users/knight05/Raw_Data/Database_samples/healthy_references/BMMC_atlas/marrow_E')

# Metadata generator

cell_metadata <- function(assay, 
                          col.prefix) {
  total.counts <- Matrix::colSums(assay)
  total.features <- Matrix::colSums(assay != 0)
  df <- as.data.frame(as.numeric(total.counts))
  df[['total.features']] <- as.numeric(total.features)
  colnames(df) <- c(paste0(col.prefix, '_total.counts'), 
                    paste0(col.prefix, '_total.features'))
  return(df)
}

feature_metadata <- function(assay, 
                             col.prefix) {
  df <- as.data.frame(as.numeric(Matrix::rowSums(assay)))
  rownames(df) <- rownames(assay)
  df$temp <- as.numeric(Matrix::rowSums(assay > 0))
  colnames(df) <- c(paste0(col.prefix,'_total.counts'), paste0(col.prefix,'_total.cells'))
  return(df)
}

filter_IBRAP <- function(object, ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be class IBRAP'))
    return(NULL)
    
  }
  
  metadata <- object@sample_metadata
  
  metadata <- subset(metadata, ...)
  
  object <- object[,rownames(metadata)]
  
  return(object)
  
}

isUnique <- function(vector){
  return(!any(duplicated(vector)))
}

# define classes 

setClass(Class = 'IBRAP', 
         representation = representation(
           methods = 'list', 
           sample_metadata = 'data.frame'
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
           alt_objects = 'list'
         ))

### add decontaminated slot

setMethod(f = 'show', signature = 'methods', definition = function(object) {
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
            rownames(x@methods[[1]]@counts)
          })

setMethod(f = 'colnames', 
          signature = 'IBRAP', 
          function(x, 
                   do.NULL = TRUE, 
                   prefix = 'row') {
            colnames(x@methods[[1]]@counts)
          })

setMethod(f = 'dim', 
          signature = 'IBRAP',
          function(x) {
            dim(x@methods[[1]]@counts)
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

setMethod(f = 'merge', signature = 'IBRAP',
          function(x, 
                   y){
            
            items <- c(x,y)
            
            for(i in items) {
              
              if(length(i@methods[[1]]) > 1) {
                
                cat(crayon::cyan('No analysis can be performed prior to merging\n'))
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
              .feature_metadata$total.cells.x <- .feature_metadata$RAW_total.cells.x + .feature_metadata$RAW_total.cells.y
              .feature_metadata$total.counts.x <- .feature_metadata$RAW_total.counts.x + .feature_metadata$RAW_total.counts.y
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
                                                     feature_metadata = .feature_metadata)
            
            ibrap <- new(Class = 'IBRAP',
                         methods = new.method, 
                         sample_metadata = .sample_metadata)
            
            return(ibrap)
            
          })

pancreas.data <- readRDS(file = "~/Raw_Data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS('~/Raw_Data/pancreas_v3_files/pancreas_metadata.rds')
pancreas.data <- as.matrix(pancreas.data)
pancreas.data <- round(pancreas.data)

celseq2 <- pancreas.data[,rownames(metadata[metadata$tech == 'celseq2',])]
metadata_celseq2 <- metadata[rownames(metadata[metadata$tech == 'celseq2',]),]
celseq <- pancreas.data[,rownames(metadata[metadata$tech == 'celseq',])]
metadata_celseq <- metadata[rownames(metadata[metadata$tech == 'celseq',]),]

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
  
  if(!is(object = counts, class2 = 'matrix')) {
    
    if (!is(object = counts, class2 = 'dgCMatrix')) {
      
      cat(crayon::cyan('counts must be in matrix or dgCMatrix format\n'))
      return(counts)
      
    }

  } 
  
  if(!is.null(total_counts)) {
    
    if(!is.numeric(total_counts)) {
      
      cat(crayon::cyan('total_counts must be numerical\n'))
      return(counts)
      
    }
    
  }
  
  if(!is.numeric(sim_doublet_ratio)) {
    
    cat(crayon::cyan('sim_doublet_ratio must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.null(n_neighbors)) {
    
    if(!is.numeric(n_neighbors)) {
      
      cat(crayon::cyan('n_neighbors must be numerical\n'))
      return(counts)
      
    }
    
  }
  
  if(!is.numeric(expected_doublet_rate)) {
    
    cat(crayon::cyan('expected_doublet_rate must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(stdev_doublet_rate)) {
    
    cat(crayon::cyan('stdev_doublet_rate must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(random_state)) {
    
    cat(crayon::cyan('random_state must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(synthetic_doublet_umi_subsampling)) {
    
    cat(crayon::cyan('synthetic_doublet_umi_subsampling must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.logical(use_approx_neighbors)) {
    
    cat(crayon::cyan('use_approx_neighbors must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.character(distance_metric)) {
    
    cat(crayon::cyan('distance_metric must be character string\n'))
    return(counts)
    
  }
  
  if(!is.logical(get_doublet_neighbor_parents)) {
    
    cat(crayon::cyan('get_doublet_neighbor_parents must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.numeric(min_counts)) {
    
    cat(crayon::cyan('min_counts must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(min_cells)) {
    
    cat(crayon::cyan('min_cells must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(min_gene_variability_pctl)) {
    
    cat(crayon::cyan('min_gene_variability_pctl must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.logical(log_transform)) {
    
    cat(crayon::cyan('log_transform must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.logical(mean_center)) {
    
    cat(crayon::cyan('mean_center must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.logical(normalize_variance)) {
    
    cat(crayon::cyan('normalize_variance must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.numeric(n_prin_comps)) {
    
    cat(crayon::cyan('n_prin_comps must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.character(svd_solver)) {
    
    cat(crayon::cyan('svd_solver must be numerical\n'))
    return(counts)
    
  }
  
  cat(crayon::cyan('Initialising scrublet\n'))
  scrublet <- reticulate::import('scrublet', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  
  scrub1 <- scrublet$Scrublet(counts_matrix = as.data.frame(as.matrix(t(counts))), 
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
  counts <- as.matrix(counts)
  counts <- counts[,!res1[[2]]]
  counts <- Matrix::Matrix(data = counts, sparse = T)
  cat(crayon::cyan('matrix scrubbed\n'))
  
  return(counts)
}

celseq <- perform.scrublet(counts = celseq)
celseq2 <- perform.scrublet(counts = celseq2)

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
                            seed = 12345) {
  
  if(!is(object = counts, class2 = 'matrix')) {
    
    if (!is(object = counts, class2 = 'dgCMatrix')) {
      
      cat(crayon::cyan('counts must be in matrix or dgCMatrix format\n'))
      return(counts)
      
    }
    
  } 
  
  if(!is.null(z)) {
    
    if(length(z) != ncol(counts)) {
      
      cat(crayon::cyan('z must have the same length as ncol: counts\n'))
      return(counts)
      
    }
    
  }
  
  if(!is.numeric(maxIter)) {
    
    cat(crayon::cyan('maxIter must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(delta)) {
    
    cat(crayon::cyan('delta must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.logical(estimateDelta)) {
    
    cat(crayon::cyan('estimateDelta must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.numeric(convergence)) {
    
    cat(crayon::cyan('convergence must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(iterLogLik)) {
    
    cat(crayon::cyan('iterLogLik must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(varGenes)) {
    
    cat(crayon::cyan('varGenes must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(dbscanEps)) {
    
    cat(crayon::cyan('dbscanEps must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(seed)) {
    
    cat(crayon::cyan('seed must be numerical\n'))
    return(counts)
    
  }
  
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
                        verbose = TRUE)
    
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
                        verbose = TRUE)
    
  }
  
  cat(crayon::cyan('Decontamination comlpleted\n'))
  
  print(celda::plotDecontXContamination(x = d))
  
  cat(crayon::cyan(paste0(formatC(sum(d$contamination)/length(d$contamination), 
                                  digits = 2), 
                          '% average contamination\n')))
  
  clean.matrix <- d$decontXcounts
  cat(crayon::cyan('Matrix isolated\n'))
  clean.matrix <- round(clean.matrix)
  zero.samples <- Matrix::colSums(as.matrix(clean.matrix)) > 0
  clean.matrix <- clean.matrix[,zero.samples]
  cat(crayon::cyan('converted to integer\n'))
  return(clean.matrix)
  
}

celseq <- perform.decontX(counts = celseq)
celseq2 <- perform.decontX(counts = celseq2)

# IBRAP object

createIBRAPobject <- function(counts, 
                              original.project, 
                              method.name = 'RAW', 
                              meta.data = NULL,
                              min.cells=NULL,
                              min.features=NULL) {
  
  if(!is.character(original.project)) {
    
    cat(crayon::cyan('original.project must be a character string\n'))
    return(counts)
    
  }
  
  if(!is.character(method.name)) {
    
    cat(crayon::cyan('method.name must be a character string\n'))
    return(counts)
    
  }
  
  if(!is(object = counts, class2 = 'matrix')) {
    
    if (!is(object = counts, class2 = 'dgCMatrix')) {
      
      cat(crayon::cyan('counts must be in matrix or dgCMatrix format\n'))
      return(counts)
      
    }
    
  } 
  
  cat(crayon::cyan(paste0('Adding ', original.project, ' as barcode prefix\n')))

  counts <- as.matrix(x = counts)

  if(!is.null(meta.data)) {
    
    if(!is.data.frame(meta.data)) {
      
      cat(crayon::cyan('meta.data must be of class data.frame'))
      
    }
    
    if(FALSE %in% (rownames(meta.data) %in% colnames(counts))) {
      
      cat(crayon::cyan('meta.data rownames must be the same as counts colnames\n'))
      
      return(counts)
      
    } else {
      
      rownames(meta.data) <- paste0(original.project, '_', rownames(meta.data))
      
    }

  }
  
  colnames(counts) <- paste0(original.project, '_', colnames(counts))
  
  if(!is.null(min.features)) {
    
    if(!is.numeric(min.features)) {
      
      cat(crayon::cyan('min.features must be numerical'))
      return(counts)
      
    }
    
    nfeatures <- Matrix::colSums(x = counts > 0)
    
    counts <- counts[, which(x = nfeatures >= min.features)]
    
  }
  
  if(!is.null(min.cells)) { 
    
    if(!is.numeric(min.cells)) {
      
      cat(crayon::cyan('min.cells must be numerical'))
      return(counts)
      
    }
    
    num.cells <- Matrix::rowSums(x = counts > 0)
    
    counts <- counts[which(x = num.cells >= min.cells), ]
    
  }
  
  meta <- as.data.frame(replicate(n = length(colnames(counts)), expr = original.project))
  
  colnames(meta) <- as.character('original.project')

  meta.2 <- cell_metadata(assay = counts, col.prefix = method.name)
  
  for(f in colnames(meta.2)) {
    
    meta[,f] <- as.numeric(meta.2[,f])
    
  }
  
  rownames(meta) <- colnames(counts)
  
  f.metadata <- feature_metadata(assay = counts, col.prefix = method.name)
  
  if(!is.null(meta.data)) {
    
    cat(crayon::cyan('Concatenating metadata\n'))
    
    l1 <- colnames(meta)
    
    l2 <- colnames(meta.data)
    
    if(isFALSE(isUnique(c(l1,l2)))) {
      
      cat(crayon::cyan('Column names from meta.data cannot be named:', 'original.project, counts_total.counts or counts_total.features\n'))
      
      return(counts)
      
    }
    
    meta <- meta[match(rownames(meta.data), rownames(meta)),]
    
    meta <- cbind(meta, meta.data)
    
    meta <- meta[match(colnames(counts), rownames(meta)),]
    
  }
  
  f.metadata <- f.metadata[match(rownames(counts), rownames(f.metadata)),]
  
  ##########################################################

  first.method <- new('methods', 
                      counts = Matrix::Matrix(counts, sparse = T),
                      feature_metadata = f.metadata)

  methods <- list()
  
  methods[[as.character(method.name)]] <- first.method
  
  IBRAP.obj <- new(Class = 'IBRAP', 
                   methods = methods, 
                   sample_metadata = meta)
  
  return(IBRAP.obj)
  
}

metadata_celseq <- metadata_celseq[colnames(celseq),]

celseq <- createIBRAPobject(counts = celseq,
                              original.project = 'celseq',
                              method.name = 'RAW',
                              min.cells = 3,
                              min.features = 200)

metadata_celseq2 <- metadata_celseq2[colnames(celseq2),]

celseq2 <- createIBRAPobject(counts = celseq2@methods$RAW@counts,
                             original.project = 'pancreas_celseq2',
                             method.name = 'RAW',
                             min.cells = 3,
                             min.features = 200)

metadata_smartseq2 <- metadata[metadata$tech=='smartseq2',]
smartseq2 <- pancreas.data[,rownames(metadata_smartseq2)]

smartseq2 <- createIBRAPobject(counts = smartseq2,
                             original.project = 'pancreas_smartseq2', 
                             meta.data = metadata_smartseq2,
                             method.name = 'RAW',
                             min.cells = 3,
                             min.features = 200)

celseq_comb <- merge(x = celseq, y = celseq2)

find_percentage_genes <- function(object, 
                                  pattern='^MT-', 
                                  assay='RAW', 
                                  slot='counts',
                                  column.name = 'RAW_percent.mt') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(pattern)) {
    
    cat(crayon::cyan('pattern must be character string\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  cat(crayon::cyan('Calculating percentage\n'))
  mat <- as.matrix(object@methods[[assay]][[slot]])
  subbed <- mat[grep(pattern = pattern, x = rownames(mat)),]
  temp <- Matrix::colSums(subbed) / Matrix::colSums(mat) * 100
  
  cat(crayon::cyan('Percentage calculated\n'))
  temp <- as.data.frame(temp)
  colnames(temp) <- column.name
  if(column.name %in% colnames(object@sample_metadata)) {
    cat(crayon::cyan('Removing old metadata column\n'))
    object@sample_metadata <- object@sample_metadata[,colnames(object@sample_metadata) != column.name]
  }
  cat(crayon::cyan('Appending new column\n'))
  colnames(temp) <- column.name
  temp <- apply(temp, 2, function(x) as.numeric(x))
  object@sample_metadata <- cbind(object@sample_metadata, temp)
  return(object)
}

celseq_comb <- find_percentage_genes(object = celseq_comb, pattern = '^MT-',
                              assay = 'RAW', slot = 'counts',
                              column.name = 'RAW_percent.mt')
smartseq2 <- find_percentage_genes(object = smartseq2, pattern = 'RPL', 
                                 assay = 'RAW', slot = 'counts',
                                 column.name = 'RAW_percent.rp')
celseq2 <- find_percentage_genes(object = celseq2, pattern = 'RPL', 
                              assay = 'RAW', slot = 'counts',
                              column.name = 'RAW_percent.rp')
celseq2 <- find_percentage_genes(object = celseq2, pattern = '^MT-',
                                   assay = 'RAW', slot = 'counts',
                                   column.name = 'RAW_percent.mt')
celseq_comb <- find_percentage_genes(object = celseq_comb, pattern = '^MT-',
                                     assay = 'RAW', slot = 'counts',
                                     column.name = 'RAW_percent.mt')
celseq_comb <- find_percentage_genes(object = celseq_comb, pattern = 'RPL', 
                                     assay = 'RAW', slot = 'counts',
                                     column.name = 'RAW_percent.rp')

plot.QC.vln <- function(object, 
                        metadata.columns=c('RAW_total.features', 
                                           'RAW_total.counts'), 
                        split.by='original.project') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  plots.list <- list()
  metadata <- object@sample_metadata
  
  for(m in metadata.columns) {
    
    if(!m %in% colnames(metadata)) {
      
      cat(crayon::cyan('Provided column names do not exist\n'))
      return(NULL)
      
    }
    
  }
  
  if(!split.by %in% colnames(object@sample_metadata)) {
    
    cat(crayon::cyan(paste0(split.by, ' does not exist\n')))
    return(object)
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  cols <- RColorBrewer::brewer.pal(n = length(metadata.columns), name = 'Pastel2')
  
  count <- 1
  
  for(o in metadata.columns) {
    
    new.metadata <- data.frame(project=as.factor(object[[split.by]]))
    new.metadata$sample <- as.factor(colnames(object))
    new.metadata$variable <- object[[o]]
    proj.length <- length(unique(new.metadata$project))
    
    if(proj.length < 3) {
      
      proj.length.new <- 3
      cols.proj <- RColorBrewer::brewer.pal(n = proj.length.new, name = 'Pastel1')
      cols.proj <- cols.proj[1:proj.length]
      
    }
    
    plots.list[[o]] <- ggplot2::ggplot(data = new.metadata, 
                                       mapping = ggplot2::aes(x=variable, y=project, fill=project)) + 
      ggplot2::geom_violin() + ggplot2::coord_flip() + ggplot2::ggtitle(o) + 
      ggplot2::xlab('') + ggplot2::ylab('project') + ggplot2::theme_classic() + 
      ggplot2::geom_boxplot(lwd = 0.6, width = 0.09, fill = cols[[count]]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(face = 'bold', angle = 45, vjust = 1, hjust=1), 
                     legend.position="none", plot.title = ggplot2::element_text(hjust=0.5)) + 
      ggplot2::scale_fill_manual(values=cols.proj)
    count <- count + 1
  }
  
  do.call(what = 'ggarrange.tmp', args = list(plots = plots.list, nrow=1, ncol=length(plots.list)))
  
}

plot.QC.vln(object = celseq_comb, 
            metadata.columns = c('RAW_total.features', 
                                 'RAW_total.counts', 
                                 'RAW_percent.mt'))
plot.QC.vln(object = celseq2, 
            metadata.columns = c('RAW_total.features', 
                                 'RAW_total.counts', 
                                 'RAW_percent.rp'))
plot.QC.vln(object = celseq_comb, 
            metadata.columns = c('RAW_total.features', 
                                 'RAW_total.counts', 
                                 'RAW_percent.rp'))

plot.QC.scatter <- function(object, 
                            x, 
                            y, 
                            split.by) {
  
  metadata <- object@sample_metadata
  
  if(!x %in% colnames(metadata)) {
    
    cat(crayon::cyan('X variable does not exist\n'))
    
  }
  
  if(!y %in% colnames(metadata)) {
    
    cat(crayon::cyan('Y variable does not exist\n'))
    
  }
  
  if(!is.null(split.by)) {
    
    if(!split.by %in% colnames(metadata)){

      cat(crayon::cyan('split.by variable does not exist\n'))
      
    }
    
  }
  
  new.df <- data.frame(as.factor(rownames(metadata)))
  new.df$x <- metadata[,x]
  new.df$y <- metadata[,y]
  new.df$project <- metadata[,split.by]
  
  proj.length <- length(unique(new.df$project))
  
  if(proj.length < 3) {
    
    proj.length.new <- 3
    cols.proj <- RColorBrewer::brewer.pal(n = proj.length.new, name = 'Pastel1')
    cols.proj <- cols.proj[1:proj.length]
    
  }
  
  p <- ggplot2::ggplot(data = new.df, mapping = ggplot2::aes(x = x, y = y, col = project)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + ggplot2::ggtitle(paste0(x,'_vs_',y)) + 
    ggplot2::ylab(y) + ggplot2::xlab(x) + ggplot2::labs(color='identifier') + ggplot2::scale_color_manual(values=cols.proj)
  
  print(p)
  
}

plot.QC.scatter(object = celseq_comb, 
                x = 'RAW_total.counts', 
                y = 'RAW_total.features', 
                split.by = 'original.project')
plot.QC.scatter(object = celseq2, 
                x = 'RAW_total.counts', 
                y = 'RAW_total.features', 
                split.by = 'original.project')
plot.QC.scatter(object = celseq_comb, 
                x = 'RAW_total.counts', 
                y = 'RAW_total.features', 
                split.by = 'original.project')

sd.value <- sd(celseq_comb$RAW_total.features)
med.value <- median(celseq_comb$RAW_total.features)
max.features <- (sd.value*3)+med.value

celseq_comb <- filter_IBRAP(object = celseq_comb, 
                            RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 8)

sd.value <- sd(marrow_merged$RAW_total.features)
med.value <- median(marrow_merged$RAW_total.features)
max.features <- (sd.value*3)+med.value

marrow_merged <- filter_IBRAP(object = marrow_merged, 
                          RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 8)

sd.value <- sd(celseq2$RAW_total.features)
med.value <- median(celseq2$RAW_total.features)
max.features <- (sd.value*3)+med.value

celseq2 <- filter_IBRAP(object = celseq2, 
                        RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 8)

add.cell.cycle <- function(object, 
                           assay,
                           slot,
                           transform, ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot must be character string\n'))
    return(object)
    
  }
  
  if(!is.logical(transform)) {
    
    cat(crayon::cyan('transform must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  r <- utils::read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/data/Homo_sapiens.csv', header = TRUE, sep = ',')
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
    df <- seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))]
    object@sample_metadata <- cbind(object@sample_metadata, df)
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
    
    df <- seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))]
    object@sample_metadata <- cbind(object@sample_metadata, df)
    cat(crayon::cyan('New metadata added\n'))
    
  }
  return(object)
}

celseq_comb <- add.cell.cycle(object = celseq_comb, 
                              assay = 'RAW', 
                              slot = 'counts', 
                              transform = TRUE)

celseq2 <- add.cell.cycle(object = celseq2, 
                          assay = 'RAW', 
                          slot = 'counts', 
                          transform = TRUE)

celseq_comb <- add.cell.cycle(object = celseq_comb, 
                              assay = 'RAW', 
                              slot = 'counts', 
                              transform = TRUE)

add.feature.score <- function(object, 
                              assay, 
                              slot,
                              transform, 
                              features, 
                              column.name,
                              ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.logical(transform)) {
    
    cat(crayon::cyan('transform must be logical: TRUE/FALSEt\n'))
    return(object)
    
  }
  
  if(!is.character(features)) {
    
    cat(crayon::cyan('features must be character string(s)\n'))
    return(object)
    
  }
  
  if(!is.character(column.name)) {
    
    cat(crayon::cyan('column.name must be character string\n'))
    return(object)
    
  }
  
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

celseq_comb <- add.feature.score(object = celseq_comb, 
                               assay = 'RAW', 
                               slot = 'counts',
                               transform = TRUE, 
                               features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
                                            'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
                                            'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
                                            'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
                                            'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
                               column.name = 'StressScore')
celseq2 <- add.feature.score(object = celseq2, 
                             assay = 'RAW', 
                             slot = 'counts',
                             transform = TRUE, 
                             features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
                                          'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
                                          'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
                                          'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
                                          'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
                             column.name = 'StressScore')

celseq_comb <- add.feature.score(object = celseq_comb, 
                                 assay = 'RAW', 
                                 slot = 'counts',
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
                                      do.scale = TRUE,
                                      do.center = TRUE,
                                      n.genes = 1500,
                                      min_cells = 3,
                                      save.seuratobject = TRUE,
                                      ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be of class IBRAP\n'))
    return(object)
    
  }
  if(!is.character(assay)) {
    
    cat(crayon::cyan('Assay must be a character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('Slot must be a character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(new.assay.name)) {
    
    cat(crayon::cyan('new.assay.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.logical(do.scale)) {
    
    cat(crayon::cyan('do.scale must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.logical(do.center)) {
    
    cat(crayon::cyan('do.center must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.numeric(n.genes)) {
    
    cat(crayon::cyan('n.genes must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(min_cells)) {
    
    cat(crayon::cyan('min_cells must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(save.seuratobject)) {
    
    cat(crayon::cyan('save.seuratobject must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  cat(crayon::cyan('Converting to Seurat object\n'))
  seuratobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]][[slot]], project = 'NA')
  cat(crayon::cyan('Initiating SCTransform\n'))
  seuratobj <- Seurat::SCTransform(object = seuratobj, do.scale = do.scale, do.center = do.center, min_cells = min_cells, variable.features.n = n.genes, ...)
  cat(crayon::cyan('SCTransform completed!\n'))
  .highly.variable.genes <- as.character(seuratobj@assays$SCT@var.features)
  .counts <- as(object = as.matrix(seuratobj@assays$SCT@counts), Class = 'dgCMatrix')
  .normalised <- as(as.matrix(seuratobj@assays$SCT@data), Class = 'dgCMatrix')
  .norm.scaled <- as.matrix(seuratobj@assays$SCT@scale.data)
  feat.meta <- feature_metadata(assay = as.matrix(.counts), col.prefix = new.assay.name)
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = .counts, 
                                          normalised = .normalised, 
                                          norm.scaled = .norm.scaled,
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feat.meta)
  if(isTRUE(save.seuratobject)) {
    
    object@methods[[new.assay.name]]@alt_objects[['seurat']] <- seuratobj
    
  }
  cat(crayon::cyan('Populated IBRAP object\n'))
  return(object)
}

# celseq <- perform.sct.normalisation(object = celseq, 
#                                     assay = 'RAW', 
#                                     slot = 'counts')
# celseq2 <- perform.sct.normalisation(object = celseq2, 
#                                      assay = 'RAW', 
#                                      slot = 'counts')
celseq_comb <- perform.sct.normalisation(object = celseq_comb, 
                                       assay = 'RAW', 
                                       slot = 'counts')
celseq2 <- perform.sct.normalisation(object = celseq2, 
                                     assay = 'RAW', 
                                     slot = 'counts')
celseq_comb <- perform.sct.normalisation(object = celseq_comb, 
                                         assay = 'RAW', 
                                         slot = 'counts')

perform.scran.normalisation <- function(object, 
                                        assay = 'RAW',
                                        slot = 'counts',
                                        batch=NULL,
                                        vars.to.regress=NULL,
                                        do.scale=TRUE,
                                        do.center=TRUE,
                                        new.assay.name = 'SCRAN',
                                        n.genes=1500,
                                        max.cluster.size = 1000,
                                        center_size_factors=TRUE) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be of class IBRAP\n'))
    return(object)
    
  }
  if(!is.character(assay)) {
    
    cat(crayon::cyan('Assay must be a character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('Assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('Slot must be a character string\n'))
    return(object)
    
  }
  
  if(!is.null(batch)) {
    
    cat(crayon::cyan('batch must be a character string of column name contained in sample_metadata\n'))
    return(object)
    
  }
  
  if(!is.null(vars.to.regress)) {
    
    if(!is.character(vars.to.regress)) {
      
      cat(crayon::cyan('vars.to.regress must be a character string(s) of column name contained in sample_metadata\n'))
      return(object)
      
    }
    
  }
  
  if(!is.logical(do.scale)) {
    
    cat(crayon::cyan('do.scale must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.logical(do.center)) {
    
    cat(crayon::cyan('do.center must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.character(new.assay.name)) {
    
    cat(crayon::cyan('new.assay.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(n.genes)) {
    
    cat(crayon::cyan('n.genes must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(max.cluster.size)) {
    
    cat(crayon::cyan('max.cluster.size must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(center_size_factors)) {
    
    cat(crayon::cyan('center_size_factors must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  assay.list <- list()
  mat <- object@methods[[assay]][[slot]]
  assay.list[[slot]] <- mat
  sce <- SingleCellExperiment::SingleCellExperiment(assay.list)
  
  clusters <- scran::quickCluster(mat)
  
  cat(crayon::cyan('quickCluster completed\n'))
  sce <- scran::computeSumFactors(sce, clusters=clusters, max.cluster.size=max.cluster.size, assay.type=slot)
  sce <- scuttle::logNormCounts(x = sce, log = F, center.size.factors=center_size_factors, exprs_values=slot)
  SummarizedExperiment::assay(sce, 'non-logged') <- SummarizedExperiment::assay(sce, 'normcounts')
  sce <- scuttle::logNormCounts(x = sce, log = T, center.size.factors=center_size_factors, exprs_values=slot)
  .counts <- SummarizedExperiment::assay(sce, 'non-logged')
  .normalised <- SummarizedExperiment::assay(sce, 'normcounts')
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
                                          feature_metadata = feat.meta)
  cat(crayon::cyan('Done\n'))
  return(object)
}


celseq_comb <- perform.scran.normalisation(object = celseq_comb, 
                                         assay = 'RAW', 
                                         slot = 'counts', 
                                         vars.to.regress = 'RAW_total.counts')
celseq2 <- perform.scran.normalisation(object = celseq2, 
                                       assay = 'RAW', 
                                       slot = 'counts', 
                                       vars.to.regress = 'RAW_total.counts')
celseq_comb <- perform.scran.normalisation(object = celseq_comb, 
                                           assay = 'RAW', 
                                           slot = 'counts', 
                                           vars.to.regress = 'RAW_total.counts')

perform.tpm.normalisation <- function(object, 
                                      assay = 'RAW', 
                                      slot = 'counts',
                                      n.genes = 1500,
                                      do.scale = TRUE,
                                      do.center = TRUE,
                                      vars.to.regress = NULL,
                                      new.assay.name = 'TPM') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be a character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be a character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.numeric(n.genes)) {
    
    cat(crayon::cyan('n.genes must be numerical\n'))
    return(object)
    
  }
  
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/data/mart_export.csv', header = TRUE, sep = ',')
  
  if(is.null(r)) {
    
    cat(crayon::cyan('cannot find gene lengths\n'))
    
  }
  
  if(!is.logical(do.scale)) {
    
    cat(crayon::cyan('do.scale must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.logical(do.center)) {
    
    cat(crayon::cyan('do.center must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.null(vars.to.regress)) {
    
    if(!is.character(vars.to.regress)) {
      
      cat(crayon::cyan('vars.to.regress must be character string\n'))
      
    }
    
  }
  
  if(!is.character(new.assay.name)) {
    
    cat(crayon::cyan('new.assay.name must be a character string\n'))
    return(object)
    
  }
  
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
                                          counts = Matrix::Matrix(.counts, sparse = T), 
                                          normalised = Matrix::Matrix(.normalised, sparse = T), 
                                          norm.scaled = as.matrix(.norm.scaled),
                                          highly.variable.genes = .highly.variable.genes)
  
  cat(crayon::cyan('Completed!\n'))
  
  return(object)
}

celseq_comb <- perform.tpm.normalisation(object = celseq_comb, 
                                       assay = 'RAW', 
                                       slot = 'counts', 
                                       n.genes = 1500,
                                       vars.to.regress = 'RAW_total.counts')
celseq2 <- perform.tpm.normalisation(object = celseq2, 
                                     assay = 'RAW', 
                                     slot = 'counts', 
                                     n.genes = 1500,
                                     vars.to.regress = 'RAW_total.counts')
celseq_comb <- perform.tpm.normalisation(object = celseq_comb, 
                                         assay = 'RAW', 
                                         slot = 'counts', 
                                         n.genes = 1500,
                                         vars.to.regress = 'RAW_total.counts')

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

  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be of class IBRAP\n'))
    return(object)
    
  }
  if(!is.character(assay)) {
    
    cat(crayon::cyan('Assay must be a character string\n'))
    return(object)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('Slot must be a character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(new.assay.name)) {
    
    cat(crayon::cyan('new.assay.name must be character string \n'))
    return(object)
    
  }
  
  if(!is.numeric(target_sum)) {
    
    cat(crayon::cyan('target_sum must be numerical \n'))
    return(object)
    
  }
  
  if(!is.logical(exclude_highly_expressed)) {
    
    cat(crayon::cyan('exclude_highly_expressed must be logical\n'))
    return(object)
    
  }
  
  if(!is.numeric(max_fraction)) {
    
    cat(crayon::cyan('max_fraction must be numerical \n'))
    return(object)
    
  }
  
  if(!is.character(key_added)) {
    
    cat(crayon::cyan('key_added must be character string \n'))
    return(object)
    
  }
  
  if(!is.numeric(n_top_genes)) {
    
    cat(crayon::cyan('n_top_genes must be numerical \n'))
    return(object)
    
  }
  
  if(!is.numeric(max_mean)) {
    
    cat(crayon::cyan('max_mean must be numerical \n'))
    return(object)
    
  }
  
  if(!is.numeric(min_mean)) {
    
    cat(crayon::cyan('min_mean must be numerical \n'))
    return(object)
    
  }
  
  if(!is.numeric(min_disp)) {
    
    cat(crayon::cyan('min_disp must be numerical \n'))
    return(object)
    
  }
  
  if(!is.numeric(span)) {
    
    cat(crayon::cyan('span must be numerical \n'))
    return(object)
    
  }
  
  if(!is.character(flavor)) {
    
    cat(crayon::cyan('flavor must be character string \n'))
    return(object)
    
  }
  
  if(!is.null(batch_key)) {
    
    if(!is.character(batch_key)) {
      
      cat(crayon::cyan('batch_key must be character string\n'))
      return(object)
      
    }
    
  }
  
  if(!is.logical(do.scale)) {
    
    cat(crayon::cyan('do.scale must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.null(batch_key)) {
    
    if(!is.numeric(batch_key)) {
      
      cat(crayon::cyan('batch_key must be character string\n'))
      return(object)
      
    }
    
  }
  
  if(!is.null(vars.to.regress)) {
    
    if(!is.character(vars.to.regress)) {
      
      cat(crayon::cyan('vars.to.regress must be character string\n'))
      return(object)
      
    }
    
  }
  
  if(!is.logical(do.scale)) {
    
    cat(crayon::cyan('do.scale must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.null(n_jobs)) {
    
    if(!is.numeric(n_jobs)) {
      
      cat(crayon::cyan('n_jobs must be numerical\n'))
      return(object)
      
    }
    
  }
  
  if(!is.null(zero_center)) {
    
    if(!is.logical(zero_center)) {
      
      cat(crayon::cyan('zero_center must be logical: TRUE/FALSE\n'))
      return(object)
      
    }
    
  }
  
  if(!is.null(max_value)) {
    
    if(!is.numeric(max_value)) {
      
      cat(crayon::cyan('max_value must be numerical\n'))
      return(object)
      
    }
    
  }
  
  if(!is.null(obsm)) {
    
    if(!is.character(obsm)) {
      
      cat(crayon::cyan('obsm must be numerical\n'))
      return(object)
      
    }
    
  }
  
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
                                          feature_metadata = feat.metadata)
  
  if(isTRUE(save.anndata)) {
    
    object@methods[[new.assay.name]]@alt_objects$anndata
    
  }
  
  
  
  return(object)
}

celseq_comb <- perform.scanpy.normalisation(object = celseq_comb, 
                                          vars.to.regress = 'RAW_total.counts')
celseq2 <- perform.scanpy.normalisation(object = celseq2, 
                                        vars.to.regress = 'RAW_total.counts')
celseq_comb <- perform.scanpy.normalisation(object = celseq_comb, 
                                            vars.to.regress = 'RAW_total.counts')

perform.pca <- function(object, 
                        assay,
                        slot='norm.scaled',
                        n.pcs=50,
                        reduction.save='pca', 
                        ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('assay: ', x, ' does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be character string'))
    return(object)
    
  }
  
  if(!is.numeric(n.pcs)) {
    
    cat(crayon::cyan('n.pcs must be numerical'))
    return(object)
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be numerical'))
    return(object)
    
  }
  
  for(t in assay) {
    
    mat <- object@methods[[t]][[slot]]
    cat(crayon::cyan('Initialising PCA for assay:', t, '\n'))
    a <- PCAtools::pca(mat = mat, center = F, scale = F)
    b <- PCAtools::findElbowPoint(a$variance)
    
    p <- PCAtools::screeplot(pcaobj = a, components = 1:sum(as.numeric(b)+10), 
                        title = paste0(assay,'_PCA_variance'), vline = b) +
      ggplot2::geom_label(ggplot2::aes(x = b, y = 50,
                     label = 'Elbow point', vjust = -1, size = 8)) +
      ggplot2::ggtitle(paste0(t,'_screeplot'))
    
    print(p)
    
    cat(crayon::cyan('PCA completed\n'))
    
    object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,n.pcs])
    
  }
  
  return(object)
  
}

celseq_comb <- perform.pca(object = celseq_comb, 
                         assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                         n.pcs = 1:50, reduction.save = 'pca')
celseq2 <- perform.pca(object = celseq2, 
                       assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                       n.pcs = 1:50, reduction.save = 'pca')
celseq_comb <- perform.pca(object = celseq_comb, 
                           assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                           n.pcs = 1:50, reduction.save = 'pca')

perform.dbmap <- function(object, 
                          assay, 
                          slot='normalised',
                          n_components = 100, 
                          n_neighbors = 15, 
                          reduction.save='dbmap',
                          save.object = TRUE) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s)\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan('assay: ', x, ' does not exist\n'))
      return(object)
      
    }
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist \n'))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    cat(crayon::cyan('n_components must be numerical \n'))
    return(object)
    
  }
  
  if(!is.numeric(n_neighbors)) {
    
    cat(crayon::cyan('n_neighbors must be numerical \n'))
    return(object)
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string(s)\n'))
    return(object)
    
  }
  
  if(!is.logical(save.object)) {
    
    cat(crayon::cyan('save.object must be logical: TRUE/FALSE \n'))
    return(object)
    
  }
  
  if(is.null(reticulate::import('scipy.sparse', convert = FALSE))) {
    
    cat(crayon::cyan('scipy.sparse is not installed\n'))
    return(object)
    
  }
  
  if(is.null(reticulate::import('dbmap', convert = FALSE))) {
    
    cat(crayon::cyan('dbmap is not installed\n'))
    return(object)
    
  }
  
  scipy.sparse <- reticulate::import('scipy.sparse', convert = FALSE)
  
  dbmap <- reticulate::import('dbmap', convert = FALSE)
  
  for(o in assay) {
    
    cat(crayon::cyan(paste0('calculating dbmap for assay: ', o,'\n')))
    
    cellnames <- colnames(object@methods[[o]][[slot]])
    
    data <- scipy.sparse$csr_matrix(reticulate::r_to_py(t(as.matrix(object@methods[[o]][[slot]]))[,object@methods[[o]]@highly.variable.genes]))
    
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

celseq_comb <- perform.dbmap(object = celseq_comb, 
                           assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'),  
                           reduction.save = 'dbmap')
celseq2 <- perform.dbmap(object = celseq2, 
                         assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'),  
                         reduction.save = 'dbmap')
celseq_comb <- perform.dbmap(object = celseq_comb, 
                             assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'),  
                             reduction.save = 'dbmap')

perform.bbknn <- function(object,
                          assay,
                          reduction,
                          graph.name,
                          batch,
                          n_pcs = NULL,
                          trim = NULL,
                          n_trees = 10,
                          use_faiss = TRUE,
                          set_op_mix_ratio = 1.0,
                          local_connectivity= 1,
                          save.object = TRUE) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan('reduction must be character string\n'))
    return(object)
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.character(graph.name)) {
    
    cat(crayon::cyan('graph.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.character(batch)) {
    
    cat(crayon::cyan('batch must be character string\n'))
    return(object)
    
  }
  
  if(!batch %in% colnames(object@sample_metadata)) {
    
    cat(crayon::cyan('batch does not exist\n'))
    return(object)
    
  }
  
  if(!is.null(trim)) {
    
    if(!is.numeric(trim)) {
      
      cat(crayon::cyan('trim must be character string\n'))
      return(object)
      
    }
    
  }
  
  if(!is.numeric(n_trees)) {
    
    cat(crayon::cyan('n_trees must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(use_faiss)) {
    
    cat(crayon::cyan('use_faiss must be logical: TRUE/FALSE \n'))
    
  }
  
  if(!is.numeric(set_op_mix_ratio)) {
    
    cat(crayon::cyan('set_op_mix_ratio must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(local_connectivity)) {
    
    cat(crayon::cyan('local_connectivity must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(save.object)) {
    
    cat(crayon::cyan('save.object must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  for(p in assay) {
    
    count <- 1
    
    for(r in reduction) {
      
      sc <- reticulate::import('scanpy')
      
      scobj <- sc$AnnData(X = object@methods[[p]]@computational_reductions[[reduction[count]]])
      
      scobj$obs_names <- as.factor(colnames(object))
      
      scobj$var_names <- as.factor(colnames(object@methods[[p]]@computational_reductions[[reduction[count]]]))
      
      scobj$obsm$update(X_pca = object@methods[[p]]@computational_reductions[[reduction[count]]])
      
      if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
        
        pd <- reticulate::import('pandas')
        
        scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))
        
      }
      
      if(is.null(n_pcs)) {
        
        cat(crayon::cyan('npcs calculated\n'))
        
        n_pcs <- as.integer(length(colnames(object@methods[[p]]@computational_reductions[[reduction[count]]])))
        
      }
      
      if(is.null(trim)) {
        
        cat(crayon::cyan('initialising bbknn\n'))
        sc$external$pp$bbknn(scobj,
                             batch_key = as.character(batch),
                             approx = as.logical(FALSE),
                             n_pcs = n_pcs,
                             n_trees = as.integer(n_trees),
                             use_faiss = as.logical(use_faiss),
                             set_op_mix_ratio = set_op_mix_ratio,
                             local_connectivity = local_connectivity)
        
      } else if (!is.null(trim)) {
        
        sc$external$pp$bbknn(scobj,
                             batch_key= as.character(batch),
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
      connectivities <- as(object = as.matrix(connectivities), Class = 'dgCMatrix')
      distances <- as(object = as.matrix(distances), Class = 'dgCMatrix')
      
      graph.list[['connectivities']] <- connectivities
      graph.list[['distances']] <- distances
      
      object@methods[[p]]@graphs[[graph.name[count]]] <- graph.list
      object@methods[[p]]@alt_objects$anndata <- scobj
      count <- count + 1
      
    }
    
  }
  
  return(object)
  
}

celseq_comb <- perform.bbknn(object = celseq_comb, 
                             assay = c('SCT', 'SCANPY', 'TPM', 'SCRAN'), 
                             reduction = c('pca', 'dbmap'),
                             graph.name = c('pca_bbknn', 'dbmap_bbknn'),
                             batch = 'original.project')

library(harmony)

perform.harmony <- function(object, 
                            assay, 
                            vars.use, 
                            reduction = 'pca', 
                            reduction.save = 'harmony',
                            dims.use = NULL, 
                            ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(vars.use)) {
    
    cat(crayon::cyan('vars.use must be character string(s)\n'))
    return(object)
    
  }
  
  if(!vars.use %in% names(object@sample_metadata)) {
    
    cat(crayon::cyan('vars.use does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan('reduction must be character string(s)\n'))
    return(object)
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string(s)\n'))
    return(object)
    
  }
  
  if(is.list(dims.use)) {
    
    for(x in dims.use) {

      if(!is.null(x)) {
        
        if(!is.numeric(x)) {
          
          cat(crayon::cyan('dims must either be numeric or NULL\n'))
          return(object)
          
        } 
        
      } else if (!is.numeric(x)) {
        
        if(!is.null(x)) {
          
          cat(crayon::cyan('dims must either be numeric or NULL\n'))
          return(object)
          
        }
        
      }
      
    }
    
  }
  
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
    
    count <- 1
    
    for(g in reduction) {
      
      red <- reduction.list[[g]]
      
      red.save <- reduction.save[[count]]
      
      dims <- dims.use[[count]]
      
      if(is.null(dims)) {
        
        dims <- 1:ncol(red)
        
      }
      
      cat(crayon::cyan('Initialising harmony\n'))
      
      harm <- harmony::HarmonyMatrix(data_mat = red[,dims], meta_data = object@sample_metadata, vars_use = vars.use, do_pca = FALSE, verbose = TRUE, ...)
      
      object@methods[[p]]@integration_reductions[[red.save]] <- harm
      
      count <- count + 1
      
    }

  }

  cat(crayon::cyan('Harmony completed\n'))
  return(object)
}

celseq_comb <- perform.harmony(object = celseq_comb, 
                               assay = c('SCRAN', 'SCT', 'TPM', 'SCANPY'), 
                               vars.use = 'original.project', 
                               reduction = c('pca', 'dbmap'), 
                               reduction.save = c('pca_harmony', 'dbmap_harmony'), 
                               dims.use = list(1:8, NULL))

perform.scanorama <- function(object, 
                              assay,
                              slot,
                              split.by, 
                              n.dims = 50, 
                              reduction.save='scanorama', 
                              batch_size = 5000, 
                              approx = TRUE, 
                              sigma = 15, 
                              alpha = 0.1, 
                              knn = 20,
                              union = FALSE,
                              seed = 12345) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan('slot must be a character string\n'))
    return(object)
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    cat(crayon::cyan('slot does not exist\n'))
    return(object)
    
  }
  
  if(!is.character(split.by)) {
    
    cat(crayon::cyan('split.by must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(n.dims)) {
    
    cat(crayon::cyan('n.dims must be numerical\n'))
    return(object)
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(batch_size)) {
    
    cat(crayon::cyan('batch_size must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(approx)) {
    
    cat(crayon::cyan('approx must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.numeric(sigma)) {
    
    cat(crayon::cyan('sigma must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(alpha)) {
    
    cat(crayon::cyan('alpha must be numerical\n'))
    return(object)
    
  }
  
  if(!is.numeric(knn)) {
    
    cat(crayon::cyan('knn must be numerical\n'))
    return(object)
    
  }
  
  if(!is.logical(union)) {
    
    cat(crayon::cyan('union must be logical: TRUE/FALSE\n'))
    return(object)
    
  }
  
  if(!is.numeric(seed)) {
    
    cat(crayon::cyan('seed must be numerical\n'))
    return(object)
    
  }
  
  cat(crayon::cyan('Initialising scanorama\n'))
  scanorama <- reticulate::import('scanorama', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  
  count <- 1
  
  for(p in assay) {
    
    cat(crayon::cyan(paste0('Correcting assay: ', p, '\n')))
    
    list.matrix <- list()
    column.names <- list()
    sep <- unique(object@sample_metadata[,split.by])
    mat <- object@methods[[p]][[slot]]
    counter <- 1
    
    for(x in sep) {
      column.names[[counter]] <- colnames(mat[,object@sample_metadata[,split.by] == x])
      list.matrix[[counter]] <- t(mat[,object@sample_metadata[,split.by] == x])
      counter <- counter + 1
    }
    
    cat(crayon::cyan('Matrices isolated\n'))
    gene.list <- list()
    
    for(x in 1:length(sep)) {
      gene.list[[x]] <- rownames(mat[,object@sample_metadata[,split.by] == x])
    }
    
    cat(crayon::cyan('Genes identified\n'))
    cat(crayon::cyan('Corrections starting\n'))
    integrated.corrected.data <- scanorama$correct(datasets_full = reticulate::r_to_py(list.matrix), 
                                                   genes_list = reticulate::r_to_py(gene.list), 
                                                   dimred = as.integer(n.dims), 
                                                   return_dimred=TRUE, 
                                                   return_dense=FALSE, 
                                                   verbose = TRUE, 
                                                   batch_size = as.integer(batch_size), 
                                                   approx = approx, 
                                                   sigma = as.integer(sigma), 
                                                   alpha = as.numeric(alpha), 
                                                   knn = as.integer(knn),
                                                   union = as.logical(union),
                                                   seed = as.integer(seed))
    
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
    object@methods[[p]]@integration_reductions[[reduction.save]] <- t(combined)
    
  }
  
  return(object)

}

celseq_comb <- perform.scanorama(object = celseq_comb, 
                                 assay = c('SCT', 'TPM', 'SCRAN', 'SCANPY'), 
                                 slot = 'norm.scaled', 
                                 split.by = 'original.project', 
                                 n.dims = 50, 
                                 reduction.save = c('scanorama'))

perform.umap <- function(object, 
                         assay,
                         reduction=NULL,
                         graph=NULL,
                         reduction.save='umap',
                         n.dims=NULL, 
                         n_components = 3, 
                         ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string \n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan('assay: ', x, 'does not exist\n'))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string \n'))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    cat(crayon::cyan('n_components must be numerical\n'))
    return(object)
    
  }
  
  if(!is.null(reduction) & !is.null(graph)) {
    
    cat(crayon::cyan('only graphs OR reductions can be provided\n'))
    return(object)
    
  }
  
  for(u in assay) {
    
    if(is.null(reduction) & !is.null(graph)) {
      
      count <- 1
      
      graphs <- object@methods[[u]]@graphs
      
      for (g in graph) {
        
        cat(crayon::cyan('Processing', g, 'for assay:', u,'\n'))
        
        red.save <- reduction.save[count]
        
        seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[u]]@counts))
        
        seuobj[['temp']] <- suppressWarnings(Seurat::as.Graph(object@methods[[u]]@graphs[[g]]$connectivities))
        
        seuobj <- suppressWarnings(Seurat::RunUMAP(object = seuobj, 
                                                   graph = 'temp', 
                                                   n_components = n_components, 
                                                   verbose = TRUE,
                                                   ...))
        
        red.iso <- Seurat::Embeddings(object = seuobj, 
                                      reduction = 'umap')
        
        dim.names <- list()
        
        for(l in 1:2) {
          
          dim.names[[l]] <- paste0('umap_', l)
          
        }
        
        colnames(red.iso) <- unlist(dim.names)
        
        object@methods[[u]]@visualisation_reductions[[red.save]] <- red.iso
        
        count <- count + 1
        
      }
      
    } else {
      
      if(!is.list(n.dims)) {
        
        cat(crayon::cyan('dimensions must be supplied in list format\n'))
        return(object)
        
      }
      
      if(!is.character(reduction)) {
        
        cat(crayon::cyan('reduction must be character string \n'))
        return(object)
        
      }
      
      for(r in reduction) {
        
        if(!r %in% c(names(object@methods[[u]]@computational_reductions), 
                     names(object@methods[[u]]@integration_reductions),
                     names(object@methods[[u]]@visualisation_reductions))) {
          
          cat(crayon::cyan(paste0('reduction:', r, ' does not exist\n')))
          return(object)
          
        }
        
      }
      
      reduction.list <- list()
      
      red.names <- c(names(object@methods[[u]]@computational_reductions), 
                     names(object@methods[[u]]@integration_reductions),
                     names(object@methods[[u]]@visualisation_reductions))
      
      for(i in red.names) {
        
        if(i %in% names(object@methods[[u]]@computational_reductions)) {
          
          reduction.list[[i]] <- object@methods[[u]]@computational_reductions[[i]]
          
        }
        
        if(i %in% names(object@methods[[u]]@integration_reductions)) {
          
          reduction.list[[i]] <- object@methods[[u]]@integration_reductions[[i]]
          
        }
        
        if(i %in% names(object@methods[[u]]@visualisation_reductions)) {
          
          reduction.list[[i]] <- object@methods[[u]]@visualisation_reductions[[i]]
          
        }
        
      }
      
      count <- 1
      
      for(i in reduction) {
        
        dim <- n.dims[[count]]
        
        red.save <- reduction.save[count]
        
        red <- reduction.list[[i]]
        
        cat(crayon::cyan('Processing', i, 'for assay:', u,'\n'))
        
        if(!is.null(dim)) {
          
          c <- uwot::umap(X = red[,dim], n_components = n_components, verbose = TRUE, ...)
          
        } else {
          
          c <- uwot::umap(X = red, n_components = n_components, verbose = TRUE, ...)
          
        }
        
        dim.names <- list()
        
        for(l in 1:n_components) {
          
          dim.names[[l]] <- paste0('umap_', l)
          
        }
        
        colnames(c) <- unlist(dim.names)
        
        object@methods[[u]]@visualisation_reductions[[red.save]] <- c
        
        count <- count + 1
        
      }
      
    }
    
  }
  
  return(object)
  
}

# celseq_comb <- perform.umap(object = celseq_comb, 
#                             assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
#                             reduction = 'pca', 
#                             reduction.save = 'pca_umap', 
#                             n.dim = 1:8, n_components = 3)
# celseq_comb <- perform.umap(object = celseq_comb, 
#                             assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
#                             reduction = 'dbmap', 
#                             reduction.save = 'dbmap_umap', 
#                             n_components = 3)

celseq_comb <- perform.umap(object = celseq_comb, 
                        assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                        reduction = c('pca', 'dbmap'), 
                        reduction.save = c('pca_umap', 'dbmap_umap'), 
                        n_components = 3, 
                        n.dims = list(1:8, NULL))

plot.reduced.dim(object = smartseq2, 
                 reduction = 'pca_umap', 
                 assay = 'SCT', 
                 clust.method = 'metadata', 
                 column = 'celltype', 
                 pt.size = 2, 
                 dimensions = 2)


perform.tsne <- function(object, 
                         assay,
                         reduction=NULL,
                         reduction.save, 
                         n.dim=NULL,
                         n_components = 3, 
                         ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Must be IBRAP class object\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan(paste0('assay must be character string\n')))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan('assay: ', x, ' does not exist \n'))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan(paste0('reduction must be character string\n')))
    return(object)
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan(paste0('reduction.save must be character string\n')))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    cat(crayon::cyan(paste0('n_components must be numerical\n')))
    return(object)
    
  }
  
  if(!is.list(n.dim)) {
    
    cat(crayon::cyan('dimensions must be supplied in list format\n'))
    return(object)
    
  }
  
  for(g in assay) {
    
    reduction.list <- list()
    
    red.names <- c(names(object@methods[[g]]@computational_reductions), 
                   names(object@methods[[g]]@integration_reductions),
                   names(object@methods[[g]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[g]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[g]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[g]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[g]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[g]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[g]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    count <- 1
    
    for(r in reduction) {
      
      red <- reduction.list[[r]]
      
      red.save <- reduction.save[count]
      
      dim <- n.dim[[count]]
      
      cat(crayon::cyan('Processing', r, 'for assay:', g,'\n'))
      
      if(!is.null(dim)) {
        
        c <- ProjectionBasedClustering::tSNE(DataOrDistances = red[,dim], 
                                             OutputDimension = n_components, Iterations = 1000, verbose = TRUE, ...)$ProjectedPoints
        
      } else {
        
        c <- ProjectionBasedClustering::tSNE(DataOrDistances = red, 
                                             OutputDimension = n_components, Iterations = 1000, verbose = TRUE, ...)$ProjectedPoints
        
      }
      
      cat(crayon::cyan('t-SNE reduction completed\n'))
      
      dim.names <- list()
      
      for(t in 1:n_components) {
        
        dim.names[[t]] <- paste0('tsne_', t)
        
      }
      
      colnames(c) <- unlist(dim.names)
      
      rownames(c) <- colnames(object)
      
      object@methods[[g]]@visualisation_reductions[[red.save]] <- c
      
      cat(crayon::cyan('t-SNE data added\n'))
      
      count <- count + 1
      
    }
    
  }
  
  return(object)
  
}

# celseq_comb <- perform.tsne(object = celseq_comb, 
#                             assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'),
#                             reduction = 'pca', 
#                             reduction.save = 'pca_tsne', 
#                             n.dim = 1:8, 
#                             n_components = 3)
# celseq_comb <- perform.tsne(object = celseq_comb, 
#                             assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
#                             reduction = 'dbmap', 
#                             reduction.save = 'dbmap_tsne', 
#                             n_components = 3)

celseq2 <- perform.tsne(object = celseq2, 
                        assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                        reduction = c('pca', 'dbmap'), 
                        reduction.save = c('pca_tsne_per10', 'dbmap_tsne_per10'), 
                        n_components = 3, 
                        n.dim = list(1:8, NULL), perplexity = 10)

perform.lvish <- function(object, 
                          assay,
                          reduction='pca',
                          reduction.save='lvish',
                          n.dim=NULL, 
                          n_components = 3, 
                          ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string\n'))
    return(object)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan('reduction must be character string\n'))
    return(object)
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.character(reduction.save)) {
    
    cat(crayon::cyan('reduction.save must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(n_components)) {
    
    cat(crayon::cyan('n_components must be numerical\n'))
    return(object)
    
  }
  
  if(!is.list(n.dim)) {
    
    cat(crayon::cyan('dimensions must be supplied in list format\n'))
    return(object)
    
  }
  
  for(u in assay) {
    
    reduction.list <- list()
    
    red.names <- c(names(object@methods[[u]]@computational_reductions), 
                   names(object@methods[[u]]@integration_reductions),
                   names(object@methods[[u]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[u]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[u]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[u]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[u]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[u]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[u]]@visualisation_reductions[[i]]
        
      }
    }
    
    count <- 1
    
    for(r in reduction) {
      
      red <- reduction.list[[r]]
      
      red.save <- reduction.save[count]
      
      dim <- n.dim[[count]]
      
      cat(crayon::cyan('Processing', r, 'for assay:', u,'\n'))
      
      if(!is.null(dim)) {
        
        c <- uwot::lvish(X = red[,dim], n_components = n_components, verbose = TRUE, ...)
        
      } else {
        
        c <- uwot::lvish(X = red, n_components = n_components, verbose = TRUE, ...)
        
      }
      
      dim.names <- list()
      
      for(l in 1:n_components) {
        
        dim.names[[l]] <- paste0('lvish_', l)
        
      }
      
      colnames(c) <- unlist(dim.names)
      
      object@methods[[u]]@visualisation_reductions[[red.save]] <- c
      
    }
    
  }
  
  return(object)
  
}

smartseq2 <- perform.lvish(object = smartseq2,
                           assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                           reduction = c('pca', 'dbmap'), 
                           reduction.save = c('pca_tsne', 'dbmap_tsne'), 
                           n_components = 3, 
                           n.dim = list(1:8, NULL))
# panc <- perform.lvish(object = panc, 
#                       assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
#                       reduction = 'dbmap', 
#                       reduction.save = 'dbmap_lvish', n_components = 3)
celseq_comb <- perform.lvish(object = celseq_comb, 
                             assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                             reduction = c('pca_harmony', 'dbmap_harmony', 'scanorama'), 
                             reduction.save = c('pca_harmony_lvish', 'dbmap_harmony_lvish', 'scanorama_lvish'), 
                             n_components = 3, n.dim = list(1:8, NULL, NULL))

perform.seurat.cluster <- function(object, 
                                   assay,
                                   reduction=NULL,
                                   graph=NULL,
                                   algorithm=1,
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
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(object)
    
  }
  
  
  
  if(!is.null(reduction)) {
    
    if(is.list(dims)) {
      
      if(length(dims) != length(reduction)) {
        
        cat(crayon::cyan('number of dims must match number of reduction(s)\n'))
        return(object)
        
      }
      
      if(length(dims) != length(assignment.df.name)) {
        
        cat(crayon::cyan('number of dims must match number of assignment.df.name(s)\n'))
        return(object)
        
      }
      
      if(!is.character(assay)) {
        
        cat(crayon::cyan('assay must be character string\n'))
        return(object)
        
      }
      
    }
    
    if(!is.list(dims)) {
      
      cat(crayon::cyan('dims must be provided in a list\n'))
      return(object)
      
    }
    
    if(!is.character(reduction)) {
      
      cat(crayon::cyan('reduction must be character string\n'))
      return(object)
      
    }
    
    for(x in reduction) {
      
      for(i in assay) {
        
        if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                           object@methods[[i]]@visualisation_reductions, 
                           object@methods[[i]]@integration_reductions))) {
          
          cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
          return(object)
          
        }
        
      }
      
    }
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('assay: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  

  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan('assignment.df.name must be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(res)) {
    
    cat(crayon::cyan('res must be numerical\n'))
    return(object)
    
  }

  if(!is.numeric(algorithm)) {
    
    cat(crayon::cyan('method must either be numeric\n'))
    return(object)
    
  }
  
  if(!algorithm %in% c(1,2,3,4)) {
    
    cat(crayon::cyan('method must either: 1 (original Louvain), 2 (Louvain with multilevel refinement), 3 (SLM) or 4 (Leiden)\n'))
    return(object)
    
  }
  
  if(!is.numeric(prune.SNN)) {
    
    cat(crayon::cyan('prune.SNN must either be numeric\n'))
    return(object)
    
  }
  
  if(!is.character(nn.method)) {
    
    cat(crayon::cyan('nn.method must either be character string\n'))
    return(object)
    
  }
  
  if(!is.character(annoy.metric)) {
    
    cat(crayon::cyan('annoy.metric must either be character string\n'))
    return(object)
    
  }
  
  if(!is.numeric(nn.eps)) {
    
    cat(crayon::cyan('nn.eps must either be numeric\n'))
    return(object)
    
  }
  
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
    
    if(!is.null(reduction)) {
      
      for(r in reduction) {
        
        if(!r %in% names(reduction.list)) {
          
          cat(crayon::cyan('reductions could not be found\n'))
          return(object)
          
        }
        
      }
      
    }
    
    if(!is.null(graph)) {
      
      for(g in graph) {
        
        if(!g %in% names(object@methods[[p]]@graphs)) {
          
          cat(crayon::cyan('graphs could not be found\n'))
          return(object)
          
        }
        
      }
      
    }
    
    count <- 1
    
    if(is.null(reduction) & !is.null(graph)) {
      
      for(t in graph) {
        
        graph.iso <- object@methods[[p]]@graphs[[graph]]$connectivities
        
        cat(crayon::cyan(paste0('Calculating Seurat clusters for assay: ', p, ' using graph: ', t, '\n')))
        
        tmp <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
        
        orig.name <- names(tmp@meta.data)
        
        tmp[['temp']] <- Seurat::as.Graph(x = graph.iso)
        
        cat(crayon::cyan('Graph added\n'))
        
        tmp <- suppressWarnings(Seurat::FindClusters(object = tmp, resolution = res, graph.name = 'temp', algorithm = algorithm, verbose = FALSE, ...))
        
        cat(crayon::cyan('Clusters identified\n'))
        
        post.name <- names(tmp@meta.data)
        
        post.name <- post.name[!post.name %in% orig.name]
        
        post.name <- post.name[1:length(post.name)-1]
        
        new.clusters <- tmp@meta.data[,post.name]
        
        new.clusters <- new.clusters[match(rownames(object@sample_metadata), rownames(new.clusters)),]
        
        object@methods[[p]]@cluster_assignments[[assignment.df.name[count]]] <- new.clusters
        
        cat(crayon::cyan('Seurat clusters added\n'))
        
        count <- count + 1
        
      }
      

      
    } else {
      
      for(g in reduction) {
        
        dim <- dims[[count]]
        red <- reduction.list[[g]]
        red.key <- reduction[count]
        
        cat(crayon::cyan(paste0('Calculating Seurat clusters for assay: ', p, ' using reduction: ', g, '\n')))
        
        tmp <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
        
        tmp@reductions[[red.key]] <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = as.matrix(red), key = red.key))
        
        orig.name <- names(tmp@meta.data)
        
        if(!is.null(dim)) {
          
          tmp <- suppressWarnings(Seurat::FindNeighbors(object = tmp, 
                                                        reduction = red.key, 
                                                        verbose = FALSE, 
                                                        dims = dim, 
                                                        compute.SNN = TRUE, 
                                                        prune.SNN = prune.SNN,
                                                        nn.method = nn.method, 
                                                        annoy.metric = annoy.metric, 
                                                        nn.eps = nn.eps))
          
        } else {
          
          tmp <- suppressWarnings(Seurat::FindNeighbors(object = tmp, 
                                                        reduction = red.key, 
                                                        verbose = FALSE, 
                                                        dims = 1:ncol(red), 
                                                        compute.SNN = TRUE, 
                                                        prune.SNN = prune.SNN,
                                                        nn.method = nn.method, 
                                                        annoy.metric = annoy.metric, 
                                                        nn.eps = nn.eps))
          
        }
        
        cat(crayon::cyan('Neighbours identified\n'))
        
        tmp <- suppressWarnings(Seurat::FindClusters(object = tmp, resolution = res, algorithm = algorithm, verbose = FALSE, ...))
        
        cat(crayon::cyan('Clusters identified\n'))
        
        post.name <- names(tmp@meta.data)
        
        post.name <- post.name[!post.name %in% orig.name]
        
        post.name <- post.name[1:length(post.name)-1]
        
        new.clusters <- tmp@meta.data[,post.name]
        
        new.clusters <- new.clusters[match(rownames(object@sample_metadata), rownames(new.clusters)),]
        
        object@methods[[p]]@cluster_assignments[[assignment.df.name[count]]] <- new.clusters
        
        cat(crayon::cyan('Seurat clusters added\n'))
        
        count <- count + 1
      
        }
    
      }
    
    }
  
  return(object)
  
}

smartseq2 <- perform.seurat.cluster(object = smartseq2, 
                                    algorithm = 2,
                                    assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                    reduction = c('pca', 'dbmap'), 
                                    assignment.df.name = c('pca_LouvainMR', 'dbmap_LouvainMR'), 
                                    dims = list(1:8, NULL))
smartseq2 <- perform.seurat.cluster(object = smartseq2, 
                                    algorithm = 4,
                                    assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                    reduction = c('pca', 'dbmap'), 
                                    assignment.df.name = c('pca_Leiden', 'dbmap_Leiden'), 
                                    dims = list(1:8, NULL))
smartseq2 <- perform.seurat.cluster(object = smartseq2, 
                                    algorithm = 3,
                                    assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                    reduction = c('pca', 'dbmap'), 
                                    assignment.df.name = c('pca_SLM', 'dbmap_SLM'), 
                                    dims = list(1:8, NULL))
smartseq2 <- perform.seurat.cluster(object = smartseq2, 
                                    algorithm = 1,
                                    assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                    reduction = c('pca', 'dbmap'), 
                                    assignment.df.name = c('pca_Louvain', 'dbmap_Louvain'), 
                                    dims = list(1:8, NULL))

celseq2 <- perform.seurat.cluster(object = celseq2, 
                                  algorithm = 2,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca', 'dbmap'), 
                                  assignment.df.name = c('pca_LouvainMR', 'dbmap_LouvainMR'), 
                                  dims = list(1:8, NULL))
celseq2 <- perform.seurat.cluster(object = celseq2, 
                                  algorithm = 4,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca', 'dbmap'), 
                                  assignment.df.name = c('pca_Leiden', 'dbmap_Leiden'), 
                                  dims = list(1:8, NULL))
celseq2 <- perform.seurat.cluster(object = celseq2, 
                                  algorithm = 3,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca', 'dbmap'), 
                                  assignment.df.name = c('pca_SLM', 'dbmap_SLM'), 
                                  dims = list(1:8, NULL))
celseq2 <- perform.seurat.cluster(object = celseq2, 
                                  algorithm = 1,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca', 'dbmap'), 
                                  assignment.df.name = c('pca_Louvain', 'dbmap_Louvain'), 
                                  dims = list(1:8, NULL))

celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                  algorithm = 2,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca_harmony', 'dbmap_harmony', 'scanorama'), 
                                  assignment.df.name = c('pca_harmony_LouvainMR', 'dbmap_harmony_LouvainMR', 'scanorama_LovainMR'), 
                                  dims = list(1:8, NULL, NULL))
celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                  algorithm = 4,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca_harmony', 'dbmap_harmony', 'scanorama'), 
                                  assignment.df.name = c('pca_harmony_Leiden', 'dbmap_harmony_Leiden', 'scanorama_Leiden'), 
                                  dims = list(1:8, NULL, NULL))
celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                  algorithm = 3,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca_harmony', 'dbmap_harmony', 'scanorama'), 
                                  assignment.df.name = c('pca_harmony_SLM', 'dbmap_harmony_SLM', 'scanorama_SLM'), 
                                  dims = list(1:8, NULL, NULL))
celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                  algorithm = 1,
                                  assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                  reduction = c('pca_harmony', 'dbmap_harmony', 'scanorama'), 
                                  assignment.df.name = c('pca_harmony_Louvain', 'dbmap_harmony_Louvain', 'scanorama_Louvain'), 
                                  dims = list(1:8, NULL, NULL))

celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                      algorithm = 1,
                                      assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                      graph = c('pca_bbknn', 'dbmap_bbknn'),
                                      assignment.df.name = c('pca_bbknn_Louvain', 'dbmap_bbknn_Louvain'), 
                                      dims = list(1:8, NULL, NULL))
celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                      algorithm = 2,
                                      assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                      graph = c('pca_bbknn', 'dbmap_bbknn'),
                                      assignment.df.name = c('pca_bbknn_LouvainMR', 'dbmap_bbknn_LouvainMR'), 
                                      dims = list(1:8, NULL, NULL))
celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                      algorithm = 3,
                                      assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                      graph = c('pca_bbknn', 'dbmap_bbknn'),
                                      assignment.df.name = c('pca_bbknn_SLM', 'dbmap_bbknn_SLM'), 
                                      dims = list(1:8, NULL, NULL))
celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
                                      algorithm = 4,
                                      assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                                      graph = c('pca_bbknn', 'dbmap_bbknn'),
                                      assignment.df.name = c('pca_bbknn_Leiden', 'dbmap_bbknn_Leiden'), 
                                      dims = list(1:8, NULL, NULL))




# celseq_comb@methods$SCANPY@clus .ter_assignments$pca_LouvainMR <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$dbmap_LouvainMR <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$pca_SLM <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$dbmap_SLM <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$pca_Leiden <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$dbmap_Leiden <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$pca_Louvain <- NULL
# celseq_comb@methods$SCANPY@cluster_assignments$dbmap_Louvain <- NULL

# 
# celseq_comb <- perform.seurat.cluster(object = smartseq2, 
#                        assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
#                        graph = 'pca_bbknn', 
#                        algorithm = 2, 
#                        assignment.df.name = c('pca_bbknn_LouvainML'))
# celseq_comb <- perform.seurat.cluster(object = celseq_comb, 
#                                       assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
#                                       graph = 'pca_bbknn', 
#                                       algorithm = 4, 
#                                       assignment.df.name = c('pca_bbknn_Leiden'))


library(SC3)
library(SingleCellExperiment)

perform.sc3.reduction.cluster <- function(object, 
                                          assay,
                                           reduction,
                                          dims,
                                          assignment.df.name,
                                          ks, 
                                          n.core = 3) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(NULL)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s) \n'))
    return(NULL)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan(paste0('assignment.df.name must be character string(s)\n')))
    return(object)
    
  }
  
  if(!is.numeric(ks)) {
    
    cat(crayon::cyan(paste0('ks must be numerical\n')))
    return(object)
    
  }
  
  if(!is.numeric(n.core)) {
    
    cat(crayon::cyan(paste0('n.core must be numerical\n')))
    return(object)
    
  }
  
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
  
    reduction.list <- reduction.list[reduction]
    
    count <- 1
    
    for(r in reduction) {

      red <- reduction.list[[r]]
      
      dimen <- dims[[count]]
      
      if(is.null(dimen)) {
        
        dimen <- 1:ncol(red)
        
      }
      
      temp.2 <- SingleCellExperiment(list('logcounts' = t(red)[dimen,]))
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
      temp.2 <- sc3_calc_dists(temp.2)
      temp.2 <- sc3_calc_transfs(temp.2)
      temp.2 <- sc3_kmeans(temp.2, ks = ks)
      temp.2 <- sc3_calc_consens(temp.2)
      object@methods[[p]]@cluster_assignments[[assignment.df.name[[count]]]] <- as.data.frame(colData(temp.2))
      cat(crayon::cyan('SC3 clustering completed\n'))
      count <- count + 1
      
    }
    
  }

  return(object)
  
}

test <- perform.sc3.reduction.cluster(object = celseq_comb, 
                              assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                              reduction = c('pca_harmony', 'dbmap_harmony', 'scanorama'), 
                              dims = list(NULL, NULL, NULL), 
                              assignment.df.name = c('pca_harmony_SC3', 'dbmap_harmony_SC3', 'scanorama_SC3'), 
                              ks = 10:16, 
                              n.core = 3)

perform.sc3.slot.cluster <- function(object, 
                                     assay,
                                     slot,
                                     HVGs,
                                     assignment.df.name,
                                     ks, 
                                     n.core = 3) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(NULL)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s) \n'))
    return(NULL)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    cat(crayon::cyan(paste0('slot must be character string\n')))
    return(object)
    
  }
  
  if(!is.logical(HVGs)) {
    
    cat(crayon::cyan(paste0('HVGs must be logical: TRUE/FALSE\n')))
    return(object)
    
  }
  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan(paste0('assignment.df.name must be character string(s)\n')))
    return(object)
    
  }
  
  if(!is.numeric(ks)) {
    
    cat(crayon::cyan(paste0('ks must be numerical\n')))
    return(object)
    
  }
  
  if(!is.numeric(n.core)) {
    
    cat(crayon::cyan(paste0('n.core must be numerical\n')))
    return(object)
    
  }
  
  cat(crayon::cyan('Initialising SC3 clustering\n'))
  
  for(p in assay) {
    
    mat <- object@methods[[p]][[slot]]
    
    if(!is.null(HVGs)) {
      
      mat <- mat[object@methods[[p]]@highly.variable.genes,]
      
    }
    

    if(is.null(HVGs)) {
      
      temp.2 <- SingleCellExperiment(list('counts' = as.matrix(object@methods[[p]]@counts[rownames(mat),]), 'logcounts' = as.matrix(mat)))
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- sc3_prepare(temp.2, gene_filter = TRUE, n_cores = n.core)
      
    } else {
      
      temp.2 <- SingleCellExperiment(list('logcounts' = as.matrix(mat)))
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
      
    }
    
    temp.2 <- sc3_calc_dists(temp.2)
    temp.2 <- sc3_calc_transfs(temp.2)
    temp.2 <- sc3_kmeans(temp.2, ks = ks)
    temp.2 <- sc3_calc_consens(temp.2)
    object@methods[[p]]@cluster_assignments[[assignment.df.name]] <- as.data.frame(colData(temp.2))
    cat(crayon::cyan('SC3 clustering completed\n'))
    
  }
  
  return(object)
  
}

library(SC3)

smartseq2 <- perform.sc3.slot.cluster(object = smartseq2, 
                                 assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                                 slot = 'normalised', 
                                 HVGs = T, 
                                 assignment.df.name = 'counts_SC3', 
                                 ks = 10:16, n.core = 3)
celseq2 <- perform.sc3.slot.cluster(object = celseq2, 
                                    assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                                    slot = 'normalised', 
                                    HVGs = T, 
                                    assignment.df.name = 'counts_SC3', 
                                    ks = 10:16, n.core = 3)

perform.reduction.kmeans <- function(object, 
                                assay,
                                reduction=NULL, 
                                dims = NULL,
                                k=NULL,
                                assignment.df.name,
                                method='kmeans',
                                ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(NULL)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s) \n'))
    return(NULL)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('reduction: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
        return(object)
        
      }
      
    }
    
  }
  
  if(!is.numeric(k)) {
    
    cat(crayon::cyan(paste0('k must be numerical\n')))
    return(object)
    
  }
  
  if(!is.character(assignment.df.name)) {
    
    cat(crayon::cyan(paste0('assignment.df.name must be character string(s)\n')))
    return(object)
    
  }
  
  if(!method %in% c('kmeans', 'pam')) {
        
        cat(crayon::cyan(paste0('method must be kmeans or pam\n')))
        return(object)
        
      }
  
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
    
    print(reduction)
    print(names(reduction.list))
    
    reduction.list <- reduction.list[reduction]
    
    count <- 1
    
    for(o in reduction) {
      
      if(length(reduction) == 1) {
        
        red <- reduction.list
        
      } else {
        
        red <- reduction.list[[o]] 
        
      }
      
      dimen <- dims[[count]]
      
      if(is.null(dimen)) {
        
        dimen <- 1:length(colnames(red))
        
      }
      
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
          clusters[,paste0('pam_clustering_K_', i)] <- as.factor(cluster::pam(x = red[,dimen], k = i, ...)$clustering)
        }
      }
      if(method == 'kmeans') {
        for(i in k) {
          clusters[[paste0('kmeans_clustering_K_', i)]] <- as.factor(kmeans(x = red[,dimen], centers = i, ...)$cluster)
        }
      } else {
        cat(crayon::cyan('Please specify method: pam or kmeans\n'))
      }
      
      clusters <- clusters[,2:length(colnames(clusters))]
      object@methods[[p]]@cluster_assignments[[assignment.df.name[[count]]]] <- clusters
      
      count <- count + 1
    }

  }
  
  return(object)
}

smartseq2 <- perform.reduction.kmeans(object = smartseq2, 
                                      assay = c('SCT','SCRAN','SCANPY','TPM'), 
                                      reduction = c('dbmap_tsne_per30', 
                                                    'pca_tsne_per30'), 
                                      dims = list(NULL, NULL), 
                                      k = 10:16, 
                                      assignment.df.name = c('dbmap_tsne_kmeans', 
                                                             'pca_tsne_kmeans'), 
                                      method = 'kmeans')

celseq2 <- perform.reduction.kmeans(object = celseq2, 
                                    assay = c('SCT','SCRAN','SCANPY','TPM'), 
                                    reduction = c('dbmap_tsne_per30', 
                                                  'pca_tsne_per30'), 
                                    dims = list(NULL, NULL), 
                                    k = 10:16, 
                                    assignment.df.name = c('dbmap_tsne_kmeans', 
                                                           'pca_tsne_kmeans'), 
                                    method = 'kmeans')

benchmark.clustering <- function(object, 
                                 assay,
                                 clustering,
                                 reduction, 
                                 n.dims = 1:3,
                                 dist.method='euclidean',
                                 ground.truth=NULL) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP \n'))
    return(NULL)
    
  }
  
  if(!is.character(assay)) {
    
    cat(crayon::cyan('assay must be character string(s) \n'))
    return(NULL)
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      cat(crayon::cyan(paste0('assay: ', x, 'does not exist\n')))
      return(object)
      
    }
    
  }
  
  if(!is.character(clustering)) {
    
    cat(crayon::cyan('clustering must be character string(s)\n'))
    return(NULL)
    
  }
  
  for(l in assay) {
    
    for(x in clustering) {
      
      if(!x %in% names(object@methods[[l]]@cluster_assignments)) {
        
        cat(crayon::cyan('clustering: ,', x, 'not present in assay: ', l, '\n'))
        return(NULL)
        
      }
      
    }
    
    for(x in reduction) {
      
      for(i in assay) {
        
        if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                           object@methods[[i]]@visualisation_reductions, 
                           object@methods[[i]]@integration_reductions))) {
          
          cat(crayon::cyan(paste0('reduction: ', x, ' does not exist\n')))
          return(object)
          
        }
        
      }
      
    }
    
    if(!is.numeric(n.dims)) {
      
      cat(crayon::cyan('n.dims must be numerical\n'))
      return(object)
      
    }
    
    if(!is.character(dist.method)) {
      
      cat(crayon::cyan('dist.method must be a character string\n'))
      return(object)
      
    }
    
    
    
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
      
      reduction_sub <- reduction.list[[reduction[count]]][,n.dims]
      
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
        
        object@methods[[l]]@benchmark_results[[k]] <- as.data.frame(results)
        
      } else {
        
        results <- cbind(sil.results, dunn.results, conn.results)
        rownames(results) <- colnames(all.clusters)
        colnames(results) <- c(paste0(k, '_sil.results'), paste0(k, '_dunn.results'), paste0(k, '_conn.results'))
        object@methods[[l]]@benchmark_results[[k]] <- as.data.frame(results)
        
      }
    }
  }
  
  return(object)
  
}

smartseq2 <- benchmark.clustering(object = smartseq2, 
                                  assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                                  clustering = c('pca_LouvainMR',
                                                 'dbmap_LouvainMR',
                                                 'pca_Leiden',
                                                 'dbmap_Leiden',
                                                 'pca_SLM',
                                                 'dbmap_SLM',
                                                 'pca_Louvain',
                                                 'dbmap_Louvain',
                                                 'counts_SC3',
                                                 'dbmap_tsne_kmeans',
                                                 'pca_tsne_kmeans'), 
                                  reduction = c('pca_umap_a2b1', 
                                                'dbmap_umap_a2b1',
                                                'pca_umap_a2b1', 
                                                'dbmap_umap_a2b1',
                                                'pca_umap_a2b1', 
                                                'dbmap_umap_a2b1',
                                                'pca_umap_a2b1', 
                                                'dbmap_umap_a2b1',
                                                'pca_umap_a2b1', 
                                                'dbmap_umap_a2b1',
                                                'pca_umap_a2b1'), 
                                  components = 1:3, 
                                  dist.method = 'euclidean', 
                                  ground.truth = smartseq2$celltype)

celseq2 <- benchmark.clustering(object = celseq2, 
                                assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                                clustering = c('pca_LouvainMR',
                                               'dbmap_LouvainMR',
                                               'pca_Leiden',
                                               'dbmap_Leiden',
                                               'pca_SLM',
                                               'dbmap_SLM',
                                               'pca_Louvain',
                                               'dbmap_Louvain',
                                               'counts_SC3',
                                               'dbmap_tsne_kmeans',
                                               'pca_tsne_kmeans'), 
                                reduction = c('pca_umap', 
                                              'dbmap_umap',
                                              'pca_umap', 
                                              'dbmap_umap',
                                              'pca_umap', 
                                              'dbmap_umap',
                                              'pca_umap', 
                                              'dbmap_umap',
                                              'pca_umap', 
                                              'dbmap_umap',
                                              'pca_umap'), 
                                components = 1:3, 
                                dist.method = 'euclidean', 
                                ground.truth = celseq2$celltype)

panc <- benchmark.clustering(object = panc, 
                             assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                             clustering = c('dbmap_tsne_kmeans', 'pca_tsne_kmeans'), 
                             reduction = c('dbmap_umap', 'pca_umap'), 
                             components = 1:3, 
                             dist.method = 'euclidean', 
                             ground.truth = panc$celltype)

panc <- benchmark.clustering(object = panc, 
                             assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                             clustering = c('counts_SC3'), 
                             reduction = c('pca_umap'), 
                             components = 1:3, 
                             dist.method = 'euclidean', 
                             ground.truth = panc$celltype)

panc <- benchmark.clustering(object = celseq_comb, 
                             assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                             clustering = c("pca_harmony_LouvainMR",
                                            "dbmap_harmony_LouvainMR",
                                            "scanorama_LovainMR",
                                            "pca_harmony_Leiden",
                                            "dbmap_harmony_Leiden",
                                            "scanorama_Leiden",
                                            "pca_harmony_Louvain",
                                            "dbmap_harmony_Louvain",
                                            "scanorama_Louvain",
                                            "pca_harmony_SLM",
                                            "dbmap_harmony_SLM",
                                            "scanorama_SLM",
                                            "pca_harmony_SC3",
                                            "dbmap_harmony_SC3",
                                            "scanorama_SC3",
                                            "pca_harmony_tsne_kmeans",
                                            "scanorama_tsne_kmeans"), 
                             reduction = c('pca_harmony_umap',
                                           'dbmap_harmony_umap',
                                           'scanorama_umap',
                                           'pca_harmony_umap',
                                           'dbmap_harmony_umap',
                                           'scanorama_umap',
                                           'pca_harmony_umap',
                                           'dbmap_harmony_umap',
                                           'scanorama_umap',
                                           'pca_harmony_umap',
                                           'dbmap_harmony_umap',
                                           'scanorama_umap',
                                           'pca_harmony_umap',
                                           'dbmap_harmony_umap',
                                           'scanorama_umap',
                                           'pca_harmony_umap',
                                           'scanorama_umap'), 
                             components = 1:3, 
                             dist.method = 'euclidean', 
                             ground.truth = celseq_comb$celltype)


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
                             dimensions) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('object must be of class IBRAP\n'))
    return(NULL)
    
  }
  
  if(!is.character(reduction)) {
    
    cat(crayon::cyan('reduction must be character string\n'))
    return(NULL)
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    cat(crayon::cyan('assay does not exist\n'))
    return(NULL)
    
  }
  
  if(!reduction %in% names(c(object@methods[[assay]]@computational_reductions, 
                             object@methods[[assay]]@visualisation_reductions, 
                             object@methods[[assay]]@integration_reductions))) {
    
    cat(crayon::cyan('reduction does not exist\n'))
    return(NULL)
    
  }
  
  if(!is.character(clust.method)) {
    
    cat(crayon::cyan('clust.method must be character string\n'))
    return(NULL)
    
  }
  
  if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
    
    if(!clust.method == 'metadata') {
      
      cat(crayon::cyan('clust.method does not exist\n'))
      return(NULL)
      
    }
    
  }
  
  if(!is.character(column)) {
    
    cat(crayon::cyan('column must be character string\n'))
    return(NULL)
    
  }
  
  if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
    
    if(!column %in% colnames(object@sample_metadata)) {
      
      cat(crayon::cyan('column:', column, ', does not exist in clust.method: ', clust.method, '\n'))
      return(NULL)
      
    }
    
  }
  
  if(!is.numeric(pt.size)) {
    
    cat(crayon::cyan('pt.size must be numerical\n'))
    return(NULL)
    
  }
  
  if(!is.numeric(dimensions)) {
    
    cat(crayon::cyan('dimensions must be numerical: 2 or 3\n'))
    return(NULL)
    
  }
  
  reduction.list <- list()
  red.names <- c(names(object@methods[[assay]]@computational_reductions), 
                 names(object@methods[[assay]]@integration_reductions),
                 names(object@methods[[assay]]@visualisation_reductions))
  
  for(i in red.names) {
    
    if(i %in% names(object@methods[[assay]]@computational_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@computational_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@integration_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@integration_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@visualisation_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@visualisation_reductions[[i]]
      
    }
    
  }
  
  if(clust.method == 'metadata') {
    
    project.met <- object@sample_metadata
    
    print(head(project.met))
    
    results <- as.data.frame(reduction.list[[reduction]])
    
    if(dimensions == 2) {
      
      if(ncol(results) < 2) {
        
        cat(crayon::cyan('Not enough dimensions present in ', reduction, '\n'))
        return(NULL)
        
      }
      
      results <- results[,1:2]
      
    } else if (dimensions == 3) {
      
      if(ncol(results) < 3) {
        
        cat(crayon::cyan('Not enough dimensions present in ', reduction, '\n'))
        return(NULL)
        
      }
      
      results <- results[,1:3]
      
    }
    
    orig.names <- colnames(results)
    
    results <- cbind(results, project.met[,column])
    
    colnames(results) <- c(orig.names, 'variable')
    
    rownames(results) <- colnames(object)
    
  } else {
    
    project.met <- object@sample_metadata
    
    assay.met <- object@methods[[assay]]@cluster_assignments[[clust.method]]
    
    assay.met <- assay.met[match(rownames(project.met), rownames(assay.met)),]
    
    meta <- cbind(project.met, assay.met)
    
    results <- as.data.frame(reduction.list[[reduction]])
    
    if(dimensions == 2) {
      
      if(ncol(results) < 2) {
        
        cat(crayon::cyan('Not enough dimensions present in ', reduction, '\n'))
        return(NULL)
        
      }
      
      results <- results[,1:2]
      
    } else if (dimensions == 3) {
      
      if(ncol(results) < 3) {
        
        cat(crayon::cyan('Not enough dimensions present in ', reduction, '\n'))
        return(NULL)
        
      }
      
      results <- results[,1:3]
      
    }
    
    orig.names <- colnames(results)
    
    results <- cbind(results, meta[,column])
    
    colnames(results) <- c(orig.names, 'variable')
    
    rownames(results) <- colnames(object)
    
  }
  
  if(dimensions == 3) {
    
    p <- plotly::plot_ly(data = results,
                         x = as.formula(paste0('~', colnames(results)[1])), 
                         y = as.formula(paste0('~', colnames(results)[2])),
                         z = as.formula(paste0('~', colnames(results)[3])), 
                         color = as.formula(paste0('~',colnames(results)[4])), 
                         colors = colorspace::qualitative_hcl(n = length(unique(results[,4])), palette = 'Dark 3'),
                         mode = "markers", 
                         marker = list(size = pt.size, width=0.5), 
                         text=as.formula(paste0('~',colnames(results)[4])), 
                         hoverinfo="text", plot_bgcolor = 'black')
    
    print(p)
    
  } else if (dimensions == 2) {
    
    p <- plotly::plot_ly(data = as.data.frame(results), 
                         x = as.formula(paste0('~', colnames(results)[1])), 
                         y = as.formula(paste0('~', colnames(results)[2])), 
                         color = as.formula(paste0('~',colnames(results)[3])),
                         colors = colorspace::qualitative_hcl(n = length(unique(results[,3])), palette = 'Dark 3'), 
                         mode = "markers", 
                         marker = list(size = pt.size, width=0.5), 
                         text=as.formula(paste0('~',colnames(results)[3])), 
                         hoverinfo="text", plot_bgcolor = 'black')
    
    print(p)
    
  }
  
}

plot.reduced.dim(object = celseq_comb, reduction = 'trajectory_inference', assay = 'SCANPY', 
                 clust.method = 'metadata', column = 'celltype', dimensions = 2, pt.size = 2)

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
                                   slot,
                                   reduction,
                                   features) {
  
  plot.list <- list()
  
  for(x in features) {
    print('1')
    results <- as.data.frame(object@methods[[assay]]@visualisation_reductions[[reduction]])[,1:2]
    orig.colnames <- colnames(object@methods[[assay]]@visualisation_reductions[[reduction]][,1:2])
    print('2')
    print(assay)
    iso <- object@methods[[assay]][[slot]][x,]
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
      theme_bw() + labs(title=x, x=orig.colnames[1], y=orig.colnames[2]) + 
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

plot.features.multiple(object = panc, 
                       assay = 'SCT', 
                       reduction = 'pca_umap', 
                       slot = 'normalised', features = c('GCG', 'MAFA', 'MAFB'))
# 
# plot.barplot <- function(object, 
#                          x.value, 
#                          y.value) {
#   object[['var']] <- y.value 
#   object[['group']] <- x.value
#   print(length(unique(object[['var']])))
#   SingleCellExperiment()
#   p <- dittoSeq::dittoBarPlot(object = object, 
#                               var = 'var', 
#                               group.by = 'group') + 
#     labs(title='') + 
#     theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
#   
#   print(p)
# }

plot.barplot(object = panc,
             y.value = metadata(panc)$clustering[['seurat:pca_harmony']][['Seurat_res_0.2']],
             x.value = panc$tech)

plot.benchmarking <- function(object, 
                              assay,
                              clustering, 
                              ARI){
  clust.bench <- object@methods[[assay]]@benchmark_results[[clustering]]
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
  
  last.fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(last.label), group = 1)) +
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
                         slot,
                         features,
                         group.by,
                         ...) {

  ass <- object@methods[[assay]][[slot]]
  cat(crayon::cyan('Isolated assay\n'))
  ass <- as.matrix(ass)
  logged <- log(ass+1)
  means <- apply(X = logged, MARGIN = 1, FUN = mean)
  standev <- apply(X = logged, MARGIN = 1, FUN = sd)
  z_scores <- logged - means / standev
  cat(crayon::cyan('z-scores calculated\n'))
  
  seuobj <- Seurat::CreateSeuratObject(counts = z_scores)
  
  if(length(group.by) != ncol(seuobj)) {
    
    cat(crayon::cyan('group.by variable does not contain the same quantity as cells.\n'))
    
  }
  
  features <- features[features %in% rownames(seuobj)]
  
  seuobj$clusters <- group.by
  
  print(Seurat::DoHeatmap(object = seuobj, features = features, group.by = 'clusters', slot = 'counts'))
  
}

plot.heatmap(object = panc, assay = 'SCT', slot = 'normalised', features = features, group.by = panc@methods$SCT@cluster_assignments$counts_SC3$sc3_13_clusters)

plot.vln <- function(object, 
                     assay,
                     slot,
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
    ass <- t(object@methods[[assay]][[slot]][x,])
    df <- data.frame(barcodes = colnames(object))
    rownames(df) <- df$barcodes
    df[,'feature'] <- t(ass)
    df[,'group'] <- group.by
    
    p <- ggplot(data = df, aes(x = group, y = feature, color = group)) + 
      geom_violin() + geom_boxplot() +
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

plot.vln(object = panc, assay = 'SCT', 
         features = c('GCG', 'MAFA', 'MAFB'), 
         slot = 'normalised', 
         group.by = panc@methods$SCT@cluster_assignments$pca_seurat$RNA_snn_res.1.5)

perform.seurat.diffexp <- function(object, 
                                   assay,
                                   test, 
                                   identity,
                                   latent.vars,
                                   ...) {
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
  seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
  
  seuobj$clusters <- identity
  Seurat::Idents(seuobj) <- 'clusters'
  
  if(!is.null(latent.vars)) {
    
    if(!is.character(latent.vars)) {
      
      cat(crayon::cyan('latent.vars must be character(s)\n'))
      return(NULL)
      
    }
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, latent.vars = latent.vars, ...)
    
  } else {
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, ...)
    
  }
  
  return(results)
  
}

results_integrated_scanpy <- perform.seurat.diffexp(object = celseq_comb,
                                            assay = 'SCANPY',
                                            test = 'MAST', 
                                            latent.vars = c('original.project', 'RAW_total.counts'),
                                            identity = celseq_comb@methods$SCANPY@cluster_assignments$pca_bbknn_Louvain$pca_bbknn_res.0.3)

results_integrated_sct <- perform.seurat.diffexp(object = celseq_comb,
                                                    assay = 'SCT',
                                                    test = 'MAST', 
                                                    latent.vars = c('original.project', 'RAW_total.counts'),
                                                    identity = celseq_comb@methods$SCT@cluster_assignments$pca_bbknn_Leiden$pca_bbknn_res.0.2)

results_integrate_groundtruth <- perform.seurat.diffexp(object = celseq_comb,
                                                        assay = 'SCRAN',
                                                        test = 'MAST', 
                                                        latent.vars = c('original.project', 'RAW_total.counts'),
                                                        identity = celseq_comb@sample_metadata$celltype)

l <- celseq_comb@methods$SCANPY@benchmark_results
m <- l[[1]]
colnames(m) <- c('sil.results', 'dunn.results', 'conn.results', 'ARI.results', 'NMI.results')
rownames(m) <- paste0(names(l[1]), '_', rownames(m))

for(x in 2:length(l)) {
  
  h <- l[[x]]
  colnames(h) <- c('sil.results', 'dunn.results', 'conn.results', 'ARI.results', 'NMI.results')
  rownames(h) <- paste0(names(l[x]), '_', rownames(h))
  m <- rbind(m,h)
  
}

m$cluster_index <- rownames(m)
m$pipline <- rownames(m)
m$norm <- 'Scanpy'
colnames(m) <- c('ASW','Dunn Index', 'Connectivity', 'ARI', 'NMI', 'Pipeline', 'Normalisation')


m$cluster_index <- factor(x = m$cluster_index, levels = m$cluster_index[order(m$ARI.results, decreasing = T)])

ggplot(data = m, mapping = aes_string(x = 'cluster_index', y = 'ARI.results', group = 1)) +
  geom_point() +
  geom_line() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = 'Integrated SCANPY')
