#' @name perform.scran
#' @aliases perform.scran
#' 
#' @title Performs Scran normalisation, hvg selection, scaling and variance stabilisation and regression. 
#'
#' @description A new method-assay is produced. Raw counts are normalised and HVGs identified using Scran 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param batch Character. Which column in the metadata defines the batches. Default = NULL
#' @param vars.to.regress Character. Which column in the metadata should be regressed. Default = NULL
#' @param do.scale Boolean. Whether to scale the features variance. Default = TRUE
#' @param do.center Boolean. Whether to centre features to zero. Default = TRUE
#' @param new.assay.suffix Character. What should be added as a suffix for SCRAN
#' @param n.genes Numerical. Top number of genes to retain when finding HVGs. Default = 1500
#' @param max.cluster.size Numerical. When performing quickCluster, what is the maximum size the clusters can be. Default = 1000
#' @param center_size_factors Boolean Should size factor variance be centred. Default = TRUE
#' @param verbose Logical Should function messages be printed?
#' @param seed Numerical What seed should be set. Default = 1234
#' @param ... Arguments to pass to Seurat::ScaleData
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#' 
#' @examples 
#' 
#' object <- perform.scran(object = object, 
#'                         assay = 'RAW', 
#'                         slot = 'counts', 
#'                         vars.to.regress = 'RAW_total.counts', do.scale = T)
#'
#' @export

perform.scran <- function(object, 
                          assay = 'RAW',
                          slot = 'counts',
                          batch=NULL,
                          vars.to.regress=NULL,
                          do.scale=TRUE,
                          do.center=TRUE,
                          new.assay.suffix = '',
                          n.genes=1500,
                          max.cluster.size = 1000,
                          center_size_factors=TRUE,
                          verbose = FALSE,
                          seed=1234,
                          ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Object must be of class IBRAP\n')
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose must be logical, TRUE/FALSE\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('Assay must be a character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('Assay does not exist\n')
    
  }
  
  if(!is.character(slot)) {
    
    stop('Slot must be a character string\n')
    
  }
  
  if(!is.null(batch)) {
    
    stop('batch must be a character string of column name contained in sample_metadata\n')
    
  }
  
  if(!is.null(vars.to.regress)) {
    
    if(!is.character(vars.to.regress)) {
      
      stop('vars.to.regress must be a character string(s) of column name contained in sample_metadata\n')
      
    }
    
  }
  
  if(!is.logical(do.scale)) {
    
    stop('do.scale must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.logical(do.center)) {
    
    stop('do.center must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.character(new.assay.suffix)) {
    
    stop('new.assay.suffix must be character string\n')
    
  }
  
  if(!is.numeric(n.genes)) {
    
    stop('n.genes must be numerical\n')
    
  }
  
  if(!is.numeric(max.cluster.size)) {
    
    stop('max.cluster.size must be numerical\n')
    
  }
  
  if(!is.logical(center_size_factors)) {
    
    stop('center_size_factors must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical\n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  assay.list <- list()
  mat <- object@methods[[assay]][[slot]]
  assay.list[[slot]] <- mat
  sce <- SingleCellExperiment::SingleCellExperiment(assay.list)
  
  start_time <- Sys.time()
  
  clusters <- scran::quickCluster(mat)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': quickCluster completed\n')))
    
  }

  sce <- scran::computeSumFactors(sce, clusters=clusters, max.cluster.size=max.cluster.size, assay.type=slot)

  sce <- scuttle::logNormCounts(x = sce, log = F, center.size.factors=center_size_factors, exprs_values=slot)

  .counts <- object@methods[[assay]][[slot]]

  SummarizedExperiment::assay(sce, 'logcounts') <- log2(as_matrix(.counts) + 1)

  .normalised <- SummarizedExperiment::assay(sce, 'logcounts')

  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': normalisation completed\n')))
    
  }
  
  feat.meta <- feature_metadata(assay = .counts, col.prefix = paste0('SCRAN', new.assay.suffix))
  if(!is.null(batch)) {
    
    dec <- scran::modelGeneVar(sce, assay.type='logcounts', block=object@sample_metadata[[batch]])
    
  } else {
    
    dec <- scran::modelGeneVar(sce, assay.type='normcounts')
    
  }
  
  top.hvgs <- scran::getTopHVGs(stats = dec, n = n.genes)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': HVGs identified\n')))
    
  }
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  
  if(!is.null(vars.to.regress)) {
    
    vars.to.regress.df <- as.data.frame(object@sample_metadata[,vars.to.regress])
    colnames(vars.to.regress.df) <- vars.to.regress
    rownames(vars.to.regress.df) <- colnames(object)
    
    vars.to.regress.df <- as.data.frame(vars.to.regress.df[match(rownames(seuobj@meta.data), rownames(vars.to.regress.df)),])
    colnames(vars.to.regress.df) <- vars.to.regress
    rownames(vars.to.regress.df) <- colnames(seuobj)
    
    seuobj@meta.data <- cbind(seuobj@meta.data,vars.to.regress.df)
    
    seuobj@assays$RNA@data <- as(.normalised, 'dgCMatrix')[top.hvgs,]
    seuobj@assays$RNA@var.features <- top.hvgs
    seuobj <- Seurat::ScaleData(object = seuobj, vars.to.regress=vars.to.regress, do.scale=do.scale, do.center=do.center, verbose = verbose, features=top.hvgs, ...)
    
  } else {
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center, verbose=verbose, ...)
    
  }
  
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = .normalised, col.prefix = paste0('SCRAN', new.assay.suffix)))
  
  print(seuobj@assays$RNA@scale.data)
  
  .norm.scaled <- seuobj@assays$RNA@scale.data
  
  if('_' %in% unlist(strsplit(x = new.assay.suffix, split = ''))) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in new.assay.suffix, replacing with - \n')))
      
    }
    
    new.assay.suffix <- sub(pattern = '_', replacement = '-', x = new.assay.suffix)
    
  }
  
  object@methods[[paste0('SCRAN', new.assay.suffix)]] <- new(Class = 'methods',
                                          counts = as(.counts, 'dgCMatrix'), 
                                          normalised = as(.normalised, 'dgCMatrix'), 
                                          norm.scaled = .norm.scaled,
                                          highly.variable.genes = top.hvgs,
                                          feature_metadata = feat.meta)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': Scran normalisation completed \n')))
    
  }
  
  end_time <- Sys.time()
  
  function_time <- end_time - start_time
  
  if(!'normalisation_method' %in% colnames(object@pipelines)) {
    
    object@pipelines <- data.frame(normalisation_method=paste0('SCRAN', new.assay.suffix), normalisation_time=function_time)
    
  } else if (paste0('SCRAN', new.assay.suffix) %in% object@pipelines[,'normalisation_method']) {
    
    object@pipelines[which(object@pipelines[,'normalisation_method']==paste0('SCRAN', new.assay.suffix)),] <- data.frame(normalisation_method=paste0('SCRAN', new.assay.suffix), normalisation_time=function_time)
    
  } else {
    
    object@pipelines <- rbind(object@pipelines, data.frame(normalisation_method=paste0('SCRAN', new.assay.suffix), normalisation_time=function_time))
    
  }
  
  return(object)
}