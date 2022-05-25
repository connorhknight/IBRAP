#' @name perform.sct
#' @aliases perform.sct
#' 
#' @title Performs SCTransform
#'
#' @description Performs SCTransform normalisation, hvg selection, scaling and variance stabilisation and regression.  
#' 
#' @param object IBRAP S4 class object
#' @param assay A character string containing indicating which assay to use
#' @param slot String indicating which slot within the assay should be sourced
#' @param new.assay.suffix Character. What should be added as a suffix for SCT
#' @param do.scale Whether to scale residuals to have unit variance; default is FALSE
#' @param do.center Whether to center residuals to have mean zero; default is TRUE
#' @param vars.to.regress Character. Which data from `object@sample_metadata` should be regressed from the dataset.
#' @param n.genes Numerical value of how many highly variable genes should be retained. Default = 1500
#' @param min_cells Numerical value of minimum cells required for a gene to not be filtered. Default = 3
#' @param verbose Logical Should function messages be printed?
#' @param seed Numerical What seed should be set. Default = 1234
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#' 
#' @examples 
#' 
#' object <- perform.sct(object = object, 
#'                       assay = 'RAW', 
#'                       slot = 'counts')
#'
#' @export

perform.sct <- function(object, 
                        assay='RAW',
                        slot='counts',
                        new.assay.suffix = '',
                        verbose = FALSE,
                        seed=1234,
                        ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose must be logical, TRUE/FALSE\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be a character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('assay does not exist\n')
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be a character string\n')
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    stop('slot does not exist\n')
    
  }
  
  if(!is.character(new.assay.suffix)) {
    
    stop('new.assay.suffix must be character string\n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical\n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': converting to Seurat object\n')))
    
  }
  
  seuratobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = as_matrix(object@methods[[assay]][[slot]]), project = 'NA'))
  seuratobj@meta.data <- cbind(seuratobj@meta.data, object@sample_metadata)
  
  start_time <- Sys.time()
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': initiating SCTransform\n')))
    
  }
  
  seuratobj <- suppressWarnings(Seurat::SCTransform(object = seuratobj, verbose = verbose, seed.use = NULL, ...))
  print('.')
  .highly.variable.genes <- as.character(seuratobj@assays$SCT@var.features)
  print('.')
  .counts <- as(object = as_matrix(seuratobj@assays$SCT@counts), Class = 'dgCMatrix')
  print('.')
  .normalised <- as(as_matrix(seuratobj@assays$SCT@data), Class = 'dgCMatrix')
  print('.')
  .norm.scaled <- seuratobj@assays$SCT@scale.data
  print('.')
  feat.meta <- feature_metadata(assay = as_matrix(.counts), col.prefix = paste0('SCT', new.assay.suffix))
  print('.')
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = as_matrix(.normalised), col.prefix = paste0('SCT', new.assay.suffix)))
  print('.')
  if('_' %in% unlist(strsplit(x = new.assay.suffix, split = ''))) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in new.assay.suffix, replacing with - \n')))
      
    }
    
    new.assay.suffix <- sub(pattern = '_', replacement = '-', x = new.assay.suffix)
    print('.')
  }
  
  object@methods[[paste0('SCT', new.assay.suffix)]] <- new(Class = 'methods',
                                          counts = .counts, 
                                          normalised = .normalised, 
                                          norm.scaled = .norm.scaled,
                                          highly.variable.genes = .highly.variable.genes,
                                          feature_metadata = feat.meta)
  print('.')
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': SCT normalisation completed\n')))
    
  }
  
  end_time <- Sys.time()
  
  function_time <- end_time - start_time
  
  if(!'normalisation_method' %in% colnames(object@pipelines)) {
    
    object@pipelines <- data.frame(normalisation_method=paste0('SCT', new.assay.suffix), normalisation_time=function_time)
    
  } else if (paste0('SCT', new.assay.suffix) %in% object@pipelines[,'normalisation_method']) {
    
    object@pipelines[which(object@pipelines[,'normalisation_method']==paste0('SCT', new.assay.suffix)),] <- data.frame(normalisation_method=paste0('SCT', new.assay.suffix), normalisation_time=function_time)
    
  } else {
    
    object@pipelines <- rbind(object@pipelines, data.frame(normalisation_method=paste0('SCT', new.assay.suffix), normalisation_time=function_time))
    
  }
  
  return(object)
  
}
