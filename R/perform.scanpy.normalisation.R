#' @name perform.scanpy
#' @aliases perform.scanpy
#' 
#' @title Performs Scanpy normalisation, hvg selection, scaling and variance stabilisation and regression. 
#'
#' @description A new method-assay is produced. Raw counts are normalised and HVGs identified using Scanpy
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param new.assay.suffix Character. What should be added as a suffix to 'SCANPY
#' @param target_sum Numerical. What should the data be scaled to. Default = 1e6
#' @param exclude_highly_expressed Boolean. Should highly expressed genes be excluded. Default = FALSE
#' @param max_fraction Numerical. If exclude_highly_expressed=True, consider cells as highly expressed that have more counts than max_fraction of the original total counts in at least one cell. Default = 0.05
#' @param key_added Character. What should the column name be that contains cell scaling factors. Default = 'scanpy_norm_factor'
#' @param n_top_genes Numerical. How many HVGs should be identified. Default = 1500
#' @param max_mean Numerical. If n_top_genes is NULL, this is the maximum mean to determine HVGs. Default = 6
#' @param min_mean Numerical. If n_top_genes is NULL, this is the minimum mean to determine HVGs. Default = 0.0125
#' @param min_disp Numerical. If n_top_genes is NULL, The minimum dispersion that should be presented in a gene for it to be considered highly varaible. Default = 0.5
#' @param span Numerical. The fraction of cells that should be subset for the LOESS fit model. Default = 0.3
#' @param n_bins Numerical. Number of bins to produce when determining HVGs
#' @param flavour Character. Choosing which HVG selection method to use when, options: 'seurat', 'cell_ranger', 'seurat_v3'. Default = 'seurat'
#' @param batch_key Character. Which column in the metadata identifies the batches of the cells. Default = NULL
#' @param vars.to.regress Character. A single or multiple columns of information in the metadata that should be regressed from the dataset. Default = NULL
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#' 
#' @examples 
#' 
#' object <- perform.scanpy(object = object, 
#'                          vars.to.regress = 'RAW_total.counts')
#'
#' @export

perform.scanpy <- function(object, 
                           assay='RAW', 
                           slot='counts', 
                           new.assay.suffix='', 
                           target_sum = 1e4, 
                           exclude_highly_expressed = FALSE,  
                           max_fraction = 0.05, 
                           key_added = 'scanpy_norm_factor',
                           log1 = FALSE,
                           
                           n_top_genes = 1500, 
                           max_mean = 6, 
                           min_mean = 0.0125, 
                           min_disp = 0.5, 
                           span = 0.3, 
                           n_bins = 20, 
                           flavor = 'seurat', 
                           batch_key = NULL,
                           vars.to.regress=NULL
) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Object must be of class IBRAP\n')
    
  }
  if(!is.character(assay)) {
    
    stop('Assay must be a character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('assay does not exist\n')
    
  }
  
  if(!is.character(slot)) {
    
    stop('Slot must be a character string\n')
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    stop('slot does not exist\n')
    
  }
  
  if(!is.character(new.assay.suffix)) {
    
    stop('new.assay.suffix must be character string \n')
    
  }
  
  if(!is.numeric(target_sum)) {
    
    stop('target_sum must be numerical \n')
    
  }
  
  if(!is.logical(exclude_highly_expressed)) {
    
    stop('exclude_highly_expressed must be logical\n')
    
  }
  
  if(!is.numeric(max_fraction)) {
    
    stop('max_fraction must be numerical \n')
    
  }
  
  if(!is.character(key_added)) {
    
    stop('key_added must be character string \n')
    
  }
  
  if(!is.numeric(n_top_genes)) {
    
    stop('n_top_genes must be numerical \n')
    
  }
  
  if(!is.numeric(max_mean)) {
    
    stop('max_mean must be numerical \n')
    
  }
  
  if(!is.numeric(min_mean)) {
    
    stop('min_mean must be numerical \n')
    
  }
  
  if(!is.numeric(min_disp)) {
    
    stop('min_disp must be numerical \n')
    
  }
  
  if(!is.numeric(span)) {
    
    stop('span must be numerical \n')
    
  }
  
  if(!is.character(flavor)) {
    
    stop('flavor must be character string \n')
    
  }
  
  if(!is.null(batch_key)) {
    
    if(!is.character(batch_key)) {
      
      stop('batch_key must be character string\n')
      
    }
    
  }
  
  if(!is.null(batch_key)) {
    
    if(!is.numeric(batch_key)) {
      
      stop('batch_key must be character string\n')
      
    }
    
  }
  
  if(!is.null(vars.to.regress)) {
    
    if(!is.character(vars.to.regress)) {
      
      stop('vars.to.regress must be character string\n')
      
    }
    
  }
  
  sc <- reticulate::import('scanpy')
  scobj <- sc$AnnData(X = t(as.matrix(object@methods[[assay]][[slot]])))
  scobj$obs_names <- as.factor(colnames(object@methods[[assay]][[slot]]))
  scobj$var_names <- as.factor(rownames(object@methods[[assay]][[slot]]))
  
  if(length(names(object@sample_metadata)) >= 1) {
    scobj$obs <- object@sample_metadata
  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': normalising counts\n')))
  
  if(!is.null(target_sum) & !is.null(key_added)) {
    sc$pp$normalize_total(adata = scobj, target_sum = as.integer(target_sum), 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction), key_added = as.character(key_added))
  } else if (!is.null(target_sum)) {
    sc$pp$normalize_total(adata = scobj, target_sum = as.integer(target_sum), 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction))
  } else if (!is.null(key_added)) {
    sc$pp$normalize_total(adata = scobj,  target_sum=target_sum,
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction))
  } else {
    sc$pp$normalize_total(adata = scobj, target_sum=target_sum)
  }
  
  .counts <- t(scobj$X)
  rownames(.counts) <- rownames(object)
  colnames(.counts) <- colnames(object)
  
  feat.metadata <- feature_metadata(assay = .counts, col.prefix = 'SCANPY')
  
  cat(crayon::cyan(paste0(Sys.time(), ': log transforming data\n')))
  
  if(isTRUE(log1)) {
    
    sc$pp$log1p(scobj)
    
  } else if(isFALSE(log1)) {
    
    scobj$X <- log2(scobj$X+1)
    
  }
  
  .normalised <- t(scobj$X)
  rownames(.normalised) <- rownames(object@methods$RAW@counts)
  colnames(.normalised) <- colnames(object@methods$RAW@counts)
  
  cat(crayon::cyan(paste0(Sys.time(), ': computing highly variable genes\n')))
  
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
  
  .highly.variable.genes <- rownames(object@methods$RAW@counts)[scobj$var[['highly_variable']]]
  
  scobj2 <- sc$AnnData(X = t(.normalised[scobj$var$highly_variable,]))
  
  pd <- reticulate::import('pandas')
  
  scobj2$var_names <- as.factor(rownames(object@methods$RAW@counts)[scobj$var$highly_variable])
  scobj2$obs_names <- as.factor(colnames(object@methods$RAW@counts))
  
  if(length(vars.to.regress) > 1) {
    
    scobj2$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata[,vars.to.regress]))
    
  } else {
    
    scobj2$obs[,vars.to.regress] <- pd$DataFrame(data = as.data.frame(object@sample_metadata[,vars.to.regress]))
    
  }
  

  cat(crayon::cyan(paste0(Sys.time(), ': regressing covaraites \n')))
  
  sc$pp$regress_out(adata = scobj2, keys = vars.to.regress)
  
  cat(crayon::cyan(paste0(Sys.time(), ': scaling data \n')))
  
  sc$pp$scale(scobj2)

  .norm.scaled <- t(scobj2$X)

  colnames(.norm.scaled) <- colnames(object@methods$RAW@counts)

  rownames(.norm.scaled) <- .highly.variable.genes
  
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = as.matrix(.normalised), col.prefix = paste0('SCANPY', new.assay.suffix)))
  
  if('_' %in% unlist(strsplit(x = new.assay.suffix, split = ''))) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in new.assay.suffix, replacing with - \n')))
    
    new.assay.suffix <- sub(pattern = '_', replacement = '-', x = new.assay.suffix)
    
  }

  object@methods[[paste0('SCANPY', new.assay.suffix)]] <- new(Class = 'methods',
                                                              counts = as(.counts, 'dgCMatrix'), 
                                                              normalised = as(.normalised, 'dgCMatrix'), 
                                                              norm.scaled = as.matrix(.norm.scaled),
                                                              highly.variable.genes = .highly.variable.genes,
                                                              feature_metadata = feat.metadata)
  
  cat(crayon::cyan(paste0(Sys.time(), ': Scanpy normalisation completed \n')))
  
  return(object)
}