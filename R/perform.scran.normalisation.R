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
#' @param do.centre Boolean. Whether to centre features to zero. Default = TRUE
#' @param new.assay.name Character. What should the new assay be called. Default = 'SCRAN'
#' @param n.genes Numerical. Top number of genes to retain when finding HVGs. Default = 1500
#' @param max.cluster.size Numerical. When performing quickCluster, what is the maximum size the clusters can be. Default = 1000
#' @param center_size_factors Boolean Should size factor variance be centred. Default = TRUE
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#'
#' @export

perform.scran <- function(object, 
                          assay = 'RAW',
                          slot = 'counts',
                          batch=NULL,
                          vars.to.regress=NULL,
                          do.scale=FALSE,
                          do.center=TRUE,
                          new.assay.name = 'SCRAN',
                          n.genes=1500,
                          max.cluster.size = 1000,
                          center_size_factors=TRUE) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Object must be of class IBRAP\n')
    
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
      
      cstop('vars.to.regress must be a character string(s) of column name contained in sample_metadata\n')
      
    }
    
  }
  
  if(!is.logical(do.scale)) {
    
    stop('do.scale must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.logical(do.center)) {
    
    stop('do.center must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.character(new.assay.name)) {
    
    stop('new.assay.name must be character string\n')
    
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
  
  assay.list <- list()
  mat <- object@methods[[assay]][[slot]]
  assay.list[[slot]] <- mat
  sce <- SingleCellExperiment::SingleCellExperiment(assay.list)
  
  clusters <- scran::quickCluster(mat)
  
  cat(crayon::cyan(paste0(Sys.time(), ': quickCluster completed\n')))
  sce <- scran::computeSumFactors(sce, clusters=clusters, max.cluster.size=max.cluster.size, assay.type=slot)
  sce <- scuttle::logNormCounts(x = sce, log = F, center.size.factors=center_size_factors, exprs_values=slot)
  .counts <- SummarizedExperiment::assay(sce, 'normcounts')
  SummarizedExperiment::assay(sce, 'logcounts') <- log2(.counts + 1)
  .normalised <- SummarizedExperiment::assay(sce, 'logcounts')
  cat(crayon::cyan(paste0(Sys.time(), ': normalisation completed\n')))
  feat.meta <- feature_metadata(assay = .counts, col.prefix = new.assay.name)
  if(!is.null(batch)) {
    
    dec <- scran::modelGeneVar(sce, assay.type='logcounts', block=object@sample_metadata[[batch]])
    
  } else {
    
    dec <- scran::modelGeneVar(sce, assay.type='normcounts')
    
  }
  
  top.hvgs <- scran::getTopHVGs(stats = dec, n = n.genes)
  
  cat(crayon::cyan(paste0(Sys.time(), ': HVGs identified\n')))
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  
  if(!is.null(vars.to.regress)) {
    
    vars.to.regress.df <- as.data.frame(object@sample_metadata[,vars.to.regress])
    colnames(vars.to.regress.df) <- vars.to.regress
    rownames(vars.to.regress.df) <- colnames(object)
    
    vars.to.regress.df <- vars.to.regress.df[match(rownames(seuobj@meta.data), rownames(vars.to.regress.df)),]
    seuobj@meta.data <- cbind(seuobj@meta.data,vars.to.regress.df)
    
    colnames(seuobj@meta.data) <- c(names(seuobj@meta.data)[1:sum(length(names(seuobj@meta.data))-length(vars.to.regress))], vars.to.regress)
    
    seuobj@assays$RNA@data <- as(.normalised, 'dgCMatrix')[top.hvgs,]
    seuobj <- Seurat::ScaleData(object = seuobj, vars.to.regress=vars.to.regress, do.scale=do.scale, do.center=do.center)
    
  } else {
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center)
    
  }
  
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = as.matrix(.counts), col.prefix = new.assay.name))
  
  .norm.scaled <- seuobj@assays$RNA@scale.data
  
  object@methods[[new.assay.name]] <- new(Class = 'methods',
                                          counts = as(.counts, 'dgCMatrix'), 
                                          normalised = as(.normalised, 'dgCMatrix'), 
                                          norm.scaled = as.matrix(.norm.scaled),
                                          highly.variable.genes = top.hvgs,
                                          feature_metadata = feat.meta)
  cat(crayon::cyan(paste0(Sys.time(), ': Scran normalisation completed \n')))
  return(object)
}