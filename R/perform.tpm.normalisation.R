#' @name perform.tpm
#' @aliases perform.tpm
#' 
#' @title Performs TPM normalisation
#'
#' @description Performs TPM normalisation, scran hvg selection, scaling and variance stabilisation and regression. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param n.genes Numerical. Top number of genes to retain when finding HVGs. Default = 1500
#' @param do.scale Boolean. Whether to scale the features variance. Default = TRUE
#' @param do.centre Boolean. Whether to centre features to zero. Default = TRUE
#' @param vars.to.regress Character. Which column in the metadata should be regressed. Default = NULL
#' @param new.assay.suffix Character. What should the new assay be called. Default = 'SCRAN'
#' @param biomart.dataset Character. Which biomart dataset should be used, this normally corresponds with the species in question, default = 'hsapiens_gene_ensembl'. Check available datasets by performing the following. ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl"), then biomaRt::listDatasets(ensembl)
#' @gene.lengths DataFrame. A dataframe containing two columns, external_gene_name and transcript_length.
#' @param verbose Logical Should function messages be printed?
#' @param ... Arguments to pass to Seurat::ScaleData
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#' 
#' @examples 
#' 
#' object <- perform.scanpy(object = object, 
#'                          vars.to.regress = 'RAW_total.counts', do.scale = T)
#'
#' @export

perform.tpm <- function(object, 
                        assay = 'RAW', 
                        slot = 'counts',
                        n.genes = 1500,
                        do.scale = FALSE,
                        do.center = TRUE,
                        vars.to.regress = NULL,
                        new.assay.suffix = '',
                        biomart.dataset='hsapiens_gene_ensembl',
                        gene.lengths=NULL,
                        verbose = FALSE,
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
  
  if(!is.numeric(n.genes)) {
    
    stop('n.genes must be numerical\n')
    
  }
  
  if(!is.logical(do.scale)) {
    
    stop('do.scale must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.logical(do.center)) {
    
    stop('do.center must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.null(vars.to.regress)) {
    
    if(!is.character(vars.to.regress)) {
      
      stop('vars.to.regress must be character string\n')
      
    }
    
  }
  
  if(!is.character(new.assay.suffix)) {
    
    stop('new.assay.suffix must be a character string\n')
    
  }
  
  if('_' %in% unlist(strsplit(x = new.assay.suffix, split = ''))) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in new.assay.suffix, replacing with - \n')))
      
    }
    
    new.assay.suffix <- sub(pattern = '_', replacement = '-', x = new.assay.suffix)
    
  }
  
  if(!is.null(gene.lengths)) {
    
    if(!is.data.frame(gene.lengths)) {
      
      stop('gene.lengths must be a dataframe\n')
      
    }
    
    if(!c('external_gene_name','transcript_length') %in% colnames(gene.lengths)) {
      
      stop('gene.lengths must contain column names external_gene_name and transcript_length\n')
      
    }
    
  }
  
  if(is.null(gene.lengths)) {
    
    ensembl <- biomaRt::useEnsembl(biomart = 'genes', dataset = biomart.dataset)
    
    query <- biomaRt::getBM(attributes = c('external_gene_name','transcript_start','transcript_end'), 
                            filters = 'external_gene_name', mart = ensembl, values = rownames(object))
    
    query <- query[order(query$external_gene_name, decreasing = T),]
    
    query <- query[!duplicated(query$external_gene_name),]
    
    query$transcript_length <- query$transcript_end - query$transcript_start
    
    gene.lengths <- query[,c('external_gene_name','transcript_length')]
    
  }
  
  start_time <- Sys.time()
  
  gene.length.subset <- gene.lengths[gene.lengths$external_gene_name %in% rownames(object),]
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': matrix subsetted\n')))
    
  }
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': rownames added\n')))
    
  }
  
  meta <- object@methods[[assay]]@feature_metadata[intersect(rownames((object@methods[[assay]]@feature_metadata)), gene.length.subset$external_gene_name),]
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': gene names interesected\n')))
    
  }
  
  mat <- as.matrix(object@methods[[assay]][[slot]])
  
  mat <- mat[intersect(rownames(mat), gene.length.subset$external_gene_name),]
  
  ordered <- gene.length.subset[match(rownames(mat), gene.length.subset$external_gene_name),]
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': matrices ordered\n')))
    
    cat(crayon::cyan(paste0(Sys.time(), ': calculated counts/feature length\n')))
    
  }
  
  calc <- sweep(mat, 1, as.numeric(ordered$transcript_length), `/`)
  
  scale.factor <- colSums(calc)/1000000
  
  .counts <- sweep(calc, 2, as.numeric(scale.factor), `/`)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': calculations completed\n')))
    
  }
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': log2(x+1) transforming\n')))
    
  }
  
  .normalised <- log2(.counts+1)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': transformation completed\n')))
    
  }
  
  dec <- scran::modelGeneVar(x = .normalised)
  .highly.variable.genes <- scran::getTopHVGs(stats = dec, n=n.genes)
  seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = mat))
  seuobj@assays$RNA@data <- .normalised[.highly.variable.genes,]
  
  if(!is.null(vars.to.regress)) {
    
    vars.to.regress.df <- as.data.frame(object@sample_metadata[,vars.to.regress])
    colnames(vars.to.regress.df) <- vars.to.regress
    rownames(vars.to.regress.df) <- colnames(object)
    
    vars.to.regress.df <- vars.to.regress.df[match(rownames(seuobj@meta.data), 
                                                   rownames(vars.to.regress.df)),]
    seuobj@meta.data <- cbind(seuobj@meta.data,vars.to.regress.df)
    
    colnames(seuobj@meta.data) <- c(names(seuobj@meta.data)[1:sum(length(names(seuobj@meta.data))-length(vars.to.regress))], vars.to.regress)
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center,vars.to.regress=vars.to.regress, verbose = verbose, ...)
    
  } else {
    
    seuobj <- Seurat::ScaleData(object = seuobj, do.scale=do.scale, do.center=do.center, verbose = verbose, ...)
    
  }
  
  .norm.scaled <- seuobj@assays$RNA@scale.data
  
  feat.meta <- feature_metadata(assay = .counts, col.prefix = paste0('SCRAN', new.assay.suffix))
  
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = as.matrix(.normalised), col.prefix = paste0('TPM', new.assay.suffix)))
  
  object@methods[[paste0('TPM', new.assay.suffix)]] <- new(Class = 'methods',
                                                           counts = Matrix::Matrix(.counts, sparse = T), 
                                                           normalised = Matrix::Matrix(.normalised, sparse = T), 
                                                           norm.scaled = as.matrix(.norm.scaled),
                                                           highly.variable.genes = .highly.variable.genes,
                                                           feature_metadata = feat.meta)
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': TPM normalisation completed\n')))
    
  }
  
  end_time <- Sys.time()
  
  function_time <- end_time - start_time
  
  if(!'normalisation_method' %in% colnames(object@pipelines)) {
    
    object@pipelines <- data.frame(normalisation_method=paste0('TPM', new.assay.suffix), normalisation_time=function_time)
    
  } else if (paste0('TPM', new.assay.suffix) %in% object@pipelines[,'normalisation_method']) {
    
    object@pipelines[which(object@pipelines[,'normalisation_method']==paste0('TPM', new.assay.suffix)),] <- data.frame(normalisation_method=paste0('TPM', new.assay.suffix), normalisation_time=function_time)
    
  } else {
    
    object@pipelines <- rbind(object@pipelines, data.frame(normalisation_method=paste0('TPM', new.assay.suffix), normalisation_time=function_time))
    
  }
  
  return(object)
}