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
                        new.assay.suffix = '') {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
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
  
  r <- read.table(text = as.character(IBRAP::mart_export$Gene.stable.ID.Gene.name.Gene.start..bp..Gene.end..bp.), sep = ',')
  colnames(r) <- c('geneID', 'geneName', 'start', 'end')
  
  if(is.null(r)) {
    
    stop('cannot find gene lengths\n')
    
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
  
  r$Gene.length <- r$end - r$start
  
  subset <- r[r$geneName %in% rownames(object),]
  
  cat(crayon::cyan(paste0(Sys.time(), ': matrix subsetted\n')))
  
  rownames(subset) <- make.unique(names = as.character(subset$geneName), '.')
  
  cat(crayon::cyan(paste0(Sys.time(), ': rownames added\n')))
  
  meta <- object@methods[[assay]]@feature_metadata[intersect(rownames((object@methods[[assay]]@feature_metadata)), rownames(subset)),]
  
  cat(crayon::cyan(paste0(Sys.time(), ': gene names interesected\n')))
  
  mat <- as.matrix(object@methods[[assay]][[slot]])
  
  mat <- mat[intersect(rownames(mat), rownames(subset)),]
  
  ordered <- subset[match(rownames(mat), rownames(subset)),]
  
  cat(crayon::cyan(paste0(Sys.time(), ': matrices ordered\n')))
  
  cat(crayon::cyan(paste0(Sys.time(), ': calculated counts/feature length\n')))
  
  calc <- sweep(mat, 1, as.numeric(ordered$Gene.length), `/`)
  
  scale.factor <- colSums(calc)/1000000
  
  .counts <- sweep(calc, 2, as.numeric(scale.factor), `/`)
  
  cat(crayon::cyan(paste0(Sys.time(), ': calculations completed\n')))
  
  cat(crayon::cyan(paste0(Sys.time(), ': log2(x+1) transforming\n')))
  .normalised <- log2(.counts+1)
  cat(crayon::cyan(paste0(Sys.time(), ': transformation completed\n')))
  
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
  
  object@sample_metadata <- cbind(object@sample_metadata, cell_metadata(assay = as.matrix(.normalised), col.prefix = new.assay.name))
  
  if('_' %in% unlist(strsplit(x = new.assay.suffix, split = ''))) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in new.assay.suffix, replacing with - \n')))
    
    new.assay.suffix <- sub(pattern = '_', replacement = '-', x = new.assay.suffix)
    
  }
  
  object@methods[[paste0('TPM', new.assay.suffix)]] <- new(Class = 'methods',
                                          counts = Matrix::Matrix(.counts, sparse = T), 
                                          normalised = Matrix::Matrix(.normalised, sparse = T), 
                                          norm.scaled = as.matrix(.norm.scaled),
                                          highly.variable.genes = .highly.variable.genes)
  
  cat(crayon::cyan(paste0(Sys.time(), ': TPM normalisation completed\n')))
  
  return(object)
}