#' @name perform.tpm.normalisation
#' @aliases perform.tpm.normalisation
#' 
#' @title Performs TPM normalisation
#'
#' @description A new method-assay is produced. Raw counts are normalised and HVGs identified using TPM and Scran, respectively.  
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param n.genes Numerical. Top number of genes to retain when finding HVGs. Default = 1500
#' @param do.scale Boolean. Whether to scale the features variance. Default = TRUE
#' @param do.centre Boolean. Whether to centre features to zero. Default = TRUE
#' @param vars.to.regress Character. Which column in the metadata should be regressed. Default = NULL
#' @param new.assay.name Character. What should the new assay be called. Default = 'SCRAN'
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#'
#' @export

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
  
  r <- utils::read.csv(system.file("data", "mart_export.csv", package = "IBRAP"), header = TRUE, sep = ',')
  
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