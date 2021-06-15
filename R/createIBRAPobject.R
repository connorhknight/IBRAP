#' @name createIBRAPobject
#' @aliases createIBRAPobject
#' 
#' @title Create an IBRAP class object 
#'
#' @description Creates and produces project metadata into an IBRAP S4 class object
#' 
#' @param counts Counts matrix
#' @param original.project Character string naming the project
#' @param method.name Character string naming the original method-assay. Default = 'RAW' `(Recommended)`
#' @param meta.data data.frame of extra metadata to append to generated dataframe. Warning: Must be in the same order and colnames
#' @param min.cells Numerical value of the minimum number of cells a gene should be present in
#' @param min.features Numerical value minimum features that should be present in a cell
#' 
#' @usage createIBRAPobject(counts = counts, original.project = 'project_1', method.name = 'RAW', meta.data = df, min.cells = 3, min.features = 200)
#' 
#' @return IBRAP S4 class object containing raw counts and metadata
#'
#' @export

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