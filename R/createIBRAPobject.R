#' @name createIBRAPobject
#' @aliases createIBRAPobject
#' 
#' @title Create an IBRAP class object 
#'
#' @description Creates and produces project metadata into an IBRAP S4 class object
#' 
#' @param counts Counts matrix
#' @param original.project Character string naming the project
#' @param meta.data data.frame of extra metadata to append to generated dataframe. Warning: Must be in the same order and colnames
#' @param min.cells Numerical value of the minimum number of cells a gene should be present in
#' @param min.features Numerical value minimum features that should be present in a cell
#' @param verbose Logical. Should function information be printed?
#' 
#' @usage createIBRAPobject(counts = counts, original.project = 'project_1', meta.data = df, min.cells = 3, min.features = 200)
#' 
#' @return IBRAP S4 class object containing raw counts and metadata
#' 
#' @examples object <- createIBRAPobject(counts = counts,
#'                                       meta.data = metadata_df,
#'                                       original.project = 'bmmc',
#'                                       min.cells = 3,
#'                                       min.features = 200)
#'
#' @export

createIBRAPobject <- function(counts, 
                              original.project, 
                              add.suffix = FALSE,
                              meta.data = NULL,
                              min.cells=NULL,
                              min.features=NULL,
                              verbose=FALSE) {
  
  if(!is.character(original.project)) {
    
    stop('original.project must be a character string \n')
    
  }
  
  
  
  if(!is(object = as(object = counts, Class = 'dgCMatrix'), class2 = 'dgCMatrix')) {
    
    stop('counts could not be converted into a dgCMatrix\n')
    
  }
 
  if(!is.null(meta.data)) {
    
    if(!is.data.frame(meta.data)) {
      
      stop('meta.data must be of class data.frame')
      
    }
    
    if(FALSE %in% (rownames(meta.data) %in% colnames(counts))) {
      
      stop('meta.data rownames must be the same as counts colnames \n')
      
    } else {
      
      if(isTRUE(add.suffix)) {
        
        rownames(meta.data) <- paste0(original.project, '_', rownames(meta.data))
        
      }

    }
    
  }
  
  if(isTRUE(add.suffix)) {
    
    colnames(counts) <- paste0(original.project, '_', colnames(counts))
    
  }
  
  if(!is.null(min.features)) {
    
    if(!is.numeric(min.features)) {
      
      stop('min.features must be numerical \n')
      
    }
    
    nfeatures <- Matrix::colSums(x = counts > 0)
    
    counts <- counts[, which(x = nfeatures >= min.features)]
    
  }
  
  if(!is.null(min.cells)) { 
    
    if(!is.numeric(min.cells)) {
      
      stop('min.cells must be numerical \n')
      
    }
    
    num.cells <- Matrix::rowSums(x = counts > 0)
    
    counts <- counts[which(x = num.cells >= min.cells), ]
    
  }
  
  if(isTRUE(any(stringr::str_detect(string = rownames(counts), pattern='_'))) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': cannot have _ in gene names, replacing with - ')))
    
    rownames(counts) <- gsub(pattern='_', replacement='-', rownames(seuobj@assays$RNA@data))
    
  }

  meta <- as.data.frame(replicate(n = length(colnames(counts)), expr = original.project))
  
  colnames(meta) <- as.character('original.project')
  
  meta.2 <- cell_metadata(assay = counts, col.prefix = 'RAW')
  
  for(f in colnames(meta.2)) {
    
    meta[,f] <- as.numeric(meta.2[,f])
    
  }
  
  rownames(meta) <- colnames(counts)

  f.metadata <- feature_metadata(assay = counts, col.prefix = 'RAW')
  
  if(!is.null(meta.data)) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': concatenating metadata \n')))
      
    }
    
    l1 <- colnames(meta)
    
    l2 <- colnames(meta.data)
    
    if(isFALSE(isUnique(c(l1,l2)))) {
      
      stop('Column names from meta.data cannot be named: original.project, counts_total.counts or counts_total.features \n')
      
    }
    
    meta <- meta[match(rownames(meta.data), rownames(meta)),]
    
    meta <- cbind(meta, meta.data)
    
    meta <- meta[match(colnames(counts), rownames(meta)),]

  }
  
  f.metadata <- f.metadata[match(rownames(counts), rownames(f.metadata)),]

  ##########################################################
  
  first.method <- new('methods', 
                      counts = as(object = counts, Class = 'dgCMatrix'),
                      normalised = as(matrix(nrow = 0, ncol = 0), 'dgCMatrix'),
                      feature_metadata = f.metadata)
  
  methods <- list()
  
  methods[[as.character('RAW')]] <- first.method
  
  IBRAP.obj <- new(Class = 'IBRAP', 
                   methods = methods, 
                   sample_metadata = meta)
  
  return(IBRAP.obj)
  
}