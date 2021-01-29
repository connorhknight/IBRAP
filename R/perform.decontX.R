#' Performs celda
#'
#' @param object A SingleCellExperiment
#' @param z Numeric or character vector. Cell cluster labels. If NULL, PCA will be used to reduce the dimensionality of the dataset initially, 'umap' from the 'uwot' package will be used to further reduce the dataset to 2 dimenions and the 'dbscan' function from the 'dbscan' package will be used to identify clusters of broad cell types. Default NULL.
#' @param batch Numeric or character vector. Batch labels for cells. If batch labels are supplied, DecontX is run on cells from each batch separately. Cells run in different channels or assays should be considered different batches. Default NULL.
#' @param ... Variables to pass on to celda::decontX
#' @keywords
#' @export
#' @examples Read10X_output(directory = 'Users/user/10x_output')
#'

perform.decontX <- function(object,
                            z = NULL,
                            batch = NULL,
                            verbose = TRUE,
                            ...) {
  
  d <- celda::decontX(x = object,
                      z = z,
                      batch = colData(object)[,batch],
                      verbose = verbose,
                      ...)
  
  cat(crayon::cyan('Decontamination comlpleted\n'))
  
  cat(crayon::cyan(paste0(formatC(sum(d$decontX_contamination)/length(d$decontX_contamination), 
                                  digits = 2), '% average contamination\n')))
  
  clean.matrix <- assay(d, 'decontXcounts')
  cat(crayon::cyan('Matrix isolated\n'))
  clean.matrix <- round(clean.matrix)
  zero.samples <- colSums(clean.matrix) > 0
  object <- object[,zero.samples]
  clean.matrix <- clean.matrix[,zero.samples]
  cat(crayon::cyan('converted to integer\n'))
  assay(object, 'decontXcounts') <- clean.matrix
  cat(crayon::cyan('Added matrix\n'))
  
  object <- metadata.generator(object = object, 
                               assay = 'decontXcounts', 
                               prefix = 'decontaminated_')
  
  cat(crayon::cyan('Finished\n'))
  return(object)
}