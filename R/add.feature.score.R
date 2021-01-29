#' Takes a list of features and produces an expression score. 
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param assay Which assay should be utilised
#' @param transform Should the matrix be logtransformed, only do this to raw counts. 
#' @param features A list of features to calculate the score for. 
#' @export
#'

add.feature.score <- function(object, 
                              assay, 
                              transform, 
                              features, 
                              ...) {
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/Homo_sapiens.csv', header = TRUE, sep = ',')
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = assay(object, assay))
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- as.Seurat(object, counts = NULL, data = assay)
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}