#' Add a cell cycle score according to a list of pre-saved list of genes. 
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param assay Which assay should be utilised
#' @param transform Should the matrix be logtransformed, only do this to raw counts. 
#' @export
#'

add.cell.cycle <- function(object, 
                           assay, 
                           transform, ...) {
  
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/Homo_sapiens.csv', header = TRUE, sep = ',')
  cat(crayon::cyan('Cell cycle genes loaded\n'))
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = assay(object, assay))
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Cell cycle scores identified\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- as.Seurat(object, counts = NULL, data = assay)
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Data transformed\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}