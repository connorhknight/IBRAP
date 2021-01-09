#' Here filters the expression matrix
#'
#' Specific metadata variables can be filtered here. This importantly omits useless data , thus reducing computation time and increasing efficency.
#' For other metadata omission, please perform manually using the outlined example
#'
#' Remember, you can view your available  metadata by specifying colData(sce_object)
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param max.features Maximum n of features that should be present in a cell. This enables use to filter potential doublets
#' @param min.features Minimum n of features that should be present in a cell. This is important in omitting empty droplets
#' @param max.mt Maximum percetage of mitochondrial genes to be retained post-filtration. This helps to evade cells displaying mitochondrial damage during analysis.
#' @param min.cell The minimum number of cells a feature should be present in to be retained.
#' @export
#'

filter_features_cells <- function(object, max.features, min.features, max.mt, min.cell) {
  rowData(object)$total.cells <- rowSums(assay(object, 'counts') != 0)
  object <- object[rownames(subset(rowData(object), total.cells > as.numeric(min.cell))),]
  object <- object[,rownames(subset(colData(object), total.counts > as.numeric(min.features)))]
  object <- object[,rownames(subset(colData(object), total.features < as.numeric(max.features)))]
  object <- object[,rownames(subset(colData(object), percent.mt < as.numeric(max.mt)))]
  return(object)
}
