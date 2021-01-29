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

filter.features.cells <- function(object, 
                                  use.assay, 
                                  max.features=NULL, 
                                  min.features=NULL, 
                                  max.mt=NULL, 
                                  min.cell=NULL, 
                                  omit.zero.features=FALSE) {
  if(omit.zero.features != NULL) {
    rowData(object)$total.cells <- rowSums(assay(object, use.assay) != 0)
    cat(crayon::cyan('Empty features removed\n'))
  } else if (omit.zero.features == NULL) {
    cat(crayon::cyan('Empty features not removed\n'))
  }
  if (max.features != NULL) {
    object <- object[rownames(subset(rowData(object), total.cells > as.numeric(min.cell))),]
    cat(crayon::cyan('Minimum cells removed\n'))
  } else if (omit.zero.features == NULL) {
    cat(crayon::cyan('Minimum cells not removed\n'))
  }
  if (min.features != NULL) {
    object <- object[,rownames(subset(colData(object), total.features < as.numeric(max.features)))]
    cat(crayon::cyan('Maximum features removed\n'))
  } else if (omit.zero.features == NULL) {
    cat(crayon::cyan('Minimum cells not removed\n'))
  }
  if (max.mt != NULL) {
    object <- object[,rownames(subset(colData(object), percent.mt < as.numeric(max.mt)))]
    cat(crayon::cyan('Maximum mitochondrial feature percentage removed\n'))
  } else if (max.mt == NULL) {
    cat(crayon::cyan('Maximum mitochondrial feature percentage not removed\n'))
  }
  cat(crayon::cyan('Filtration complete!\n'))
  return(object)
}
