#' Here performs Seurat::FindVariableFeatures
#'
#' This method identifieds the most variable genes with a selection of three tests: vst. mvp, and disp
#'
#' @import Seurat
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param nfeatures How many features do you want to retain in as the most variable.
#' @param ... Specify parameters to be passed to the FindVariableFeatures function.
#' @export
#' @examples find_HVGs(object = sce_object, nfeatures = 2000)
#'

find_HVGs <- function(object, nfeatures, feat.to.omit, ...) {
  if(typeof(object) != 'S4') {
    return('Must be an S4 SCE object')
  } else {
    exp.names <- altExpNames(object)
    for(o in exp.names) {
      y <- altExp(object, o)
      print(paste0('calculating highly variable genes for ', o, '.'))
      tmp <- as.Seurat(y, counts = 'counts', data = 'data')
      f <- FindVariableFeatures(object = tmp[!(rownames(tmp) %in% feat.to.omit)], nfeatures = nfeatures, ...)
      tmp <- f@assays$RNA@var.features
      metadata(y)$HVGs <- tmp
      altExp(object, o) <- y
    }

    return(object)
  }
}
