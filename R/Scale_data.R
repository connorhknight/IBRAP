#' Here scales, centres, and regresses unwanted confounders
#'
#' @import Seurat
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param unwanted.variance Which column metadata you wish to regress, can be listed as c('factor_1', 'factore_2')
#' @param scale TRUE/FALSE whether to scale features
#' @param centre TRUE/FALSE whether to centre transform the matrix to make the mean equal 0
#' @param ... Specify parameters to be passed to the ScaleData function.
#' @export
#' @examples Scale_data(object = sce_object, unwanted.variance = 'percent.mt', scale = TRUE, centre = TRUE)
#'

Scale_data <- function(object, unwanted.variance, scale, centre, ...) {
  exp.names <- altExpNames(object)
  for(o in exp.names) {
    print(paste0('Transforming data for ', o, '.'))
    y <- altExp(object, o)
    tmp <- as.Seurat(y, counts = 'counts', data = 'data')
    tmp <- ScaleData(tmp, vars.to.regress = unwanted.variance, do.scale = scale, do.center = scale, ...)
    assay(y, 'scaled') <- as.matrix(tmp@assays$RNA@scale.data)
    altExp(object, o) <- y
  }
  return(object)
}
