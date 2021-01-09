#' Here scales, centres, and regresses unwanted confounders
#'
#' @import harmony
#' @import SingleCellExperiment
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param unwanted.variance Which column metadata you wish to regress, can be listed as c('factor_1', 'factor_2')
#' @param scale TRUE/FALSE whether to scale features
#' @param centre TRUE/FALSE whether to centre transform the matrix to make the mean equal 0
#' @param ... Specify parameters to be passed to the ScaleData function.
#' @export
#' @examples Scale_data(object = sce_object, unwanted.variance = 'percent.mt', scale \
#' 

apply.harmony <- function(object, norm.methods, reduction, n.components, vars.to.regress, , ...) {
  
  for(t in norm.methods) {
    print(paste0('harmony processing:', as.character(t)
    temp <- altExp(object, as.character(t))
    red <- reducedDim(object, as.character(reduction))
    met <- colData(temp)
    embeddings <- HarmonyMatrix(data_mat = red, meta_data = met, return_object = FALSE, npcs = as.numeric(n.components), 
                  vars_use = vars.to.regress, do_pca = FALSE, verbose = TRUE, ...)
    reducedDim(temp, 'harmony') <- as.data.frame(embeddings)
    altExp(object, as.character(t)) <- temp
  }
  
}