#' Specific feature filtration
#'
#' A defined feature list is omitted from the sce object
#'
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param features Here specify the features you want to remove, i.e. c('MTOR', 'IGLC1')
#' @examples filter_other_features(object = sce_object, features = c('MTOR', 'IGLC1'))
#' @export

filter_other_features <- function(object, features=list()) {
  object <- object[!rownames(object) %in% features, ]
  return(object)
}
