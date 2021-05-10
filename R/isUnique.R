#' @name isUnique
#' @aliases isUnique
#' 
#' @title Are vector items unique
#'
#' @description Identifies if a vector is unique
#' 
#' @param vector
#' 
#' @usage isUnique(vector = gene.names)
#'
#' @return Boolean result whether vector is unique
#'
#' @export

isUnique <- function(vector){
  return(!any(duplicated(vector)))
}