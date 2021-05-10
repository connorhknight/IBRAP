#' @name .convert_subset_index
#' @aliases .convert_subset_index
#' 
#' @param x index
#' @param y converts index into this varaible, must be the same length as supplied indexing 
#' 
#' @export

.convert_subset_index <- function(x, 
                                  y) {
  
  return(y[x])
  
}