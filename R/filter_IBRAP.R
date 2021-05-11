#' @name filter_IBRAP
#' @aliases filter_IBRAP
#' 
#' @title Filters object according to cell metadata
#'
#' @description Filters cells according to thresholds applied in the columns of object@metadata
#' 
#' @param assay IBRAP class object
#' @param ... 
#' 
#' @usage filter_IBRAP(assay = counts, ... = RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 8)
#'
#' @return Filtered IBRAP object
#'
#' @export

filter_IBRAP <- function(object, ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    cat(crayon::cyan('Object must be class IBRAP'))
    return(NULL)
    
  }
  
  metadata <- object@sample_metadata
  
  metadata <- subset(metadata, ...)
  
  object <- object[,rownames(metadata)]
  
  return(object)
  
}