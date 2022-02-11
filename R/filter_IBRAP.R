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
    
    stop('Object must be class IBRAP')
    
  }
  
  metadata <- object@sample_metadata
  
  pre_filt <- nrow(metadata)
  
  metadata <- subset(metadata, ...)
  
  post_filt <- nrow(metadata)
  
  filt_loss <- pre_filt - post_filt
  
  object <- object[,rownames(metadata)]
  
  cat(crayon::cyan(paste0(Sys.time(), ': a total of ', 
                          filt_loss, ' cells were filtered from ', 
                          pre_filt, ', ', post_filt, ' cells remain\n')))
  
  return(object)
  
}