#' @name splitIBRAP
#' @aliases splitIBRAP
#' 
#' @title Split the IBRAP object up
#'
#' @description This breaks the IBRAP object up by a specified IBRAP object. 
#' 
#' @param object IBRAP S4 class object
#' @param split.by Character. Which column contained within object.sample_metadata should be used to split the object up 
#' 
#' @return A list of the separated IBRAP objects 
#' 
#' @examples list.of.objects <- splitIBRAP(object = object, 
#'                                         split.by = 'original.project')
#'
#' @export

splitIBRAP <- function(object, split.by) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be IBRAP class object \n')
    
  }
  
  if(!is.character(split.by)) {
    
    stop('split.by must be a character string \n')
    
  } else if(is.character(split.by)) {
    
    if(!split.by %in% colnames(object@sample_metadata)) {
      
      stop('split.by is not a column name in object@sample_metadata')
      
    }
    
  }
  
  object.list <- list()
  
  count <- 1
  
  for(x in unique(object@sample_metadata[,split.by])) {
    
    object.list[[count]] <- object[,object@sample_metadata[,split.by]==x]
    
    count <- count + 1
    
  }
  
  return(object.list)
  
}