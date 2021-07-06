remove.hvgs <- function(object, assay, hvgs.omit) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character(s)')
    
  } else {
    
    for(x in assay) {
      
      if(!x %in% names(object@methods)) {
        
        stop('assay is not contained within object@methods \n')
        
      }
      
    }
    
  }
  
  if(!is.list(hvgs.omit)) {
    
    stop('hvgs.omit must be in list format \n')
    
  } else {
    
    if(length(hvgs.omit) != length(assay)) {
      
      stop('an equal number of hvgs.omit items and assay must be supplied \n')
      
    }
    
    for(x in 1:length(hvgs.omit)) {
      
      tmp <- hvgs.omit[[x]]
      
      for(x in tmp) {
        
        if(!is.character(x)) {
          
          stop('hvgs must be supplied as character strings \n')
          
        }
        
      }
      
    }
    
  }
  
  count <- 1
  
  for(x in assay) {
    
    genes <- as.character(object@methods[[x]]@highly.variable.genes)
    
    genes <- genes[!genes %in% as.character(hvgs.omit[[count]])]
    
    object@methods[[x]]@highly.variable.genes <- genes
    
    count <- count + 1  
      
  }
  
  return(object)
  
}