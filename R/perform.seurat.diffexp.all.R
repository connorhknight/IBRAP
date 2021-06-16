#' @name perform.seurat.diffexp.all
#' @aliases perform.seurat.diffexp.all
#' 
#' @title Perform differential expression one cluster vs all
#'
#' @description Performs differential expression on the assay differentiating the supplied identities, some methods enable variable regression, including: LR, negbinom, poisson or MAST. This function compares one cluster vs all others
#' 
#' @param object An IBRAP class object
#' @param assay Character. Which assay within the IBRAP object to access. Default = NULL
#' @param test Character. Which test to use. Can be either: wilcox, bimod, roc, t, negbinom, poisson, LR, MAST, DESeq2. Please refer to Seurat::FindMarkers for more information.
#' @param identity Vector. A vector of cell identifiers to distinguish cells. Default = NULL
#' @param latent.vars Character. String(s) identifying which variables contained within the metadata to regress from the cells 
#' @param ... arguments to pass to Seurat::FindMarkers
#' 
#' @return A dataframe containing differentially expression genes and other information
#'
#' @export
#' 

perform.seurat.diffexp.all <- function(object, 
                                   assay = NULL,
                                   test = 'wilcox', 
                                   identity = NULL,
                                   latent.vars = NULL,
                                   ...) {
  
  if(!is(object, 'IBRAP')) {
    
    cat(crayon::cyan('Object must be IBRAP class \n'))
    return(NULL)
    
  }
  
  if(!is.null(assay)) {
    
    if(!is.character(assay)) {
      
      cat(crayon::cyan('Assay must be character string \n'))
      return(NULL)
      
    }
    
  } else if (is.null(assay)) {
    
    cat(crayon::cyan('Please indicate which assay to access \n'))
    return(NULL)
    
  }
  
  if(!is.character(test)) {
    
    cat(crayon::cyan('Test must be character string \n'))
    return(NULL)
    
  }
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
  seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
  
  if(!is.null(identity)) {
    
    if(!is.vector(identity)) {
      
      cat(crayon::cyan('Identity must be a vector \n'))
      return(NULL)
      
    } else if(length(identity) != ncol(object@methods[[assay]]@counts)) {
      
      cat(crayon::cyan('Identity length does not match the number of cells \n'))
      return(NULL)
      
    }
    
  } else if (is.null(identity)) {
    
    cat(crayon::cyan('Please provide the cell identities \n'))
    return(NULL)
    
  }

  
  if(!is.null(latent.vars)) {
    
    if(!is.character(latent.vars)) {
      
      cat(crayon::cyan('Latent.vars must be character(s)\n'))
      return(NULL)
      
    } else if(!latent.vars %in% names(object@sample_metadata)) {
      
      cat(crayon::cyan('Latent.vars do not exist in object@sample_metadata \n'))
      return(NULL)
      
    }
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, latent.vars = latent.vars, ...)
    
  } else {
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, ...)
    
  }
  
  return(results)
  
}