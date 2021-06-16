#' @name perform.seurat.diffexp
#' @aliases perform.seurat.diffexp
#' 
#' @title Perform differential expression
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

perform.seurat.diffexp <- function(object, 
                                   assay = NULL,
                                   test = 'wilcox', 
                                   identity = NULL,
                                   ident.1 = NULL,
                                   ident.2 = NULL,
                                   cells.1 = NULL,
                                   cells.2 = NULL,
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
  
  if(!is.null(cells.1)) {
    
    if(!is.vector(cells.1)) {
      
      cat(crayon::cyan('cells.1 must be a vector \n'))
      return(NULL)
      
    } else if(!cells.1 %in% colnames(object@methods[[assay]][[counts]] )) {
      
      cat(crayon::cyan('cells.1 are not contained within the assay \n'))
      return(NULL)
      
    }
    
  } else if (is.null(cells.1)) {
    
    cat(crayon::cyan('Please provide the cells.1 identities \n'))
    return(NULL)
    
  }
  
  if(!is.null(cells.2)) {
    
    if(!is.vector(cells.2)) {
      
      cat(crayon::cyan('cells.2 must be a vector \n'))
      return(NULL)
      
    } else if(!cells.2 %in% colnames(object@methods[[assay]][[counts]] )) {
      
      cat(crayon::cyan('cells.2 are not contained within the assay \n'))
      return(NULL)
      
    }
    
  } else if (is.null(cells.2)) {
    
    cat(crayon::cyan('Please provide the cells.2 identities \n'))
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
  
  if(!is.null(ident.1) && !is.null(identity)) {
    
    if(!is.vector(ident.1)) {
      
      cat(crayon::cyan('ident.1 must be a vector \n'))
      return(NULL)
      
    } else if (!ident.1 %in% identity) {
      
      cat(crayon::cyan('ident.1 is not contained within identity \n'))
      return(NULL)
      
    } else {
      
      seuobj$clusters <- identity
      Seurat::Idents(seuobj) <- 'clusters'
      
    }
    
  }
  
  if(!is.null(ident.2) && !is.null(identity)) {
    
    if(!is.vector(ident.2)) {
      
      cat(crayon::cyan('ident.2 must be a vector \n'))
      return(NULL)
      
    } else if (!ident.2 %in% identity) {
      
      cat(crayon::cyan('ident.2 is not contained within identity \n'))
      return(NULL)
      
    }
    
  }
  
  if(!is.null(latent.vars)) {
    
    if(!is.character(latent.vars)) {
      
      cat(crayon::cyan('Latent.vars must be character(s)\n'))
      return(NULL)
      
    } else if (!latent.vars %in% names(object@sample_metadata)) {
      
      cat(crayon::cyan('Latent.vars do not exist in object@sample_metadata \n'))
      return(NULL)
      
    }
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    if(!is.null(identity) && !is.null(ident.1) && !is.null(ident.2)) {
      
      results <- Seurat::FindMarkers(object = seuobj, ident.1 = ident.1, ident.2 = ident.2,  test.use = test, latent.vars = latent.vars, ...)
      
    } else if (!is.null(cells.1) && !is.null(cells.2)) {
      
      results <- Seurat::FindMarkers(object = seuobj, cells.1 = cells.1, cells.2 = cells.2,  test.use = test, latent.vars = latent.vars, ...)
      
    }
    
  } else {
    
    if(!is.null(identity) && !is.null(ident.1) && !is.null(ident.2)) {
      
      results <- Seurat::FindMarkers(object = seuobj, ident.1 = ident.1, ident.2 = ident.2,  test.use = test, ...)
      
    } else if (!is.null(cells.1) && !is.null(cells.2)) {
      
      results <- Seurat::FindMarkers(object = seuobj, cells.1 = cells.1, cells.2 = cells.2,  test.use = test, ...)
      
    }
    
  }
  
  return(results)
  
}