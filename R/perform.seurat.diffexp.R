#' @name perform.diffexp
#' @aliases perform.diffexp
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

perform.diffexp <- function(object, 
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
    
    stop('Object must be IBRAP class \n')
    
  }
  
  if(!is.null(assay)) {
    
    if(!is.character(assay)) {
      
      stop('Assay must be character string \n')
      
    }
    
  } else if (is.null(assay)) {
    
    stop('Please indicate which assay to access \n')
    
  }
  
  if(!is.character(test)) {
    
    stop('Test must be character string \n')
    
  }
  
  if(!is.null(cells.1)) {
    
    if(!is.vector(cells.1)) {
      
      stop('cells.1 must be a vector \n')
      
    } else if(!cells.1 %in% colnames(object@methods[[assay]][[counts]] )) {
      
      stop('cells.1 are not contained within the assay \n')
      
    }
    
  } else if (is.null(cells.1)) {
    
    stop('Please provide the cells.1 identities \n')
    
  }
  
  if(!is.null(cells.2)) {
    
    if(!is.vector(cells.2)) {
      
      stop('cells.2 must be a vector \n')
      
    } else if(!cells.2 %in% colnames(object@methods[[assay]][[counts]] )) {
      
      stop('cells.2 are not contained within the assay \n')
      
    }
    
  } else if (is.null(cells.2)) {
    
    stop('Please provide the cells.2 identities \n')
    
  }
  
  seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
  seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
  seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
  
  if(!is.null(identity)) {
    
    if(!is.vector(identity)) {
      
      stop('Identity must be a vector \n')
      
    } else if(length(identity) != ncol(object@methods[[assay]]@counts)) {
      
      stop('Identity length does not match the number of cells \n')
      
    }
    
  } else if (is.null(identity)) {
    
    stop('Please provide the cell identities \n')
    
  }
  
  if(!is.null(ident.1) && !is.null(identity)) {
    
    if(!is.vector(ident.1)) {
      
      stop('ident.1 must be a vector \n')
      
    } else if (!ident.1 %in% identity) {
      
      stop('ident.1 is not contained within identity \n')
      
    } else {
      
      seuobj$clusters <- identity
      Seurat::Idents(seuobj) <- 'clusters'
      
    }
    
  }
  
  if(!is.null(ident.2) && !is.null(identity)) {
    
    if(!is.vector(ident.2)) {
      
      stop('ident.2 must be a vector \n')
      
    } else if (!ident.2 %in% identity) {
      
      stop('ident.2 is not contained within identity \n')
      
    }
    
  }
  
  if(!is.null(latent.vars)) {
    
    if(!is.character(latent.vars)) {
      
      stop('Latent.vars must be character(s)\n')
      
    } else if (!latent.vars %in% names(object@sample_metadata)) {
      
      stop('Latent.vars do not exist in object@sample_metadata \n')
      
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
