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
  
  if(!is.null(identity)) {
    
    if(!is.vector(identity)) {
      
      stop('Identity must be a vector \n')
      
    } else if(length(identity) != ncol(object@methods[[assay]]@counts)) {
      
      stop('Identity length does not match the number of cells \n')
      
    }
    
  } else if (is.null(identity)) {
    
    stop('Please provide the cell identities \n')
    
  }

  
  if(!is.null(latent.vars)) {
    
    if(!is.character(latent.vars)) {
      
      stop('Latent.vars must be character(s)\n')
      
    } else if(!latent.vars %in% names(object@sample_metadata)) {
      
      stop('Latent.vars do not exist in object@sample_metadata \n')
      
    }
    
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
    seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
    seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
    
    seuobj$cluster <- identity
    Seurat::Idents(seuobj) <- 'cluster'
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met

    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, latent.vars = latent.vars, ...)
    
  } else {
    
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
    seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
    seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
    
    seuobj$cluster <- identity
    Seurat::Idents(seuobj) <- 'cluster'
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, ...)
    
  }
  
  return(results)
  
}