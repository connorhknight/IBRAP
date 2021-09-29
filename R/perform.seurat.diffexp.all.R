#' @name perform.diffexp.all
#' @aliases perform.diffexp.all
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
#' @examples
#' 
#' SCT_DE <- perform.seurat.diffexp.all(object = object, 
#'                                      assay = 'SCT', 
#'                                      test = 'MAST', 
#'                                      identity = object@sample_metadata$celltype, 
#'                                      latent.vars = 'original.project')
#'                                      
#' SCRAN_DE <- perform.seurat.diffexp.all(object = object, 
#'                                        assay = 'SCRAN', 
#'                                        test = 'MAST', 
#'                                        identity = object@sample_metadata$celltype, 
#'                                        latent.vars = 'original.project')
#'                                        
#' SCANPY_DE <- perform.seurat.diffexp.all(object = object, 
#'                                         assay = 'SCRAN', 
#'                                         test = 'MAST', 
#'                                         identity = object@sample_metadata$celltype, 
#'                                         latent.vars = 'original.project')
#'
#' @export
#' 

perform.diffexp.all <- function(object, 
                                assay = NULL,
                                clust.method = NULL,
                                column = NULL,
                                test = 'wilcox', 
                                latent.vars = NULL,
                                ...) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('Object must be IBRAP class \n')
    
  }
  
  if(!is.null(assay)) {
    
    if(!is.character(assay)) {
      
      stop('Assay must be character string \n')
      
    } else if (is.character(assay)) {
      
      if(!assay %in% names(object@methods)) {
        
        stop('assay not founnd in object@methods \n')
        
      }
      
    }
    
  } else if (is.null(assay)) {
    
    stop('Please indicate which assay to access \n')
    
  }
  
  if(!is.character(test)) {
    
    stop('Test must be character string \n')
    
  }
  
  if(!is.null(clust.method)) {
    
    if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
      
      stop('clust.method is not present in cluster_assignments \n')
      
    }
    
  } else if (is.null(clust.method)) {
    
    stop('please provide a clust.method from cluster_assignments \n')
    
  }
  
  if(!is.null(column)) {
    
    if(!column %in% names(object@methods[[assay]]@cluster_assignments[[clust.method]])) {
      
      stop('column is not present in clust.method \n')
      
    }
    
  } else if (is.null(column)) {
    
    stop('please provide a colum in the dataframe of your clust.method \n')
    
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
    
    if(clust.method != 'metadata') {
      
      seuobj$cluster <- object@methods[[assay]]@cluster_assignments[[clust.method]][,column]
      
    } else if(clust.method == 'metadata') {
      
      seuobj$cluster <- object@sample_metadata[,column]
      
    }
    
    Seurat::Idents(seuobj) <- 'cluster'
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, latent.vars = latent.vars, ...)
    
  } else {
    
    seuobj <- Seurat::CreateSeuratObject(counts = object@methods[[assay]]@counts)
    seuobj@assays$RNA@data <- object@methods[[assay]]@normalised
    seuobj@assays$RNA@scale.data <- object@methods[[assay]]@norm.scaled
    
    if(clust.method != 'metadata') {
      
      seuobj$cluster <- object@methods[[assay]]@cluster_assignments[[clust.method]][,column]
      
    } else if(clust.method == 'metadata') {

      seuobj$cluster <- object@sample_metadata[,column]

    }
    
    Seurat::Idents(seuobj) <- 'cluster'
    
    met <- merge(seuobj@meta.data, object@sample_metadata, by = 0)
    rownames(met) <- colnames(seuobj)
    seuobj@meta.data <- met
    
    results <- Seurat::FindAllMarkers(object = seuobj, test.use = test, ...)
    
  }
  
  return(results)
  
}



