#' @name perform.paga.trajectory
#' @aliases perform.paga.trajectory
#' 
#' @title Create an IBRAP class object 
#'
#' @description Creates and produces project metadata into an IBRAP S4 class object
#' 
#' @usage createIBRAPobject(counts = counts, original.project = 'project_1', method.name = 'RAW', meta.data = df, min.cells = 3, min.features = 200)
#' 
#' @return IBRAP S4 class object containing raw counts and metadata
#'
#' @export

perform.paga.trajectory <- function(object, 
                                    assay, 
                                    reduction, 
                                    clust.method = NULL, 
                                    column, 
                                    n_neighbors = 10,
                                    n_components = NULL,
                                    paga.model = 'v1.2',
                                    ...) {
  
  sc <- reticulate::import('scanpy', convert = T)
  
  scobj <- sc$AnnData(X = t(as.matrix(object@methods[[assay]]@counts)))
  scobj$obs_names <- as.factor(colnames(object))
  scobj$var_names <- as.factor(rownames(object))
  
  if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
    
    pd <- reticulate::import('pandas')
    
    scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))
    
  }
  
  if(!is.null(clust.method)) {
    
    scobj$obs[['clusters']] <- object@methods[[assay]]@cluster_assignments[[clust.method]][[column]]
    
  }
  
  reduction.list <- list()
  red.names <- c(names(object@methods[[assay]]@computational_reductions), 
                 names(object@methods[[assay]]@integration_reductions),
                 names(object@methods[[assay]]@visualisation_reductions))
  
  for(i in red.names) {
    
    if(i %in% names(object@methods[[assay]]@computational_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@computational_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@integration_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@integration_reductions[[i]]
      
    }
    
    if(i %in% names(object@methods[[assay]]@visualisation_reductions)) {
      
      reduction.list[[i]] <- object@methods[[assay]]@visualisation_reductions[[i]]
      
    }
  }
  
  scobj$obsm$update(X_pca = as.matrix(reduction.list[[reduction]]))
  
  if(is.null(n_components)) {
    
    n_components <- ncol(reduction.list[[reduction]])

  }
  
  sc$pp$neighbors(adata = scobj, n_neighbors=as.integer(n_neighbors), n_pcs=as.integer(n_components))
  
  if(is.null(clust.method)) {
    
    sc$tl$paga(adata = scobj, groups = column, model = paga.model)
    
  } else if (!is.null(clust.method)) {
    
    sc$tl$paga(adata = scobj, groups = 'clusters', model = paga.model)
    
  }
  
  tmp <- scobj$uns[['paga']]
  
  if(is.null(clust.method)) {

    tmp$groups <- scobj$obs[[column]]

  } else if (!is.null(clust.method)) {
    
    tmp$groups <- scobj$obs[['clusters']]

  }
  
  return(tmp)
  
}

paga.result <- perform.paga.trajectory(object = marrow, assay = 'SCANPY', reduction = 'scanorama_umap', clust.method = 'scanorama_Louvain', column = 'RNA_snn_res.1')


