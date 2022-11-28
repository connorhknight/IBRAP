#' @name perform.celltalker
#' @aliases perform.celltalker
#' 
#' @title Infers cell-to-cell communications.
#'
#' @description Find ligand-receptor communcations between cell types. 
#' 
#' @param object IBRAP S4 class object
#' @param assay A character string containing indicating which assay to use
#' @param slot String. indicating which slot within the assay should be sourced
#' @param clust.method String. can be either metadata or a dataframe stored in cluster_assignments
#' @param column String. name of the column to access in the clust.method
#' @param number_cells_required Number of cells per group required to perform
#' analysis of ligand/receptor interactions. Defaults to 100.
#'
#' @param min_expression Minimum expression in counts to consider a ligand or
#' receptor for interactions analysis. A sensible default is set to 1000, but is
#' dataset dependent. This is meant to filter out lowly expressed ligands and
#' receptors.
#'
#' @param max_expression Maxmium expression in counts to consider a ligand or
#' receptor for interactions analysis. A sensible default is set to 20000, but is
#' dataset dependent. This is meant to filter out ubiquitously expressed ligands
#' and receptors.
#'
#' @param scramble_times Number of times to scamble ligand/receptor interactions to
#' create a background distribution for statistical comparison.
#' 
#' @usage add.cell.cycle(object = obj, assay = 'RAW', slot = 'counts')
#' 
#' @return dataframe of inferred cell-to-cell communications
#'
#' @export

perform.celltalker <- function(object, assay='RAW', slot='counts', clust.method, column,
                               number_cells_required = 100, min_expression = 1000,
                               max_expression = 20000, scramble_times = 10) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('Object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('Assay must be character string\n')
    
  }
  
  if(!assay %in% names(object@methods)) {
    
    stop('Assay not contained in the supplied IBRAP object \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('Slot must be character string\n')
    
  }
  
  if(!slot %in% c('counts','normalised','norm.scaled')) {
    
    stop('Slot must either be: counts, normalised or norm.scaled. However, counts is highly recommended! \n')
    
  }
  
  if(!is.character(clust.method)) {
    
    stop('clust.method must be character string\n')
    
  }
  
  if(!is.character(column)) {
    
    stop('column must be character string\n')
    
  }
  
  if(clust.method == 'metadata') {
    
    if(!column %in% names(object@sample_metadata)) {
      
      stop('specified column not contained in metadata \n')
      
    }
    
  } else if(!clust.method %in% names(object[[assay]]@cluster_assignments)) {
    
    stop('could not find the clust.method, it must either be metadata or one of the cluster assignments dataframes stored in clister_assignments \n')
    
  } else if(clust.method %in% names(object[[assay]]@cluster_assignments)) {
    
    stop('specified column not contained in the supplied clusters_assignment \n')
    
  }
  
  if(clust.method == 'metadata') {
    
    variable <- object[[column]]
    
  } else {
    
    variable <- object[[assay]]@cluster_assignments[[column]]
    
  }
  
  seuobj <- Seurat::CreateSeuratObject(counts = object[[assay]][[slot]])
  
  seuobj@meta.data[[column]] <- variable 
  
  if(slot == 'counts') {
    
    seuobj@assays$RNA@data <- Seurat::LogNormalize(data = seuobj@assays$RNA@counts)
    
  } else {
    
    seuobj@assays$RNA@data <- object[[assay]][[slot]]
  }
  
  
  results <- celltalker::celltalk(input_object=seuobj, metadata_grouping = column, ligand_receptor_pairs = celltalker::ramilowski_pairs, 
                                  number_cells_required = number_cells_required, min_expression = min_expression, max_expression = max_expression, 
                                  scramble_times = scramble_times)
  
  return(results)
  
}
