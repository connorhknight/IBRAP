#' @title An S4 class object of method-assays
#'
#' @description An S4 class object that contains data matrices, HVGs, feature metadata, neighbourhood graphs, reduction embeddings, cluster assignments, and benchmarking results.
#' 
#' @slot counts Raw counts matrix in dgCMatrix format.
#' @slot normalised Normalised matrix in dgCMatrix format. 
#' @slot norm.scaled HVG subset of normalised and scaled matrix in matrix format.
#' @slot highly.variable.genes list of highly variable genes in character format. 
#' @slot feature_metadata data frame containing feature-level metadata. 
#' @slot graphs list of generated neighbourhood graphs (distances)
#' @slot computational_reductions list of reductions to improve computational efficiency, i.e. PCA or dbMAP
#' @slot integration_reductions list of integration output in reduction format
#' @slot visualisation_reductions list of visualisation reductions, i.e. t-SNE, UMAP, etc. 
#' @slot cluster_assignments list of cluster assignment output. A dataframe is assigned to each method containing subsequent parameter changes. 
#' @slot benchmarking_results list of benchmarking results corresponding to the cluster assignment dataframes. 
#' @slot alt_objects objects used in the analysis derived from alternative packages, i.e. SingleCellExperiment, Seurat, Anndata, etc.

setClass(Class = 'methods',
         representation = representation(
           counts = 'dgCMatrix', 
           normalised = 'dgCMatrix', 
           norm.scaled = 'matrix',
           highly.variable.genes = 'character',
           feature_metadata = 'data.frame',
           graphs = 'list',
           computational_reductions = 'list',
           integration_reductions = 'list',
           visualisation_reductions = 'list',
           cluster_assignments = 'list',
           benchmark_results = 'list',
           alt_objects = 'list'
         ))