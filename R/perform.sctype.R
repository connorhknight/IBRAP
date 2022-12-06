#' @name perform.sctype
#' @aliases perform.sctype
#' 
#' @title Performs sctype 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String defining which slot in the assay to supply to Scanorama. Default = NULL
#' @param db Character. String indicating where to collect database from. So far we only use the standard scType reference.
#' @param tissue Character.String which tissue to subset.
#' 
#' @return cell type annotation in relation to the clusterign categories that were provided
#'
#' @export

perform.sctype <- function(object, assay='RAW', slot='counts', clust.method, column, scaled=T,
                           db="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", 
                           tissue) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('reduction: ', x, 'does not exist\n'))
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be a character string\n')
    
  }
  
  if(!slot %in% c('counts', 'normalised', 'norm.scaled')) {
    
    stop('slot does not exist\n')
    
  }
  
  if(!is.character(db)) {
    
    stop('db must be character string \n')
    
  }
  
  if(!is.character(tissue)) {
    
    stop('tissue must be character string \n')
    
  }
  
  gene_sets_prepare <- function(path_to_db_file, cell_type){
    
    cell_markers = openxlsx::read.xlsx(path_to_db_file)
    cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    # correct gene symbols from the given DB (up-genes)
    cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        suppressMessages({markers_all = unique(na.omit(HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol))})
        paste0(markers_all, collapse=",")
      } else {
        ""
      }
    })
    
    # correct gene symbols from the given DB (down-genes)
    cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        suppressMessages({markers_all = unique(na.omit(HGNChelper::checkGeneSymbols(markers_all)$Suggested.Symbol))})
        paste0(markers_all, collapse=",")
      } else {
        ""
      }
    })
    
    cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
    cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
    gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
    
    list(gs_positive = gs, gs_negative = gs2)
  }
  
  sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0){
    
    # check input matrix
    if(!is.matrix(scRNAseqData)){
      warning("scRNAseqData doesn't seem to be a matrix")
    } else {
      if(sum(dim(scRNAseqData))==0){
        warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
      }
    }
    
    # marker sensitivity
    marker_stat = sort(table(unlist(gs)), decreasing = T); 
    marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                    gene_ = names(marker_stat), stringsAsFactors = !1)
    
    # convert gene names to Uppercase
    if(gene_names_to_uppercase){
      rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
    }
    
    # subselect genes only found in data
    names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
    gs = lapply(1:length(gs), function(d_){ 
      GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
    gs2 = lapply(1:length(gs2), function(d_){ 
      GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
    names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
    cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
    
    # z-scale if not
    if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
    
    # multiple by marker sensitivity
    for(jj in 1:nrow(cell_markers_genes_score)){
      Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
    }
    
    # subselect only with marker genes
    Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
    
    # combine scores
    es = do.call("rbind", lapply(names(gs), function(gss_){ 
      sapply(1:ncol(Z), function(j) {
        gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
        sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
        if(is.na(sum_t2)){
          sum_t2 = 0;
        }
        sum_t1 + sum_t2
      })
    })) 
    
    dimnames(es) = list(names(gs), colnames(Z))
    es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
    
    es.max
  }
  
  if(clust.method != 'metadata') {
    
    clusters <- object[[assay]][[clust.method]][[column]]
    
  } else {
    
    clusters <- object[[column]]
    
  }

  gs_list = gene_sets_prepare(path_to_db_file = db, cell_type = tissue)

  if(isTRUE(scaled)) {

    es.max = sctype_score(scRNAseqData = object[[assay]][[slot]], scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

  } else if(isFALSE(scaled)) {
    
    es.max = sctype_score(scRNAseqData = object[[assay]][[slot]], scaled = F, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
    
  }

  cL_resutls = do.call("rbind", lapply(unique(clusters), function(cl){
    
    es.max.cl = sort(rowSums(es.max[,rownames(object@sample_metadata[clusters==cl,])]),decreasing=!0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(clusters==cl)), 10)
    
  }))
  
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  
  object@sample_metadata[,paste0('scType_', assay, '_', slot)] = ""
  
  for(j in unique(sctype_scores$cluster)){
    
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    object@sample_metadata[,paste0('scType_', assay, '_', slot)][clusters == j] = as.character(cl_type$type[1])
    
  }
  
  return(object)
  
}
