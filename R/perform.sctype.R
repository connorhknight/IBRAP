perform.sctype(object, assay='RAW', slot='counts', c, db_=db_, tissue=tissue)) {
  
  gene_sets_prepare <- function(path_to_db_file, cell_type){
    
    cell_markers <- 
    cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
    cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
    
    # correct gene symbols from the given DB (up-genes)
    cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
      
      markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
      markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
      markers_all = sort(markers_all)
      
      if(length(markers_all) > 0){
        markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
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
        markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
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
  
  sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
    
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
  
  sc_type <- function(obj, clusters, db_=db_, tissue=tissue) {
    
    gs_list = gene_sets_prepare(path_to_db_file = , tissue)
    
    es.max = sctype_score(scRNAseqData = obj@methods$SCRAN@norm.scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
    
    cL_resutls = do.call("rbind", lapply(unique(obj@sample_metadata$clusters), function(cl){
      print(cl)
      es.max.cl = sort(rowSums(es.max[,rownames(obj@sample_metadata[obj@sample_metadata$clusters==cl,])]),decreasing=!0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(obj@sample_metadata$clusters==cl)), 10)
    }))
    
    sctype_scores = cL_resutls
    sctype_scores = group_by(sctype_scores, cluster) 
    top_n(sctype_scores, n = 1, wt = scores)  
    
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    
    obj@sample_metadata$customclassif = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]; 
      obj@sample_metadata$customclassif[obj@sample_metadata$clusters == j] = as.character(cl_type$type[1])
    }
    
    return(obj)

}