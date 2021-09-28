#' @name perform.GO.enrichment
#' @aliases perform.GO.enrichment
#' 
#' @title Gene Ontology enrichment 
#'
#' @description Performs gene ontology enrichment for individual cluster differential expression results. 
#' 
#' @param result A database containing the differential expression results
#' @param whichOnto Character. specifying one of the three GO ontologies, namely: "BP", "MF", "CC". Default = 'BP'
#' @param feasibleGenes Character vector. vector containing a subset of gene identifiers. Only these genes will be used to annotate GO terms. Default value is NULL which means that there are no genes filtered.
#' @param mapping Character. The name of the Bioconductor package containing the gene mappings for a specific organism. For example: mapping = "org.Hs.eg.db".
#' @param ID Character. Specify the gene identifier to use. Currently only the following identifiers can be used: c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene")
#' @param nodeSize Numerical. Minimum number of genes required to consider a GO term. Default = 5
#' @param statistic Character. Which statistic to use when testing for significant GO terms, options: 'fisher', 'ks', 't', 'globaltest', 'sum', 'ks.ties'. Default = 'ks'
#' @param algorithm Character. Which algorithm to use when testing for significant GO terms, options: 'classic', 'elim', 'weight', 'weight01', 'lea', 'parentchild'. Default = 'classic'
#' @param rank.cutoff Numerical. Which cut off to apply for pathway significance, this value will change according to the statistic applied. Default = 0.001
#' @param gene.col Character. Which column name within differential expression results contains the genes. Default = 'gene.col'
#' @param pval.col Character. Which column name within differential expression results contains the p values. Default = 'p_val'
#' @param cluster.col Character. Which column name within differential expression results contains the cluster assignments. Default = 'cluster'
#' @param n.top.pathways Numerical. How many top pathways per group should be retained. Default = 10
#' 
#' @return A dataframe containing the top enriched pathways for each cluster
#' 
#' @examples 
#' 
#' SCT_DE <- perform.seurat.diffexp.all(object = object, assay = 'SCT', test = 'MAST', identity = object@sample_metadata$celltype, latent.vars = 'original.project')
#' 
#' SCT_DE_GO <- perform.GO.enrichment(result = SCT_DE)
#' 
#' plot.GO.output(result = SCT_DE_GO) + ggplot2::ggtitle(label = 'SCT')
#'
#' @export
#' 

perform.GO.enrichment <- function(result, 
                                  whichOnto = 'BP',
                                  feasibleGenes = NULL,
                                  mapping = 'org.Hs.eg.db',
                                  ID = 'symbol', 
                                  nodeSize = 5, 
                                  algorithm = 'classic',
                                  statistic = 'ks',
                                  rank.cutoff = 0.001, 
                                  gene.col = 'gene', 
                                  pval.col = 'p_val', 
                                  cluster.col = 'cluster',
                                  n.top.pathways = 10) {
  
  if(!is.data.frame(result)) {
    
    stop('results must be in data.frame format \n')
    
  } 
  
  if(!is.character(gene.col)) {
    
    stop('Gene.col must be a character string \n')
    
  } else if (!gene.col %in% colnames(result)) {
    
    stop('gene.col is not in colnames(results) \n')
    
  }
  
  if(!is.character(pval.col)) {
    
    stop('pval.col must be a character string \n')
    
  } else if (!pval.col %in% colnames(result)) {
    
    stop('pval.col is not in colnames(results) \n')
    
  }
  
  if(!is.character(cluster.col)) {
    
    stop('cluster.col must be a character string \n')
    
  } else if (!cluster.col %in% colnames(result)) {
    
    stop('cluster.col is not in colnames(results) \n')
    
  }
  
  if(!is.numeric(nodeSize)) {
    
    stop('nodeSize must be numerical\n')
    
  }
  
  if(!is.numeric(rank.cutoff)) {
    
    stop('rank.cutoff must be numerical\n')
    
  }
  
  if(!is.numeric(n.top.pathways)) {
    
    stop('n.top.pathways must be numerical\n')
    
  }
  
  require(topGO)
  
  selection <- function(x) TRUE
  allGO2genes <- topGO::annFUN.org(whichOnto=whichOnto, feasibleGenes=feasibleGenes, mapping=mapping, ID=ID)
  
  GO_outputs <- list()
  
  for(x in unique(result[,cluster.col])) {
    
    cluster_x <- result[result[,cluster.col] == x,]
    
    genes <- setNames(cluster_x[,pval.col], cluster_x[,gene.col])

    GOdata <- new("topGOdata", ontology="BP", allGenes=genes,
                  annot=annFUN.GO2genes, GO2genes=allGO2genes,
                  geneSel=selection, nodeSize=nodeSize)

    results.ks <- runTest(GOdata, algorithm=algorithm, statistic=statistic)

    goEnrichment <- GenTable(GOdata, rank=results.ks, topNodes=n.top.pathways)
    
    temp <- goEnrichment$rank
    
    for(t in 1:length(temp)) {
      
      if(stringr::str_detect(string = temp[t], pattern = '< ')) {
        
        temp[t] <- strsplit(x = temp[t], split = ' ')[[1]][2]

      } else if(stringr::str_detect(string = temp[t], pattern = '<')) {
        
        temp[t] <- strsplit(x = temp[t], split = '<')[[1]][2]

      }
      
    }
    
    goEnrichment$rank <- as.numeric(temp)

    goEnrichment <- goEnrichment[goEnrichment$rank<rank.cutoff,]

    goEnrichment <- goEnrichment[,c("GO.ID","Term","rank")]
    
    if(nrow(goEnrichment) > 0) {
      
      goEnrichment[,'cluster'] <- x
      
      if(is.numeric(x)) {
        
        x <- x + 1
        
      }
      
      GO_outputs[[x]] <- goEnrichment
      
    }
    
  }
  
  cat_df <- GO_outputs[[1]]
  
  for(x in 2:length(GO_outputs)) {
    
    if(!is.null(GO_outputs[[x]])) {
      
      cat_df <- merge(cat_df, GO_outputs[[x]], all = T)
      
    }

  }
  
  cat(crayon::cyan(paste0(Sys.time(), ': completed gene ontology enrichment \n')))
  
  return(cat_df)
  
}
