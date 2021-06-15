#' @name perform.GO.enrichment
#' @aliases perform.GO.enrichment
#' 
#' @title Plots slingshot results
#'
#' @description Plots the results of the slingshot analysis using ggplots2
#' 
#' @return A ggplot of the reduced cellular embbedings and trajectories. 
#'
#' @export
#' 

perform.GO.enrichment <- function(result, 
                                  pval.cutoff = 0.001, 
                                  gene.col = 'gene', 
                                  pval.col = 'p_val', 
                                  cluster.col = 'cluster', 
                                  nodeSize = 5, 
                                  n.top.pathways = 10) {
  
  selection <- function(x) TRUE
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Hs.eg.db", ID="symbol")
  
  GO_outputs <- list()
  
  for(x in unique(result[,cluster.col])) {
    
    cluster_x <- result[result[,cluster.col] == x,]
    
    genes <- setNames(cluster_x[,pval.col], cluster_x[,gene.col])
    
    GOdata <- new("topGOdata", ontology="BP", allGenes=genes,
                  annot=annFUN.GO2genes, GO2genes=allGO2genes,
                  geneSel=selection, nodeSize=nodeSize)
    
    results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
    goEnrichment <- GenTable(GOdata, rank=results.ks, topNodes=n.top.pathways)
    
    temp <- goEnrichment$rank
    
    for(t in 1:length(temp)) {
      
      if(stringr::str_detect(string = temp[t], pattern = '< ')) {
        
        temp[t] <- strsplit(x = temp[t], split = ' ')[[1]][2]
        
      }
      
    }
    
    goEnrichment$rank <- as.numeric(temp)
    goEnrichment <- goEnrichment[goEnrichment$rank<0.05,]
    goEnrichment <- goEnrichment[,c("GO.ID","Term","rank")]
    
    goEnrichment[,'cluster'] <- x
    
    if(is.numeric(x)) {
      
      x <- x + 1
      
    }
    
    GO_outputs[[x]] <- goEnrichment
    
  }
  
  cat_df <- GO_outputs[[1]]
  
  for(x in 2:length(GO_outputs)) {
    
    cat_df <- merge(cat_df, GO_outputs[[x]], all = T)
    
  }
  
  return(cat_df)
  
}
