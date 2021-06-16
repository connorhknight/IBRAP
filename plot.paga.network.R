
plot.paga.network <- function(result, title='', layout = "layout_nicely") {
  
  tmp <- as.matrix(reticulate::py_to_r(result$connectivities_tree))
  tmp[tmp > 0] = 1
  tmp <- igraph::graph.adjacency(adjmatrix = tmp)
  
  tmp1 <- as.matrix(reticulate::py_to_r(result$connectivities_tree))
  
  net <- igraph::get.edgelist(tmp)
  
  stren <- list()
  
  for(i in 1:nrow(net)) {
    
    stren[[i]] <- tmp1[net[i,1],net[i,2]]
    
  }
  
  net <- cbind(net, unlist(stren))
  
  colnames(net) <- c('from', 'to', 'width')

  node <- data.frame(unique(as.character(result$groups)))
  node <- cbind(node, as.integer(node[,1]))
  colnames(node) <- c('label', 'id')
  node[,'shape'] <- 'circle'
  node[,'font.size'] <- 50
  node[,'color'] <- colorspace::qualitative_hcl(n = length(unique(node[,2])), palette = 'Dark 3')
  
  edge <- as.data.frame(net)
  edge$width <- as.numeric(edge$width)*5
  edge[,'color'] <- 'black'
  edge$from <- edge$from - 1
  edge$to <- edge$to - 1
  
  node$label <- as.character(as.numeric(node$label)-1)
  node$id <- node$id-1
  
  p <- visNetwork::visNetwork(nodes = node, edges = edge, main = title) %>% 
    visNetwork::visIgraphLayout(layout = layout) %>% 
    visNetwork::visEdges(arrows = "middle")
  
  print(p)
  
  }

