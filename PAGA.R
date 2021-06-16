sc <- reticulate::import(module = 'scanpy', convert = T)

scobj <- sc$AnnData(X = t(as.matrix(marrow_subset@methods$RAW@counts)))
scobj$obs_names <- as.factor(colnames(marrow_subset))
scobj$var_names <- as.factor(rownames(marrow_subset))

if(length(colnames(as.data.frame(marrow_subset@sample_metadata))) >= 1) {
  
  pd <- reticulate::import('pandas')
  
  scobj$obs <- pd$DataFrame(data = as.data.frame(marrow_subset@sample_metadata))
  
}

scobj$obs[['clusters']] <- marrow_subset@methods$SCT@cluster_assignments$scanorama_Louvain$RNA_snn_res.1

scobj$obsm$update(X_pca = as.matrix(marrow_subset@methods$SCT@integration_reductions$scanorama))

sc$pp$neighbors(adata = scobj, n_neighbors=as.integer(50), n_pcs=as.integer(50))

sc$tl$paga(adata = scobj, groups = 'clusters')














g <- scobj$uns[['paga']]

tmp <- as.matrix(reticulate::py_to_r(g$connectivities_tree))

tmp[tmp > 0] = 1

library(igraph)

m <- graph.adjacency(adjmatrix = tmp)

tmp1 <- as.matrix(reticulate::py_to_r(g$connectivities_tree))

net <- get.edgelist(m)

stren <- list()

for(i in 1:nrow(net)) {
  
  stren[[i]] <- tmp1[net[i,1],net[i,2]]
  
}

net <- cbind(net, unlist(stren))

colnames(net) <- c('from', 'to', 'width')

node <- data.frame(unique(as.integer(marrow_subset@methods$SCANPY@cluster_assignments$scanorama_Louvain$RNA_snn_res.1)))
node <- cbind(node, as.character(node[,1]))
colnames(node) <- c('id', 'label')
node[,'shape'] <- 'circle'
node[,'font.size'] <- 50
node[,'color'] <- colorspace::qualitative_hcl(n = length(unique(node[,2])), palette = 'Dark 3')

edge <- as.data.frame(net)
edge$width <- edge$width*5
edge[,'color'] <- 'black'



colorspace::qualitative_hcl(n = length(unique(node[,2])), palette = 'Dark 3')







