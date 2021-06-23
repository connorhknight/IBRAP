sc <- reticulate::import(module = 'scanpy', convert = T)

scobj <- sc$AnnData(X = t(as.matrix(traj@methods$SCT@normalised)))
scobj$obs_names <- as.factor(colnames(traj@methods$SCT@normalised))
scobj$var_names <- as.factor(rownames(traj@methods$SCT@normalised))

if(length(colnames(as.data.frame(traj@sample_metadata))) >= 1) {
  
  pd <- reticulate::import('pandas')
  
  scobj$obs <- pd$DataFrame(data = as.data.frame(traj@sample_metadata))
  
}

# scobj$obs[['clusters']] <- traj@sample_metadata$celltype

scobj$obsm$update(X_pca = as.matrix(traj@methods$SCT@visualisation_reductions$pca))

sc$pp$neighbors(adata = scobj, n_neighbors=as.integer(5), n_pcs=as.integer(10))

sc$tl$paga(adata = scobj, groups = 'celltype')

sc$pl$paga(adata = scobj, color='celltype',edge_width_scale=0.2, threshold=0)

sc$pl$paga(adata = scobj, color='celltype',threshold=0,
           solid_edges='connectivities_tree', root = 5,
           dashed_edges='connectivities',layout='rt',
           node_size_scale=0.5,node_size_power=0.9,
           max_edge_width=0.7,fontsize=8)

sc$tl$umap(adata = scobj, min_dist = 0.3, maxiter = as.integer(1000), init_pos = 'paga')

sc$pl$umap(adata = scobj, color='celltype', legend_loc='on data')








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







