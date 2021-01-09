require(Matrix)
library(scran)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(ggplot2)
library(umap)
library(gmodels)
library(Seurat)
library(RColorBrewer)
library(plotly)
library(MASS)
library(Rtsne)
library(PCAtools)
library(uwot)
library(tibble)
library(cluster)
library(clValid)
library(SC3)

Read10X_output <- function(directory) {
  count.file <- paste0(directory, '/matrix.mtx')
  genes.file <- paste0(directory, '/genes.tsv')
  barcodes.files <- paste0(directory, '/barcodes.tsv')
  genes_ensembl <- read.table(file = genes.file, sep = '\t', header = FALSE)
  barcodes <- read.table(file = barcodes.files, sep = '\t', header = FALSE)
  mm <- readMM(count.file)
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  mm <- as.matrix(mm)
  return(mm)
}

create_sce_experiment <- function(counts, project.name) {
  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$original.project <- project.name
  colData(sce)$total.counts <- colSums(assay(sce, 'counts'))
  colData(sce)$total.features <- colSums(assay(sce, 'counts') != 0)
  sce <- sce[!duplicated(rownames(assay(sce, 'counts'))), ]
  return(sce)
}

find_percentage_genes <- function(object, pattern='^MT', column.name) {
  temp <- data.frame(tmp=colSums(assay(object, 'counts')[grep(as.character(pattern),rownames(assay(object, 'counts'))),]) / colData(object)$total.counts * 100)
  colnames(temp) <- c(column.name)
  colData(object) <- cbind(colData(object), temp)
  return(object)
}

scQC_density <- function(object, column, title, cutoff){
  metadata <- as.data.frame(colData(object))
  tmp <- as.data.frame(metadata[ , grepl(column, names(metadata)) ])
  rownames(tmp) <- colnames(object)
  colnames(tmp) <- 'selected_feature'
  tmp$original.project <- metadata$original.project
  plot <- ggplot(tmp, aes(color = original.project, x = selected_feature, fill = original.project)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  ylab("Cell density") +
  xlab(column) +
  ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = (15))) +
  geom_vline(xintercept = cutoff)
  print(plot)
}

scQC_correlate <- function(object, x, y, z){
  metadata <- as.data.frame(colData(object))
  tmp <- data.frame(x=as.numeric(metadata[ , grepl(x, names(metadata)) ]))
  tmp$y <- as.numeric(metadata[ , grepl(y, names(metadata)) ])
  tmp$z <- as.numeric(metadata[ , grepl(z, names(metadata)) ])
  tmp$original.project <- object$original.project
  rownames(tmp) <- colnames(object)
  plot <- ggplot(tmp, aes(x=x, y=y, color=z)) +
    geom_point() +
    scale_colour_gradient(low = "blue", high = "red") +
    stat_smooth(method=lm, formula = y ~ x + z) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 3000) +
    geom_hline(yintercept = 1000) +
    facet_wrap(~original.project) +
    labs(x=x, y=x, color=z)
  print(plot)
}

scQC_scatter <- function(object, x, y, title, xlab, ylab) {
  tmp <- colData(object)
  tmp <- data.frame('one'=tmp[,grepl(x,names(tmp))],'two'=tmp[,grepl(y,names(tmp))])
  plot <- ggplot(tmp, aes(x = one, y = two)) + geom_point(color='red') + ggtitle(title) + xlab(xlab) + ylab(ylab) + theme(plot.title = element_text(hjust = 0.5))
  print(plot)
}

filter_features_cells <- function(object, max.features, min.features, max.mt, min.cell) {
  rowData(object)$n_features <- rowSums(assay(object, 'counts') != 0)
  object <- object[rownames(subset(rowData(object), n_features > as.numeric(min.cell))),]
  object <- object[,rownames(subset(colData(object), total.counts > as.numeric(min.features)))]
  object <- object[,rownames(subset(colData(object), total.features < as.numeric(max.features)))]
  object <- object[,rownames(subset(colData(object), percent.mt < as.numeric(max.mt)))]
  return(object)
}

filter_other_features <- function(object, features=list()) {
  object <- object[!rownames(object) %in% features, ]
  return(object)
}

#####################################################################

normalisation_scran <- function(object) {
  clusters <- quickCluster(object)
  sce.scran <- computeSumFactors(object, clusters=clusters)
  log <- logNormCounts(sce.scran)
  assay(object, 'data') <- assay(log, 'logcounts')
  return(object)
}

normalisation_tpm <- function(object, ...) {

}

normalisation_cpm <- function(object, ...) {
  assay(object, 'data') <- cpm(assay(object, 'counts'), log = FALSE, ...)
  return(object)
}

normalisation_SCT <- function(object, ...) {
  temp <- object
  tmp <- as.Seurat(x = object, counts = 'counts', data = NULL)
  tmp <- SCTransform(object = tmp, do.scale = FALSE, do.center = FALSE,  ...)
  common_rows <- intersect(rownames(tmp), rownames(temp))
  temp <- temp[common_rows,]
  tmp <- tmp[common_rows,]
  assay(temp, 'data') <- as.matrix(tmp@assays$SCT@data)
  return(temp)
}

normalisation_SeuratLog <- function(object, ...) {
  tmp <- as.Seurat(x = object, counts = 'counts', data = NULL)
  results <- NormalizeData(object = tmp, ...)
  assay(object, 'data') <- as.matrix(results@assays$RNA@data)
  return(object)
}

find_HVGs <- function(object, nfeatures, ...) {
  if(typeof(object) != 'S4') {
    return('Must be an S4 SCE object')
  } else {
    tmp <- as.Seurat(object, counts = 'counts', data = 'data')
    f <- FindVariableFeatures(object = tmp, nfeatures = nfeatures, ...)
    object <- object[f@assays$RNA@var.features,]
    return(object)
  }
}

Scale_data <- function(object, unwanted.variance, scale, center, ...) {
  tmp <- as.Seurat(object, counts = 'counts', data = 'data')
  tmp <- ScaleData(tmp, vars.to.regress = unwanted.variance, do.scale = scale, do.center = scale, ...)
  assay(object, 'scaled') <- tmp@assays$RNA@scale.data
  return(object)
}

performPCA <- function(object, ...) {
  t <- assay(object, 'scaled')
  a <- pca(mat = assay(object, 'scaled'), ...)
  reducedDim(object, 'pca') <- as.matrix(a$rotated[,1:50])
  return(object)
}

benchmark_norm <- function(object, DR.assay, distance.method, n.k=2:25, method.name) {
  kmeans_clusters <- data.frame(barcodes=colnames(object))
  DR <- reducedDim(object, as.character(DR.assay))[,1:3]
  dist.matrix <- dist(x = DR, method = distance.method)
  for(v in n.k){
    print(paste0('Clustering, k=', v))
    tmp <- as.data.frame(pam(x = dist.matrix, k = v, diss = TRUE, metric = distance.method, cluster.only = TRUE, do.swap = FALSE))
    tmp <- cbind(rownames(tmp), tmp)
    kmeans_clusters <- cbind(kmeans_clusters, tmp[,2])
    new.names <- c(colnames(kmeans_clusters)[1:length(colnames(kmeans_clusters))-1], paste0(v))
    colnames(kmeans_clusters) <- new.names
  }
  rownames(kmeans_clusters) <- kmeans_clusters$barcodes
  kmeans_clusters <- kmeans_clusters[2:length(colnames(kmeans_clusters))]

  sil.results <- data.frame(average_silhoeutte=NA)
  for (v in colnames(kmeans_clusters)[1:length(colnames(kmeans_clusters))]) {
    print(paste0('Calculating silhouette for k=', v))
    tmp <- silhouette(x = as.numeric(x = as.factor(x = kmeans_clusters[,v])), dist = dist.matrix)
    average <- sum(tmp[,3])/length(tmp[,3])
    sil.results[v,] <- average
  }
  sil.results <- sil.results[complete.cases(sil.results),]
  max.AS <- max(sil.results)

  dunn.results <- data.frame(dunn.index=NA)
  for (p in colnames(kmeans_clusters)[1:length(colnames(kmeans_clusters))]){
    print(paste0('Calculating dunn index for k=', p))
    dunn.results[p,] <- dunn(distance = dist.matrix, clusters = kmeans_clusters[,p])
  }
  dunn.results <- dunn.results[complete.cases(dunn.results),]
  max.dunn <- max(dunn.results)

  print(paste0('Maximum average silhouette width: ', max.AS))
  print(paste0('Maximum dunn index: ', max.dunn))
  df <- data.frame(metrics=c(max.AS, max.dunn))
  rownames(df) <- c('max.AS', 'max.dunn')
  results <- list('max_values'=df, 'pam_clusters'=kmeans_clusters, 'silhouette_results'=sil.results, 'dunn_results'=dunn.results)
  metadata(object)$benchmarking$normalisation <- results
  return(object)
}

scNormalise <- function(object, method = 'scran', scale, center, unwanted.variance, n.features=2000, benchmark, distance.method, ...) {
  if(method == 'scran'){
    object <- normalisation_scran(object, ...)
    object <- find_HVGs(object = object, nfeatures = n.features)
    object <- Scale_data(object = object, unwanted.variance = unwanted.variance, scale = scale, center = center)
    print('Norm method: Scran, completed.')
  }
  if(method == 'logtransform'){
    object <- normalisation_SeuratLog(object, ...)
    object <- find_HVGs(object = object, nfeatures = n.features)
    object <- Scale_data(object = object, unwanted.variance = unwanted.variance, scale = scale, center = center)
    print('Norm method: Seurat LogNorm, completed.')
  }
  if(method == 'SCTransform'){
    object <- normalisation_SCT(object, ...)
    object <- find_HVGs(object = object, nfeatures = n.features)
    object <- Scale_data(object = object, unwanted.variance = unwanted.variance, scale = scale, center = center)
    print('Norm method: SCTransform, completed.')
  }
  if(method == 'cpm'){
    object <- normalisation_cpm(object, ...)
    object <- find_HVGs(object = object, nfeatures = n.features)
    object <- Scale_data(object = object, unwanted.variance = unwanted.variance, scale = scale, center = center)
    print('Norm method: CPM, completed.')
  }
  temp <- object
  if(method == 'all'){
    scran <- normalisation_scran(object = object)
    scran <- find_HVGs(object = scran, nfeatures = n.features)
    scran <- Scale_data(object = scran, unwanted.variance = unwanted.variance, scale = scale, center = center)
    altExp(temp, 'scran') <- scran
    print('Norm method: Scran, completed.')
    seuratlog <- normalisation_SeuratLog(object)
    seuratlog <- find_HVGs(object = seuratlog, nfeatures = n.features)
    seuratlog <- Scale_data(object = seuratlog, unwanted.variance = unwanted.variance, scale = scale, center = center)
    altExp(temp, 'seurat_log') <- seuratlog
    print('Norm method: Seurat LogNorm, completed.')
    SCTransform <- normalisation_SCT(object)
    SCTransform <- find_HVGs(object = SCTransform, nfeatures = n.features)
    SCTransform <- Scale_data(object = SCTransform, unwanted.variance = unwanted.variance, scale = scale, center = center)
    altExp(temp, 'SCT') <- SCTransform
    print('Norm method: SCTransform, completed.')
    cpm_norm <- normalisation_cpm(object)
    cpm_norm <- find_HVGs(object = cpm_norm, nfeatures = n.features)
    cpm_norm <- Scale_data(object = cpm_norm, unwanted.variance = unwanted.variance, scale = scale, center = center)
    altExp(temp, 'cpm') <- cpm_norm
    print('Norm method: CPM, completed.')
  }

  ### PCA Analysis ###
  if(method == 'all') {
    names.of.experiments <- altExpNames(temp)
    for(u in names.of.experiments) {
      print(paste0('Calculating PCs for ', u, ' normalisation method.'))
      temp.2 <- altExp(temp, u)
      temp.2 <- performPCA(object = temp.2)
      altExp(temp, u) <- temp.2
      print(paste0('Finished calculating PCs for ', u, ' normalisation method.'))
    }
  } else{
    object <- performPCA(object = object)
  }

  ### benchmarking ###
  if(benchmark==TRUE) {
    if(method == 'all') {
      names.of.experiments <- altExpNames(temp)
      results <- list()
      for(u in names.of.experiments) {
        print(paste0('Calculating benchmarking metrices for ', u, ' normalisation method.'))
        temp.3 <- altExp(temp, u)
        bench.results <- benchmark_norm(object = temp.3, DR.assay = 'pca', distance.method = distance.method)
        altExp(temp, u) <- bench.results
        print(paste0('Finished calculating benchmarking metrices for ', u, ' normalisation method.'))
      }
    } else{
      object <- benchmark_norm(object = object, DR.assay = 'pca', distance.method = distance.method)
      }
    }
  if(method=='all') {
    return(temp)
  }
  else(object)
}

Seurat_cluster <- function(object, multi.method, res=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5), dims=1:12, prune.SNN=0,
                           nn.method='annoy', annoy.metric='euclidean', nn.eps=0.0, ...) {
  mult.object <- object
  if(multi.method==TRUE){
    names.of.experiments <- altExpNames(mult.object)
    for(g in names.of.experiments) {
      print(paste0('Calculating Seurat clusters for ', g, ' normalisation method.'))
      temp <- altExp(mult.object, g)
      orig.names <- colnames(colData(object))
      tmp <- as.Seurat(x = temp, counts='counts', data='data')
      tmp <- FindNeighbors(object = tmp, reduction = 'pca', verbose = TRUE, dims = dims, compute.SNN = TRUE, prune.SNN = prune.SNN,
                           nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
      tmp <- FindClusters(object = tmp, resolution = res, ...)
      new.names <- colnames(tmp@meta.data)
      sep.names <- new.names[!(new.names %in% orig.names)]
      sep.names <-sep.names[1:length(sep.names)-1]
      new.clusters <- tmp@meta.data[,sep.names]
      colnames(new.clusters) <- c('Seurat_res_0.1','Seurat_res_0.2','Seurat_res_0.3','Seurat_res_0.4','Seurat_res_0.5',
                                  'Seurat_res_0.6','Seurat_res_0.7','Seurat_res_0.8','Seurat_res_0.9','Seurat_res_1',
                                  'Seurat_res_1.1','Seurat_res_1.2','Seurat_res_1.3','Seurat_res_1.4','Seurat_res_1.5')
      metadata(temp)$seurat <- new.clusters
      print(metadata(temp))
      altExp(mult.object, g) <- temp
    }
  } else if (multi.method==FALSE) {
    orig.names <- colnames(colData(object))
    tmp <- as.Seurat(x = object, counts='counts', data='data')
    tmp <- FindNeighbors(object = tmp, reduction = 'pca', verbose = TRUE, dims = dims, compute.SNN = TRUE, prune.SNN = prune.SNN,
                         nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
    tmp <- FindClusters(object = tmp, resolution = res, ...)
    new.names <- colnames(tmp@meta.data)
    sep.names <- new.names[!(new.names %in% orig.names)]
    sep.names <-sep.names[1:length(sep.names)-1]
    new.clusters <- tmp@meta.data[,sep.names]
    colnames(new.clusters) <- c('Seurat_res_0.1','Seurat_res_0.2','Seurat_res_0.3','Seurat_res_0.4','Seurat_res_0.5',
                                'Seurat_res_0.6','Seurat_res_0.7','Seurat_res_0.8','Seurat_res_0.9','Seurat_res_1.0',
                                'Seurat_res_1.1','Seurat_res_1.2','Seurat_res_1.3','Seurat_res_1.4','Seurat_res_1.5')
    metadata(object)$clustering$seurat <- new.clusters
    return(object)
  }

  if(multi.method==TRUE) {
    print('Seurat multi-normalisation clustering complete')
    return(mult.object)
  } else if (multi.method==FALSE) {
    print('Seurat single normalisation clustering complete.')
    return(object)
  }
}

SC3_cluster <- function(object, multi.method, ks, n.core=2) {
  multi.object <- object
  if(multi.method==TRUE) {
    names.of.experiments <- altExpNames(multi.object)
    for(n in names.of.experiments) {
      print(paste0('Calculating SC3 clusters for ', n, ' normalisation method.'))
      temp <- altExp(multi.object, n)
      temp.2 <- temp
      logcounts(temp.2) <- assay(temp.2, 'data')
      rowData(temp.2)$feature_symbol <- rownames(temp.2)
      temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
      temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
      if(ks=='estimate') {
        temp.2 <- sc3_estimate_k(temp.2)
        n.of.k <- metadata(temp.2)$sc3$k_estimation
        temp.2 <- sc3_calc_dists(temp.2)
        temp.2 <- sc3_calc_transfs(temp.2)
        temp.2 <- sc3_kmeans(temp.2, ks = n.of.k)
        temp.2 <- sc3_calc_consens(temp.2)
        orig.names <- colnames(colData(temp))
        new.names <- colData(temp.2)
        sep.names <- new.names[!(new.names %in% orig.names)]
        new.clusters <- colData(temp.2)[,sep.names]
        metadata(temp)$SC3_clusters <- new.clusters
        altExp(multi.object, n) <- temp
        print('SC3 multi-normalisation clustering complete')
      } else {
        temp.2 <- sc3_calc_dists(temp.2)
        temp.2 <- sc3_calc_transfs(temp.2)
        temp.2 <- sc3_kmeans(temp.2, ks = ks)
        temp.2 <- sc3_calc_consens(temp.2)
        orig.names <- colnames(colData(temp))
        new.names <- colnames(colData(temp.2))
        sep.names <- new.names[!(new.names %in% orig.names)]
        new.clusters <- colData(temp.2)[,sep.names]
        metadata(temp)$SC3_clusters <- new.clusters
        altExp(multi.object, n) <- temp
        print('SC3 multi-normalisation clustering complete')
      }
    }
  } else {
    sing.temp <- object
    if (ks=='estimate') {
      logcounts(sing.temp) <- assay(sing.temp, 'data')
      rowData(sing.temp)$feature_symbol <- rownames(sing.temp)
      sing.temp <- sing.temp[!duplicated(rowData(sing.temp)$feature_symbol), ]
      sing.temp <- sc3_prepare(sing.temp, gene_filter = FALSE, n_cores = n.core)
      if(ks=='estimate') {
        sing.temp <- sc3_estimate_k(sing.temp)
        n.of.k <- metadata(sing.temp)$sc3$k_estimation
        sing.temp <- sc3_calc_dists(sing.temp)
        sing.temp <- sc3_calc_transfs(sing.temp)
        sing.temp <- sc3_kmeans(sing.temp, ks = n.of.k)
        sing.temp <- sc3_calc_consens(sing.temp)
        orig.names <- colnames(colData(temp))
        new.names <- colData(sing.temp)
        sep.names <- new.names[!(new.names %in% orig.names)]
        new.clusters <- colData(sing.temp)[,sep.names]
        colData(object)$SC3_clusters <- new.clusters
        return(object)
        print('SC3 multi-normalisation clustering complete')
      } else {
        sing.temp <- sc3_calc_dists(sing.temp)
        sing.temp <- sc3_calc_transfs(sing.temp)
        sing.temp <- sc3_kmeans(sing.temp, ks = ks)
        sing.temp <- sc3_calc_consens(sing.temp)
        orig.names <- colnames(colData(temp))
        new.names <- colnames(colData(sing.temp))
        sep.names <- new.names[!(new.names %in% orig.names)]
        new.clusters <- colData(sing.temp)[,sep.names]
        colData(object)$SC3_clusters <- new.clusters
        return(object)
        print('SC3 multi-normalisation clustering complete')
      }
    }
  }
  if (multi.method==TRUE) {
    return(multi.object)
  } else {
    return(object)
  }
}

produce_umap <- function(object, multi.norm, multi.cluster) {
  if(multi.norm==TRUE) {
    multi.object <- object
    if(multi.cluster==TRUE) {

    }
    if(multi.cluster==FALSE) {
      names.of.experiments <- altExpNames(multi.object)
      for(n in names.of.experiments) {
        print(paste0('Calculating umap for ', n, ' normalisation method.'))
        temp <- altExp(multi.object, n)
        c <- uwot::umap(X = reducedDim(temp, 'pca'), verbose = TRUE)
        colnames(c) <- c('umap_1', 'umap_2')
        reducedDim(temp, 'umap') <- c
        altExp(multi.object, n) <- temp
      }
      return(multi.object)
    }
  }
}

benchmarking_clustering <- function(object, components, reduction, dist.method, multi.norm, multi.cluster, sing.clust.method) {
  if(multi.norm==TRUE) {
    multi.object <- object
    names.of.experiments <- altExpNames(multi.object)
    if(multi.cluster==TRUE) {
      for(n in names.of.experiments) {}

    } else if(multi.cluster==FALSE) {
      print(names.of.experiments)
      for(n in names.of.experiments) {
        print(paste0('Calculating benchmark metrics for ', n, ' normalisation method.'))
        temp <- altExp(multi.object, n)
        dims <- reducedDim(temp, reduction)[,1:7]
        dist.matrix <- dist(x = dims, method = dist.method)
        all.clusters <- metadata(temp)[[sing.clust.method]]
        sil.results <- data.frame(average_silhoeutte=NA)
        for (v in colnames(all.clusters)[1:length(colnames(all.clusters))]) {
          print(paste0('Calculating silhouette for ', v))
          tmp <- silhouette(x = as.numeric(x = as.factor(x = all.clusters[,v])), dist = dist.matrix)
          average <- sum(tmp[,3])/length(tmp[,3])
          sil.results[v,] <- average
        }
        sil.results <- sil.results[complete.cases(sil.results),]
        max.AS <- max(sil.results)
        print(max.AS)

        dunn.results <- data.frame(dunn.index=NA)
        for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
          print(paste0('Calculating dunn index for ', p))
          dunn.results[p,] <- dunn(distance = dist.matrix, clusters = as.numeric(x = as.factor(x = all.clusters[,p])))
        }
        dunn.results <- dunn.results[complete.cases(dunn.results),]
        max.dunn <- max(dunn.results)

        conn.results <- data.frame(connectivity=NA)
        for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
          print(paste0('Calculating connectivity for ', p))
          conn.results[p,] <- connectivity(distance = dist.matrix, clusters = all.clusters[,p])
        }
        conn.results <- conn.results[complete.cases(conn.results),]
        max.conn <- max(conn.results)
        results <- cbind(sil.results, dunn.results, conn.results)
        rownames(results) <- colnames(all.clusters)
        metadata(temp)$benchmarking$clustering <- results
        altExp(multi.object, n) <- temp
      }
      return(multi.object)

    }

    }
  }


plot_cluster_dr <- function(object, reduction='umap', pt.size=0.5, multi.norm, multi.cluster, norm.method, cluster.method, cluster.column, seurat_res, kmeans.k) {
  if(multi.norm==TRUE){
    isolated <- altExp(object, as.character(norm.method))
  }
  if(multi.cluster==TRUE){
    all.clusters <- metadata(isolated)[[cluster.method]]
  } else {
    isolated <- object
    all.clusters <- metadata(isolated)[[cluster.method]]
  }
  a <- as.data.frame(reducedDim(isolated, reduction))
  a <- cbind(colnames(isolated), a)
  colnames(a) <- c('Cell.Barcodes', 'reduction_1', 'reduction_2')
  tmp <- as.data.frame(c(colData(isolated), a))
  clust.column <- names(all.clusters[grepl(cluster.column, x = names(all.clusters))])
  cluster <- all.clusters[,clust.column]
  tmp <- cbind(tmp, cluster)
  tmp$cluster <- cluster
  colourCount = length(unique(cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  plot <- ggplot(data = tmp, aes(ident = `Cell.Barcodes`, x = `reduction_1`,
                                y = `reduction_2`, percent.mt = percent.mt,
                                color = cluster)) +
  geom_point(shape = 16, size = pt.size) +
    scale_color_hue() +
    guides(fill=FALSE, alpha=FALSE, size=FALSE) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x=paste0(reduction, '_1'), y=paste0(reduction, '_2'))
    z <- ggplotly(plot)
    return(z)
}
