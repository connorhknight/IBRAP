################################################################################################################################################################
#IBRAP integration
library(IBRAP)
library(reticulate)

#Check to see if all packages are installed.
prepare.reticulate()


#Load Sample
sample_1 <- read.csv('R33651_ec_CRv4.csv',
                     header = T, row.names = 1, sep = ',')
sample_2 <- read.csv('R34021_ec_CRv4.csv',
                     header = T, row.names = 1, sep = ',')
sample_3 <- read.csv('R34382_ec_CRv4.csv',
                     header = T, row.names = 1, sep = ',')




#Remove Doublets
sample_1 <- perform.scrublet(counts = as.matrix(sample_1))
sample_2 <- perform.scrublet(counts = as.matrix(sample_2))
sample_3 <- perform.scrublet(counts = as.matrix(sample_3))




#Ambient RNA decontamination (optional)
sample_1 <- perform.decontX(counts = sample_1)
sample_2 <- perform.decontX(counts = sample_2)
sample_3 <- perform.decontX(counts = sample_3)




#Create Objects
sample_1_obj <- createIBRAPobject(counts = sample_1,
                             original.project = 'R33651',
                             method.name = 'RAW',
                             min.cells = 3,
                             min.features = 200)

sample_2_obj <- createIBRAPobject(counts = sample_2,
                                  original.project = 'R34021',
                                  method.name = 'RAW',
                                  min.cells = 3,
                                  min.features = 200)

sample_3_obj <- createIBRAPobject(counts = sample_3,
                                  original.project = 'R34382',
                                  method.name = 'RAW',
                                  min.cells = 3,
                                  min.features = 200)




#Merge Objects
rRLN <- merge(x = sample_1_obj, y = c(sample_2_obj, sample_3_obj))
rRLN





#Plots QC
rRLN <- find_percentage_genes(object = rRLN, pattern = '^MT-',
                                   assay = 'RAW', slot = 'counts',
                                   column.name = 'RAW_percent.mt')
plot.QC.vln(object = rRLN,
             metadata.columns = c('RAW_total.features',
                                  'RAW_total.counts',
                                  'RAW_percent.mt'))
#ggpubr::annotate_figure(p = rRLN, top = ggpubr::text_grob(label = 'rRLN', size = 14, family = 'Arial'))
plot.QC.scatter(object = rRLN,
                 x = 'RAW_total.counts',
                 y = 'RAW_total.features',
                 split.by = 'original.project')
plot.QC.scatter(object = rRLN,
                 y = 'RAW_total.counts',
                 x = 'RAW_percent.mt',
                 split.by = 'original.project')




#Filter outliers
sd.value <- sd(rRLN$RAW_total.features)
med.value <- median(rRLN$RAW_total.features)
max.features <- (sd.value*3)+med.value
rRLN
rRLN <- filter_IBRAP(object = rRLN,
                       RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 10)
rRLN





#Run Normalisation, scaling, variance stabilisation, highly variable genes identification
rRLN <- perform.sct(object = rRLN,
                      assay = 'RAW',
                      slot = 'counts',
                      vars.to.regress = 'RAW_percent.mt',
                      variable.features.n = 1500)

rRLN <- perform.scran(object = rRLN,
                        assay = 'RAW',
                        slot = 'counts',
                        vars.to.regress = c('RAW_total.counts', 'RAW_percent.mt'), do.scale = T)

rRLN <- perform.scanpy(object = rRLN,
                         vars.to.regress = c('RAW_total.counts', 'RAW_percent.mt'))





#Remove any unwanted genes from highly variable gene list (optional)
unwanted2 <- as.character(rownames(rRLN)[grepl("^MT-" ,x=rownames(rRLN))])
unwanted <- c('IGK','IGL','TCRA-VDJ','TCRB-VDJ','TCRG-VDJ','TCRD-VDJ','MT-ATP6','MT-ATP8','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND3','MT-ND4','MT-ND4L','MT-ND6','MT-ND2','MT-ND5','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L7','MTRNR2L5','MTRNR2L8','MTRNR2L4','MTRNR2L1','MTRNR2L3','MTRNR2L10',unwanted2)
unwanted
rRLN <- remove.hvgs(object = rRLN, assay = c('SCT', 'SCRAN', 'SCANPY'), hvgs.omit = list(unwanted,unwanted,unwanted))





#Run PCA on unintegrated data
rRLN <- perform.pca(object = rRLN,
                        assay = c('SCT', 'SCRAN', 'SCANPY'),
                        n.pcs = 50, reduction.save = 'pca')




#Run UMAP on unintegrated data
rRLN <- perform.umap(object = rRLN,
                         assay = c('SCT', 'SCRAN','SCANPY'),
                         reduction = c('pca'),
                         n_components = 2,
                         n.dims = list(20))





#Plot Unintegrated data for exploration
plot1 <- plot.reduced.dim(object = rRLN, reduction = 'pca_umap', assay = 'SCT',
                          clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
ggplot2::ggtitle('SCT') +
ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


plot2 <- plot.reduced.dim(object = rRLN, reduction = 'pca_umap', assay = 'SCRAN',
                          clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
ggplot2::ggtitle('SCRAN') +
ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


plot3 <- plot.reduced.dim(object = rRLN, reduction = 'pca_umap', assay = 'SCANPY',
                          clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCANPY') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot1 + plot2 + plot3





#Perform BBKNN Integration
rRLN <- perform.bbknn(object = rRLN,
                          assay = c('SCT', 'SCANPY', 'SCRAN'),
                          reduction = c('pca'), batch = 'original.project', n_pcs = list(0))





#Perform Harmony Integration
rRLN <- perform.harmony(object = rRLN,
                            assay = c('SCRAN', 'SCT', 'SCANPY'),
                            vars.use = 'original.project',
                            reduction = c('pca'),
                            max.iter.harmony = 100,
                            dims.use = list(0))





#Perform Scanorama Integration
rRLN <- perform.scanorama(object = rRLN,
                              assay = c('SCT', 'SCRAN', 'SCANPY'),
                              slot = 'norm.scaled',
                              split.by = 'original.project',
                              n.dims = 50)




#Perform Seurat Integration
rRLN <- perform.seurat.integration(object = rRLN,
                                         assay = c('SCT','SCRAN','SCANPY'),
                                         normalisation.method = c('perform.sct','perform.scran','perform.scanpy'),
                                         batch = 'original.project')




#Run neighbour finding for pca
rRLN <- perform.nn(object = rRLN, assay = c('SCT', 'SCRAN', 'SCANPY'),
                          reduction = c('pca_harmony','scanorama', 'seurat'),
                          dims = list(20,20,20))





#Integration UMAP
rRLN <- perform.umap(object = rRLN,
                         assay = c('SCT', 'SCRAN', 'SCANPY'),
                         reduction = c('pca_harmony', 'scanorama', 'seurat'),
                         n_components = 2,
                         n.dims = list(0, 0, 0))



#BBKNN UMAP
rRLN <- perform.umap(object = rRLN, assay = c('SCT', 'SCRAN', 'SCANPY'), graph = 'pca_bbknn_bbknn')





#Plot Batch correction results
plot.list <- list()

plot.list[[1]] <- plot.reduced.dim(object = rRLN, reduction = 'pca_harmony_umap', assay = 'SCT',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCT_harmony') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[2]] <- plot.reduced.dim(object = rRLN, reduction = 'pca_harmony_umap', assay = 'SCRAN',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCRAN_harmony') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[3]] <- plot.reduced.dim(object = rRLN, reduction = 'pca_harmony_umap', assay = 'SCANPY',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCANPY_harmony') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[4]] <- plot.reduced.dim(object = rRLN, reduction = 'scanorama_umap', assay = 'SCT',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCT_scanorama') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[5]] <- plot.reduced.dim(object = rRLN, reduction = 'scanorama_umap', assay = 'SCRAN',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCRAN_scanorama') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[6]] <- plot.reduced.dim(object = rRLN, reduction = 'scanorama_umap', assay = 'SCANPY',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCANPY_scanorama') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[7]] <- plot.reduced.dim(object = rRLN, reduction = 'pca_bbknn_bbknn:umap', assay = 'SCT',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCT_bbknn') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[8]] <- plot.reduced.dim(object = rRLN, reduction = 'pca_bbknn_bbknn:umap', assay = 'SCRAN',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCRAN_bbknn') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[9]] <- plot.reduced.dim(object = rRLN, reduction = 'pca_bbknn_bbknn:umap', assay = 'SCANPY',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCANPY_bbknn') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[10]] <- plot.reduced.dim(object = rRLN, reduction = 'seurat_umap', assay = 'SCT',
                                   clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCT_seurat') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[11]] <- plot.reduced.dim(object = rRLN, reduction = 'seurat_umap', assay = 'SCRAN',
                                    clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCRAN_seurat') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

plot.list[[12]] <- plot.reduced.dim(object = rRLN, reduction = 'seurat_umap', assay = 'SCANPY',
                                    clust.method = 'metadata', column = 'original.project', pt.size = 0.1) +
  ggplot2::ggtitle('SCANPY_seurat') +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

egg::ggarrange(plots = plot.list, nrow = 4, ncol = 3)




#Featureplots
plot.list <- list()

plot.list[[1]] <- plot.features(object = rRLN, assay = 'SCT', slot = 'normalised',
              reduction = 'pca_harmony_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[2]] <- plot.features(object = rRLN, assay = 'SCRAN', slot = 'normalised',
              reduction = 'pca_harmony_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[3]] <- plot.features(object = rRLN, assay = 'SCANPY', slot = 'normalised',
              reduction = 'pca_harmony_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[4]]<- plot.features(object = rRLN, assay = 'SCT', slot = 'normalised',
              reduction = 'scanorama_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[5]] <- plot.features(object = rRLN, assay = 'SCRAN', slot = 'normalised',
              reduction = 'scanorama_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[6]] <- plot.features(object = rRLN, assay = 'SCANPY', slot = 'normalised',
              reduction = 'scanorama_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[7]] <- plot.features(object = rRLN, assay = 'SCT', slot = 'normalised',
              reduction = 'seurat_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[8]] <- plot.features(object = rRLN, assay = 'SCRAN', slot = 'normalised',
              reduction = 'seurat_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)

plot.list[[9]] <- plot.features(object = rRLN, assay = 'SCANPY', slot = 'normalised',
              reduction = 'seurat_umap', features = c('CD79A','CD79B',"CD3G","COL1A1"), pt_size = 0.6)




#Run clustering
rRLN <- perform.graph.cluster(object = rRLN, assay = c('SCT', 'SCRAN', 'SCANPY'),
                                  neighbours = c("pca_harmony_nn","scanorama_nn","seurat_nn","pca_bbknn_bbknn"),
                                  algorithm = 1)
#Plot cluster membership
plot.reduced.dim(object = rRLN, reduction = 'pca_harmony_umap', assay = 'SCRAN',
                 clust.method = 'pca_harmony_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'scanorama_umap', assay = 'SCRAN',
                 clust.method = 'scanorama_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'pca_bbknn_bbknn:umap', assay = 'SCRAN',
                 clust.method = 'pca_bbknn_bbknn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'seurat_umap', assay = 'SCRAN',
                 clust.method = 'seurat_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")

plot.reduced.dim(object = rRLN, reduction = 'pca_harmony_umap', assay = 'SCT',
                 clust.method = 'pca_harmony_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'scanorama_umap', assay = 'SCT',
                 clust.method = 'scanorama_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'pca_bbknn_bbknn:umap', assay = 'SCT',
                 clust.method = 'pca_bbknn_bbknn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'seurat_umap', assay = 'SCT',
                 clust.method = 'seurat_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")

plot.reduced.dim(object = rRLN, reduction = 'pca_harmony_umap', assay = 'SCANPY',
                 clust.method = 'pca_harmony_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'scanorama_umap', assay = 'SCANPY',
                 clust.method = 'scanorama_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'pca_bbknn_bbknn:umap', assay = 'SCANPY',
                 clust.method = 'pca_bbknn_bbknn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")
plot.reduced.dim(object = rRLN, reduction = 'seurat_umap', assay = 'SCANPY',
                 clust.method = 'seurat_nn:louvain', pt.size = 0.1, column = "neighbourhood_graph_res.0.1")




#Benchmarking for clustering
#rRLN_bench <- benchmark.clustering(object = rRLN, assay = c('SCT', 'SCRAN', 'SCANPY'), clustering = c("pca_harmony_nn.v1:louvain","pca_harmony_nn.v2:louvain"), reduction = c('pca_harmony_umap','pca_harmony_umap'), n.dims = 1:2)




#Benchmarking for batch correction
rRLN <- benchmark.intergation(object = rRLN, batch = 'original.project', assays = c('SCT','SCRAN','SCANPY'), reduction = c('pca_umap', 'pca_harmony_umap', 'scanorama_umap','seurat_umap'), result.names = c('uncorrected', 'harmony', 'scanorama', 'seurat'), n.components = 2)
plot.integration.benchmarking(object = rRLN, assay = c('SCT','SCRAN','SCANPY'))




#Perform Differential Expression
SCT_DE <- perform.diffexp.all(object = rRLN, assay = 'SCT', test = 'wilcox', clust.method = 'pca_harmony_nn:louvain', column = "neighbourhood_graph_res.0.1")



#Perform Enrichment Analysis
SCT_DE_GO <- perform.GO.enrichment(result = SCT_DE)
SCRAN_DE_GO <- perform.GO.enrichment(result = SCRAN_DE)
SCANPY_DE_GO <- perform.GO.enrichment(result = SCANPY_DE)
plot.GO.output(result = SCT_DE_GO) + ggplot2::ggtitle(label = 'SCT')
plot.GO.output(result = SCRAN_DE_GO) + ggplot2::ggtitle(label = 'SCRAN')
plot.GO.output(result = SCANPY_DE_GO) + ggplot2::ggtitle(label = 'SCANPY')
