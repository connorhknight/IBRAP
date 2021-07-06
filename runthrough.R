
library(IBRAP)

sample_1 <- read.csv('/Users/knight05/Downloads/SP3tFLBcells_dfall_decontx_ec.csv', 
                     header = T, row.names = 1, sep = ',')

sample <- createIBRAPobject(counts = as.matrix(sample), original.project = 'FL', 
                            min.cells = 3, min.features = 200)


sample <- find_percentage_genes(object = sample, pattern = '^MT-',
                                assay = 'RAW', slot = 'counts',
                                column.name = 'RAW_percent.mt')

plot.QC.vln(object = sample, 
            metadata.columns = c('RAW_total.features', 
                                 'RAW_total.counts', 
                                 'RAW_percent.mt'))

plot.QC.scatter(object = sample, 
                x = 'RAW_total.counts', 
                y = 'RAW_total.features', 
                split.by = 'original.project')

plot.QC.scatter(object = sample, 
                y = 'RAW_total.counts', 
                x = 'RAW_percent.mt', 
                split.by = 'original.project')

plot.QC.scatter(object = sample, 
                y = 'RAW_total.features', 
                x = 'RAW_percent.mt', 
                split.by = 'original.project')

sd.value <- sd(sample$RAW_total.features)
med.value <- median(sample$RAW_total.features)
max.features <- (sd.value*3)+med.value

sample <- filter_IBRAP(object = sample, 
                       RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 10)

sample <- add.cell.cycle(object = sample, 
                         assay = 'RAW', 
                         slot = 'counts', 
                         transform = TRUE)

sample <- add.feature.score(object = sample, 
                            assay = 'RAW', 
                            slot = 'counts',
                            transform = TRUE, 
                            features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
                                         'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
                                         'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
                                         'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
                                         'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
                            column.name = 'StressScore')

sample <- perform.sct.normalisation(object = sample, 
                                    assay = 'RAW', 
                                    slot = 'counts', vars.to.regress = 'RAW_percent.mt')

sample <- perform.scran.normalisation(object = sample, 
                                      assay = 'RAW', 
                                      slot = 'counts', 
                                      vars.to.regress = c('RAW_total.counts', 'RAW_percent.mt'), do.scale = T)

sample <- perform.scanpy.normalisation(object = sample, 
                                       vars.to.regress = c('RAW_total.counts', 'RAW_percent.mt'), 
                                       do.scale = T, log1 = F)

# sample <- remove.hvgs(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'), hvgs.omit = list(unwanted,
#                                                                                              unwanted,
#                                                                                              unwanted))    

sample <- perform.pca(object = sample, 
                      assay = c('SCT', 'SCRAN', 'SCANPY'), 
                      n.pcs = 1:50, 
                      reduction.save = 'pca')

sample <- perform.seurat.neighbours(object = sample, 
                                    assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                    reduction = c('pca'), 
                                    dims = list(20), 
                                    neighbour.name = c('pca_seurat'))

sample <- perform.scanpy.neighbours(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                    reduction = c('pca'), 
                                    neighbour.name = c('pca_scanpy'), 
                                    ims = list(20), 
                                    generate.diffmap = T,
                                    diffmap.name = c('pca_diffmap'))

sample <- perform.scanpy.neighbours(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                    reduction = c('pca_diffmap'), 
                                    neighbour.name = c('pca_diffmap_scanpy'), 
                                    dims = list(20))

sample <- perform.seurat.neighbours(object = sample, 
                                    assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                    reduction = c('pca_diffmap'), 
                                    dims = list(0), 
                                    neighbour.name = c('pca_diffmap_seurat'))

sample <- perform.umap(object = sample, 
                       
                       assay = c('SCT', 
                                 'SCRAN', 
                                 'SCANPY'), 
                       
                       reduction = c('pca', 
                                     'pca_diffmap'), 
                       
                       reduction.save = c('pca_umap', 
                                          'pca_diffmap_umap'), 
                       
                       n_components = 2, 
                       
                       n.dims = list(1:20, 
                                     NULL))

sample <- perform.seurat.cluster(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                 neighbours = c("pca_seurat",
                                                "pca_scanpy",
                                                "pca_diffmap_seurat",
                                                "pca_diffmap_scanpy"), 
                                 algorithm = 1, 
                                 cluster.df.name = c("pca_seurat_louvain",
                                                     "pca_scanpy_louvain",
                                                     "pca_diffmap_seurat_louvain",
                                                     "pca_diffmap_scanpy_louvain"))

plot.reduced.dim(object = sample, 
                 reduction = 'pca_umap', 
                 assay = 'SCANPY', 
                 clust.method = 'pca_seurat_louvain', 
                 column = 'neighbourhood_graph_res.0.6')

tmp <- perform.slingshot.trajectory(object = sample, 
                                    reduction = 'pca_umap', 
                                    assay = 'SCT', 
                                    clust.method = 'pca_seurat_louvain', 
                                    column = 'neighbourhood_graph_res.0.6')

plot.slingshot(result = tmp, 
               relevant = F, 
               Pseudotime = T)
