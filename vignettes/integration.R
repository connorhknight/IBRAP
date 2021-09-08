# IBRAP marrow script

library(IBRAP)

celseq2 <- createIBRAPobject(counts = celseq2,
                             meta.data = metadata_celseq2,
                             original.project = 'celseq2',
                             method.name = 'RAW',
                             min.cells = 3,
                             min.features = 200)

celseq <- createIBRAPobject(counts = celseq, 
                            meta.data = metadata_celseq,
                            original.project = 'celseq',
                            method.name = 'RAW',
                            min.cells = 3,
                            min.features = 200)

pancreas <- merge(x = celseq2, y = celseq)

pancreas <- perform.singleR.annotation(object = pancreas, ref = smartseq2, ref.labels = metadata_smartseq2$celltype)

rm(metadata_celseq, metadata_celseq2, metadata, pancreas.data, celseq, celseq2)

rm(smartseq2, metadata_smartseq2)

fluidigmc1 <- find_percentage_genes(object = fluidigmc1, pattern = '^MT-',
                                    assay = 'RAW', slot = 'counts',
                                    column.name = 'RAW_percent.mt')

pancreas <- find_percentage_genes(object = pancreas, pattern = 'RPL', 
                                  assay = 'RAW', slot = 'counts',
                                  column.name = 'RAW_percent.rp')

plot.QC.vln(object = fluidigmc1, 
            metadata.columns = c('RAW_total.features', 
                                 'RAW_total.counts', 
                                 'RAW_percent.mt'))

ggpubr::annotate_figure(p = p, top = ggpubr::text_grob(label = 'Pancreas', size = 14, family = 'Arial'))

plot.QC.scatter(object = fluidigmc1, 
                x = 'RAW_total.counts', 
                y = 'RAW_total.features', 
                split.by = 'original.project')

plot.QC.scatter(object = fluidigmc1, 
                y = 'RAW_total.counts', 
                x = 'RAW_percent.mt', 
                split.by = 'original.project')

plot.QC.scatter(object = pancreas, 
                y = 'RAW_total.features', 
                x = 'RAW_percent.rp', 
                split.by = 'original.project')

sd.value <- sd(fluidigmc1$RAW_total.features)
med.value <- median(fluidigmc1$RAW_total.features)
max.features <- (sd.value*3)+med.value

fluidigmc1 <- filter_IBRAP(object = fluidigmc1, 
                           RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 25)

pancreas <- add.cell.cycle(object = pancreas, 
                           assay = 'RAW', 
                           slot = 'counts', 
                           transform = TRUE)

pancreas <- add.feature.score(object = pancreas, 
                              assay = 'RAW', 
                              slot = 'counts',
                              transform = TRUE, 
                              features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
                                           'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
                                           'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
                                           'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
                                           'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
                              column.name = 'StressScore')

fluidigmc1 <- perform.sct(object = fluidigmc1, 
                          assay = 'RAW', 
                          slot = 'counts')

pancreas <- perform.scran(object = pancreas, 
                          assay = 'RAW', 
                          slot = 'counts', 
                          vars.to.regress = 'RAW_total.counts', do.scale = T)

pancreas <- perform.scanpy(object = pancreas, 
                           vars.to.regress = 'RAW_total.counts', do.scale = T)

# pancreas <- perform.tpm.normalisation(object = pancreas, 
#                                       vars.to.regress = 'RAW_total.counts', do.scale = T)

pancreas <- perform.pca(object = pancreas, 
                        assay = c('SCT', 'SCRAN', 'SCANPY'), 
                        n.pcs = 50, reduction.save = 'pca')



pancreas <- perform.bbknn(object = pancreas, 
                          assay = c('SCT', 'SCANPY', 'SCRAN'), 
                          reduction = c('pca'),
                          batch = 'tech', generate.diffmap = T)

pancreas <- perform.harmony(object = pancreas, 
                            assay = c('SCRAN', 'SCT', 'SCANPY'), 
                            vars.use = 'original.project', 
                            reduction = c('pca'),  
                            max.iter.harmony = 100,
                            dims.use = list(NULL))

pancreas <- perform.scanorama(object = pancreas, 
                              assay = c('SCT', 'SCRAN', 'SCANPY'), 
                              slot = 'norm.scaled', 
                              split.by = 'original.project', 
                              n.dims = 50)



pancreas <- perform.nn.v1(object = pancreas, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                          reduction = c('pca_harmony','scanorama')
                          , dims = list(0,0), generate.diffmap = T)

pancreas <- perform.nn.v1(object = pancreas, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                          reduction = c('pca_bbknn_bbknn:diffmap','pca_harmony_nn.v1:diffmap', 'scanorama_nn.v1:diffmap'), 
                          dims = list(0,0,0))

pancreas <- perform.nn.v2(object = pancreas, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                          reduction = c('pca_harmony','scanorama','pca_bbknn_bbknn:diffmap',
                                        'pca_harmony_nn.v1:diffmap', 'scanorama_nn.v1:diffmap'), 
                          dims = list(0,0,0,0,0))

pancreas <- perform.umap(object = pancreas, 
                         assay = c('SCT', 'SCRAN', 'SCANPY'), 
                         reduction = c('pca_harmony', 'scanorama', 'pca_bbknn_bbknn:diffmap', 'pca_harmony_nn.v1:diffmap', 'scanorama_nn.v1:diffmap'), 
                         n_components = 3, 
                         n.dims = list(1:50, 1:50, NULL, NULL, NULL))

pancreas <- perform.umap(object = pancreas, assay = c('SCT', 'SCRAN', 'SCANPY'), graph = 'pca_bbknn_bbknn')

pancreas <- perform.seurat.cluster(object = pancreas, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                   neighbours = c("pca_bbknn","pca_harmony_seurat","scanorama_seurat","pca_harmony_scanpy","scanorama_scanpy",
                                                  "pca_bbknn_diffmap_scanpy","pca_harmony_diffmap_scanpy","scanorama_diffmap_scanpy"), 
                                   algorithm = 1, cluster.df.name = c("pca_bbknn_louvain","pca_harmony_seurat_louvain","scanorama_seurat_louvain",
                                                                      "pca_harmony_scanpy_louvain","scanorama_scanpy_louvain",
                                                                      "pca_bbknn_diffmap_scanpy_louvain","pca_harmony_diffmap_scanpy_louvain",
                                                                      "scanorama_diffmap_scanpy_louvain"))

pancreas_bench <- benchmark.clustering(object = pancreas, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                       clustering = c("pca_harmony_seurat_louvain",
                                                      "scanorama_seurat_louvain",
                                                      "pca_harmony_scanpy_louvain",
                                                      "scanorama_scanpy_louvain",
                                                      "pca_bbknn_diffmap_scanpy_louvain",
                                                      "pca_harmony_diffmap_scanpy_louvain",
                                                      "scanorama_diffmap_scanpy_louvain"), 
                                       reduction = c('pca_harmony_umap',
                                                     'scanorama_umap',
                                                     'pca_harmony_umap',
                                                     'scanorama_umap',
                                                     'pca_bbknn_diffmap_umap',
                                                     'pca_harmony_diffmap_umap',
                                                     'scanorama_diffmap_umap'
                                       ), 
                                       n.dims = 1:2, ground.truth = pancreas@sample_metadata$celltype)

pancreas_bench <- benchmark.clustering(object = pancreas_bench, assay = c('SCT', 'SCRAN', 'SCANPY'), 
                                       clustering = c("pca_bbknn_louvain"), 
                                       reduction = c('pca_bbknn_umap'
                                       ), 
                                       n.dims = 1:2, ground.truth = pancreas@sample_metadata$celltype)

SCT_DE <- perform.seurat.diffexp.all(object = pancreas, assay = 'SCT', test = 'MAST', identity = pancreas@sample_metadata$celltype, latent.vars = 'original.project')
SCRAN_DE <- perform.seurat.diffexp.all(object = pancreas, assay = 'SCT', test = 'MAST', identity = pancreas@sample_metadata$celltype, latent.vars = 'original.project')
SCANPY_DE <- perform.seurat.diffexp.all(object = pancreas, assay = 'SCT', test = 'MAST', identity = pancreas@sample_metadata$celltype, latent.vars = 'original.project')

SCT_DE_GO <- perform.GO.enrichment(result = SCT_DE)
SCRAN_DE_GO <- perform.GO.enrichment(result = SCRAN_DE)
SCANPY_DE_GO <- perform.GO.enrichment(result = SCANPY_DE)

plot.GO.output(result = SCT_DE_GO) + ggplot2::ggtitle(label = 'SCT')
plot.GO.output(result = SCRAN_DE_GO) + ggplot2::ggtitle(label = 'SCRAN')
plot.GO.output(result = SCANPY_DE_GO) + ggplot2::ggtitle(label = 'SCANPY')




