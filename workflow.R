library('SeuratData')
library('Seurat')
library('IBRAP')

setwd( "/Users/knight05/work/Results/scRNA-seq/IBRAP_development/IBRAPwithdecontX/R")
files.sources = list.files()
sapply(files.sources, source)

options(future.globals.maxSize = 4000 * 1024^2)

panc8 <- LoadData('panc8')

smartseq2.counts <- panc8@assays$RNA@counts[,panc8@meta.data$dataset=='smartseq2']
smartseq2.meta <- panc8@meta.data[panc8@meta.data$dataset=='smartseq2',]

celseq2.counts <- panc8@assays$RNA@counts[,panc8@meta.data$dataset=='celseq2']
celseq2.meta <- panc8@meta.data[panc8@meta.data$dataset=='celseq2',]

rm(panc8)

smartseq2 <- createIBRAPobject(counts = smartseq2.counts, 
                               original.project = 'smartseq2', 
                               add.suffix = T, 
                               meta.data = smartseq2.meta, min.cells = 3, min.features = 200, verbose = T)

celseq2 <- createIBRAPobject(counts = celseq2.counts, 
                             original.project = 'celseq2', 
                             add.suffix = T, 
                             meta.data = celseq2.meta, min.cells = 3, min.features = 200, verbose = T)

pancreas <- merge(x = smartseq2, celseq2)

find_percentage_genes(object = smartseq2, 
                      pattern = '^MT-', 
                      assay = 'RAW', 
                      slot = 'counts', 
                      column.name = 'mito.percent')

find_percentage_genes(object = smartseq2, 
                      pattern = '^MT-', 
                      assay = 'RAW', 
                      slot = 'counts', 
                      column.name = 'mito.percent', 
                      verbose = T)

pancreas <- add.cell.cycle(object = pancreas, 
                           assay = 'RAW', 
                           slot = 'counts', 
                           transform = T)

pancreas <- add.cell.cycle(object = pancreas, 
                           assay = 'RAW', 
                           slot = 'counts', 
                           transform = T, 
                           verbose = T)

pancreas <- perform.sct(object = pancreas, verbose = T, conserve.memory=T)

pancreas <- perform.scran(object = pancreas, vars.to.regress = 'RAW_total.counts', verbose = T)

pancreas <- perform.scanpy(object = pancreas, vars.to.regress = 'RAW_total.counts', verbose = T)

pancreas <- perform.tpm(object = pancreas, vars.to.regress = 'RAW_total.counts', verbose = T)

plot.QC.vln(object = pancreas, 
            metadata.columns = c("RAW_total.counts",
                                 "SCT_total.counts",
                                 "SCRAN_total.counts",
                                 "SCANPY_total.counts", 
                                 "TPM_total.counts"))

plot.QC.scatter(object = pancreas, 
                x = 'SCT_total.counts', 
                y = 'SCT_total.features', 
                split.by = 'Phase')

pancreas <- perform.pca(object = pancreas, assay = c('SCT','SCRAN','SCANPY'), print.variance = F)

# pancreas <- perform.bbknn(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), 
#                           reduction = 'PCA', batch = 'original.project', verbose = T)

pancreas <- perform.scanorama(object = pancreas, assay = c('SCT','SCRAN','SCANPY'), 
                              batch = 'original.project', verbose = T)

pancreas <- perform.harmony(object = pancreas, assay = c('SCT','SCRAN','SCANPY'), 
                            reduction = 'PCA', batch = 'original.project', verbose = T)

# list.of.objects <- splitIBRAP(object = pancreas, split.by = 'original.project')
# 
# list.of.objects <- lapply(X = list.of.objects, FUN = 'perform.sct', verbose = T)
# 
# list.of.objects <- lapply(X = list.of.objects, FUN = 'perform.scran', vars.to.regress = 'RAW_total.counts', verbose = T)
# 
# list.of.objects <- lapply(X = list.of.objects, FUN = 'perform.scanpy', vars.to.regress = 'RAW_total.counts', verbose = T)
# 
# list.of.objects <- lapply(X = list.of.objects, FUN = 'perform.tpm', vars.to.regress = 'RAW_total.counts', verbose = T)
# 
# tmp <- perform.seurat.integration(object = pancreas, 
#                                   object.list = list.of.objects, 
#                                   assay = c('SCT','SCRAN','SCANPY','TPM'), 
#                                   nfeatures = 1500, print.variance = T, 
#                                   verbose = T)

pancreas <- perform.nn(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), 
                       reduction = c('SCANORAMA','PCA_HARMONY','PCA'))

pancreas <- perform.graph.cluster(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), 
                                  neighbours = c("PCA_BBKNN_BBKNN","SCANORAMA_NN","PCA_HARMONY_NN",'PCA_NN'), 
                                  algorithm = 1, verbose = T)

pancreas <- perform.graph.cluster(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), 
                                  neighbours = c("PCA_BBKNN_BBKNN","SCANORAMA_NN","PCA_HARMONY_NN", "PCA_NN"), 
                                  algorithm = 2, verbose = T)

pancreas <- perform.umap(object = pancreas, assay = c('SCT','SCRAN','SCANPY'), 
                         reduction = c('SCANORAMA','PCA_HARMONY','PCA'), verbose = F)

pancreas <- perform.umap(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), 
                         graph = c('PCA_BBKNN_BBKNN'), verbose = F)

tmp <- benchmark.clustering(object = pancreas, 
                            assay = c('SCT','SCRAN','SCANPY','TPM'), 
                            clustering = c("PCA_NN:LOUVAIN",
                                           "PCA_BBKNN_BBKNN:LOUVAIN",
                                           "SCANORAMA_NN:LOUVAIN",
                                           "PCA_HARMONY_NN:LOUVAIN"),
                            reduction = c('PCA_UMAP',
                                          'PCA_BBKNN_BBKNN:UMAP',
                                          'SCANORAMA_UMAP',
                                          "PCA_HARMONY_UMAP"), 
                            ground.truth.column = 'celltype', 
                            verbose = T)


