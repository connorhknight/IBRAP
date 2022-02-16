library('SeuratData')
library('Seurat')
library('IBRAP')

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

pancreas <- perform.sct(object = pancreas, verbose = T)

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
                x = 'RAW_total.counts', 
                y = 'RAW_total.features', 
                split.by = 'Phase')

pancreas <- perform.pca(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), print.variance = T)

pancreas <- perform.bbknn(object = pancreas, assay = c('SCT','SCRAN','SCANPY','TPM'), reduction = 'pca', batch = 'original.project', verbose = T)

