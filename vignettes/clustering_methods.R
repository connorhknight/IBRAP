library(IBRAP)

samp <- Read10X_output(directory = '/Users/knight05/Raw_Data/Database_samples/healthy_references/BMMC_atlas/marrow_A', matrix.file = 'matrix.mtx', 
                       genes.file = 'genes.tsv', barcodes.file = 'barcodes.tsv')

# raw cells = 2994

samp <- createIBRAPobject(counts = samp, original.project = 'BMMC', min.cells = 3, min.features = 200)

samp <- find_percentage_genes(object = samp)

plot.QC.vln(object = samp, 
            metadata.columns = c('RAW_total.features', 
                                 'RAW_total.counts', 
                                 'RAW_percent.mt'))

plot.QC.scatter(object = samp, 
                y = 'RAW_total.counts', 
                x = 'RAW_percent.mt', 
                split.by = 'original.project')

sd.value <- sd(samp$RAW_total.features)
med.value <- median(samp$RAW_total.features)
max.features <- (sd.value*3)+med.value

samp <- filter_IBRAP(object = samp, 
                     RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 8)

samp <- perform.sct(object = samp, 
                    assay = 'RAW', 
                    slot = 'counts')

samp <- perform.scran(object = samp, 
                      assay = 'RAW', 
                      slot = 'counts', 
                      vars.to.regress = 'RAW_total.counts', do.scale = T)

samp <- perform.scanpy(object = samp, 
                       vars.to.regress = 'RAW_total.counts', do.scale = T)

samp <- perform.pca(object = samp, 
                    assay = c('SCT', 'SCRAN', 'SCANPY'), 
                    n.pcs = 50, reduction.save = 'pca')



