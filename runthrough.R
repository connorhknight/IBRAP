##IBRAP
library(IBRAP)
library(reticulate)

#Check to see if all packages are installed.
#prepare.reticulate()

#Load Sample
sample_1 <- read.csv('TGFB.csv',
                     header = T, row.names = 1, sep = ',')

#Remove Doublets
sample_1 <- perform.scrublet(counts = as.matrix(sample_1))

#Ambient RNA decontamination (optional)
sample_1 <- perform.decontX(counts = sample_1)

#Create Seurat Object
sample <- createIBRAPobject(counts = as.matrix(sample_1), original.project = 'FL',
                            min.cells = 3, min.features = 200)

#Find each cell proportional of expression MT genes
sample <- find_percentage_genes(object = sample, pattern = '^MT-',
                                assay = 'RAW', slot = 'counts',
                                column.name = 'RAW_percent.mt')
#QC plots
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

#Find possible multiplets
sd.value <- sd(sample$RAW_total.features)
med.value <- median(sample$RAW_total.features)
max.features <- (sd.value*3)+med.value

#Filter outlier cells
sample <- filter_IBRAP(object = sample,
                       RAW_total.features < max.features & RAW_total.counts > 200 & RAW_percent.mt < 10)

#Add cell cycle score (optional)
sample <- add.cell.cycle(object = sample,
                         assay = 'RAW',
                         slot = 'counts',
                         transform = TRUE)

#Add any geneset of interest score (optional)
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

#Run Normalisation, scaling, variance stabilisation, highly variable genes identification
sample <- perform.sct.normalisation(object = sample,
                                    assay = 'RAW',
                                    slot = 'counts',
                                    vars.to.regress = 'RAW_percent.mt',
                                    variable.features.n = 1500)

sample <- perform.scran.normalisation(object = sample,
                                      assay = 'RAW',
                                      slot = 'counts',
                                      vars.to.regress = c('RAW_total.counts', 'RAW_percent.mt'), do.scale = T)

sample <- perform.scanpy.normalisation(object = sample,
                                       vars.to.regress = c('RAW_total.counts', 'RAW_percent.mt'),
                                       do.scale = T, log1 = F)

#Remove any unwanted genes from highly variable gene list (optional)
unwanted2 <- as.character(rownames(sample)[grepl("^MT-" ,x=rownames(sample))])
unwanted <- c('IGK','IGL',unwanted2)
sample <- remove.hvgs(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'), hvgs.omit = list(unwanted,unwanted,unwanted))

#Run PCA
sample <- perform.pca(object = sample,
                      assay = c('SCT', 'SCRAN', 'SCANPY'),
                      n.pcs = 1:50,
                      reduction.save = 'pca')

#Run neighbour finding
sample <- perform.seurat.neighbours(object = sample,
                                    assay = c('SCT', 'SCRAN', 'SCANPY'),
                                    reduction = c('pca'),
                                    dims = list(0),
                                    neighbour.name = c('pca_seurat'))

sample <- perform.scanpy.neighbours(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'),
                                    reduction = c('pca'),
                                    neighbour.name = c('pca_scanpy'),
                                    dims = list(0),
                                    generate.diffmap = T,
                                    diffmap.name = c('pca_diffmap'))

sample <- perform.scanpy.neighbours(object = sample, assay = c('SCT', 'SCRAN', 'SCANPY'),
                                    reduction = c('pca_diffmap'),
                                    neighbour.name = c('pca_diffmap_scanpy'),
                                    dims = list(0))

sample <- perform.seurat.neighbours(object = sample,
                                    assay = c('SCT', 'SCRAN', 'SCANPY'),
                                    reduction = c('pca_diffmap'),
                                    dims = list(0),
                                    neighbour.name = c('pca_diffmap_seurat'))

#Run UMAP
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

#Run clustering
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
plot.list <- list()

plot.list[[1]] <- plot.reduced.dim(object = sample,
                                   reduction = 'pca_umap',
                                   assay = 'SCT',
                                   clust.method = 'pca_seurat_louvain',
                                   column = 'neighbourhood_graph_res.0.6', pt.size = 0.6) +
  ggplot2::theme(legend.position = 'none')

plot.list[[1]]

plot.list[[2]] <- plot.reduced.dim(object = sample,
                                   reduction = 'pca_umap',
                                   assay = 'SCRAN',
                                   clust.method = 'pca_seurat_louvain',
                                   column = 'neighbourhood_graph_res.0.6', pt.size = 0.6) +
  ggplot2::theme(legend.position = 'none')

plot.list[[2]]

plot.list[[3]] <- plot.reduced.dim(object = sample,
                                   reduction = 'pca_umap',
                                   assay = 'SCANPY',
                                   clust.method = 'pca_seurat_louvain',
                                   column = 'neighbourhood_graph_res.0.6', pt.size = 0.6) +
  ggplot2::theme(legend.position = 'none')

plot.list[[3]]

#Run featureplots
plot.features(object = sample, assay = 'SCT', slot = 'normalised',
              reduction = 'pca_umap', features = c('CD79A', 'BIRC5', 'TOP2A'), pt_size = 0.6)

#Run Trajectory analysis
traject_SC <- perform.slingshot.trajectory(object = sample,
                                      reduction = 'pca_umap',
                                      assay = 'SCT',
                                      clust.method = 'pca_seurat_louvain',
                                      column = 'neighbourhood_graph_res.0.6')

#Plot trajectory
plot.slingshot(result = traject_SC,
               relevant = F,
               Pseudotime = T)
