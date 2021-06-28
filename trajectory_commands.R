traj_data <- readRDS('/Users/knight05/Downloads/pancreatic-alpha-cell-maturation_zhang.rds')

meta <- as.data.frame(traj_data$grouping)

traj <- createIBRAPobject(counts = example_dataset, 
                          original.project = '', 
                          meta.data = meta, 
                          min.cells = 0, 
                          min.features = 0)

colnames(traj@sample_metadata) <- c("original.project","RAW_total.counts","RAW_total.features","celltype")

traj <- perform.sct.normalisation(object = traj, 
                                  assay = 'RAW', 
                                  slot = 'counts')

traj <- perform.scran.normalisation(object = traj, 
                                    assay = 'RAW', 
                                    slot = 'counts', 
                                    vars.to.regress = 'RAW_total.counts')

traj <- perform.scanpy.normalisation(object = traj, 
                                     vars.to.regress = 'RAW_total.counts')

traj <- perform.pca(object = traj, 
                    assay = c('SCT', 'SCRAN', 'SCANPY'), 
                    n.pcs = 1:50, reduction.save = 'pca')

traj <- perform.dbmap(object = traj, 
                      assay = c('SCT', 'SCRAN', 'SCANPY'),  
                      reduction.save = 'dbmap')

traj <- perform.umap(object = traj, 
                     assay = c('SCT', 'SCRAN', 'SCANPY'), 
                     reduction = c('pca', 'dbmap'), 
                     reduction.save = c('pca_umap', 'dbmap_umap'), 
                     n_components = 3, 
                     n.dims = list(1:20, NULL))

traj <- perform.seurat.cluster(object = traj, assay = c('SCT', 'SCRAN', 'SCANPY'), reduction = c('pca', 'dbmap'), dims = list(1:20, NULL), assignment.df.name = c('pca_louvain', 'dbmap_louvain'))

traj_red <- perform.slingshot.trajectory(object = traj, reduction = 'pca_umap', 
                                         assay = 'SCT', clust.method = 'metadata', 
                                         column = 'celltype')

alpha_lineages <- plot.slingshot(result = traj_red, object = traj, assay = 'SCT', slot = 'normalised', feature = 'Blvra', relevant = F, Pseudotime = T)

alpha_res <- list(object = traj, slingshot = traj_red, plots = alpha_lineages)
