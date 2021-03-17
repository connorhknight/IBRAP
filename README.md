# IBRAP
 **I**ntegrated **B**enchmarking single-cell **R**NA-sequencing **A**utomated **P**ipeline endeavours to make optimal bioinformatic pipelines preconstructed into an ease-of-use tool, alongside providing gold-standard metrics to assess their effectiveness. 

## Tutorial

### Pre-processing for droplet-based technology:

Since droplet-based technology is becoming the majorly used technology in scRNA-seq we have adopted 2 novel methods that aim to rectify their specific problems. WARNING: Do not use these methods if you are not using a drople-based method.  

Method one:

Droplets are designed to capture a singular cell. However, infrequently (a small number in a sample) 2 or more cells are captured by a single droplet and thus, do not represet a true cell; therefore requiring their omission. We have incorporated scrublet, a python-based module that identifies multiplets/doublets by simualting multiplets through combining singlets from our sample and finding similar profiles in our observations.

```
pancreas.data <- perform.scrublet(counts = pancreas.data, 
                                  expected_doublet_rate = 0.025)
```

A simple graph is displayed when scrublet is performed - left corresponds with observed doublet scores and right displays simulated scores. A large number of highly scored doublets in the observed graph indicates a higher number of doublets in the sample.

![scrublet_fig](/figures/doublet_scores.png)

Droplets that reach the doublet score threshold are automatically omitted. 

Method two:

Another condundrum in scRNA-seq droplet-based platform is driven by ambient RNA. Thought to derive from cell stress or apoptosis, ambient RNA can cross-contaminate droplets, become barcoded and are amplified along with the native RNA; this inherently influences downstream analyses negatively. DecontX, a component of the celda R package is a leading technique for ambient RNA removal. 

```
pancreas.data <- perform.decontX(counts = pancreas.data)
```
A UMAP projection of the data is produced to display contamination scores across cells, this is displayed in a graph during decontamination. 

![decontX_fig](/figures/decontx_umap.png)

### IBRAP object creation: 

Now we have reduced the technical effects of droplet-based platforms we can proceed with producing our IBRAP class object. Remember, if you are not using droplet-based technology this should be the first component of your pipeline. 

```
panc <- createIBRAPobject(counts = pancreas.data,
                          original.project = 'pancreas_celseq', 
                          method.name = 'RAW',
                          min.cells = 3,
                          min.features = 200)
```
Initial filtering of the count matrix is important, cells with little information are likely incomplete, cells that are not spread across more than a few cells will also not add value to the analysis. Therefore, we remove them upon initiation. 

### Quality Control and further filtering: 

A common conudnrum shared across all platforms of scRNA-seq are mitochondrial RNA counts. During cell preparation, cell lysis can cause the mitochondria to leak mitochondrial RNA, which is not a true representation of the cells transcriptome and can influence downstream analyses. Thus, it is important to quanitify the amount of mitochondrial RNA in a cell. This can also be applied to ribosomal RNA. 

```
panc <- find_percentage_genes(object = panc, pattern = '^MT-', 
                              assay = 'RAW', slot = 'counts',
                              column.name = 'RAW_percent.mt')
panc <- find_percentage_genes(object = panc, pattern = 'RPL', 
                              assay = 'RAW', slot = 'counts',
                              column.name = 'RAW_percent.rp')
```
We can now visualise important metadata such as total counts per cell, total features per cell and metadata that we just generated: percentage mitochondrial RNA.

Violin plots enable us to see proportion distribution:

![vln_qc](/figures/QC.png)

Scatter plots allow us to visualise any correlations between the metadata: 

![scatter_qc](/figures/scatter_QC.png)

Finally, we remove any troublesome cells from our data object. 

```
sd.value <- sd(panc$RAW_total.features)
med.value <- median(panc$RAW_total.features)
max.features <- (sd.value*3)+med.value
```
The above calculation we have provided is commonly used to determine maximum number of features to retain, in this example we multiple SD by 3, however this can range between 2 for more stringent and 3 for more lenient. 

```
panc <- filter_IBRAP(object = panc, RAW_total.features < max.features & RAW_percent.mt < 8)
```

percentage mitochondrial is biologically specific, i.e cells in the lungs contain a higher percentage of mitochondrial counts than other tissues. 

Another common effect witnessed in scRNA-seq is the cell cycle effect where clustering captures the cell cycle which in most cases, is not desirable. Therefore, we can quanitify this effect to be regressed out later on. However, this is not always required in an analysis. 

```
panc <- add.cell.cycle(object = panc, assay = 'RAW', slot = 'counts', transform = TRUE)
```

Now that we have performed QC and filtration we can move on to pre-processing. 

### Pre-processing: 

Uneven sequencing depth is typically witnessed between cells inter-/intra-sample. Therefore, we have adopted 4 techniques to account for this: Scanpy (CPM), Scran, SCTransform for 3' and 5' end-bias transcripts, and TPM for full-length transcripts. The difference between the end-bias and full-length transcript methods are distinctly that TPM accounts for the length of the genes whereas the others are focussed purely on solving inequeal depth. 
Next, to reduce computational time and improve the accuracy of downstream analyses, we must identify highly variable genes within our sample. Each normalisation has it's own discovery method except TPM, which we have adopted Seurats method.
Finally, as a pre-requisite for Principal Component Analysis (PCA), we must scale and centre our data. Alongside this, we enable the regression of confounding factors in our data that we previously quanitifed (i.e. cell cycle score, mitochondrial RNA percentage, and more.). WARNING: this regression step is not designed to tackle batch effects, this will be dealt with at a different point. 

You may either use a single normalisation technique or compare multiple, since each method will be stored in a different method assay within our IBRAP object.

```
panc <- perform.scanpy.normalisation(object = panc, vars.to.regress = 'RAW_total.counts')

panc <- perform.scran.normalisation(object = panc, assay = 'RAW', slot = 'counts', 
                                    vars.to.regress = 'RAW_total.counts')

panc <- perform.sct.normalisation(object = panc, assay = 'RAW', slot = 'counts')

panc <- perform.tpm.normalisation(object = panc, assay = 'RAW', slot = 'counts', 
                                  vars.to.regress = 'RAW_total.counts')
```
Each normalisation technique is stored in an individual assay within our IBRAP object, these can be seen here:

```
panc@methods
```
henceforth, downstream functions will be applied to our normalisation assays, when we supply our assay as c('SCT', 'SCRAN', 'SCANPY', 'TPM'), the functions will iterate through these assays and perform their algorithms. 

Once your normalisation has been performed, the next step is to reduce our dataset dimension to a more manageable quantitiy, here we have adopted either PCA or diffusion-based Manifold Approximation and Projection (dbMAP). These methods differ slightly, where PCA requires data that has been scaled, whereas dbMAP functions well with either counts or normalised counts, but works poorly with scaled data. 

```
panc <- perform.pca(object = panc, 
                    assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                    n.pcs = 1:50, 
                    reduction.save = 'pca')
                    
panc <- perform.dbmap(object = panc, 
                      assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'),  
                      reduction.save = 'dbmap')
```

Now that we have reduced our datasets to a more manageable dimension, we can proceed to our next stage of the analysis, clustering & visualisation.

### Clustering, Visualisation & Benchmarking: 

Firstly, lets produce our visualisation reductions, these are required in single-cell since the count matrices are spare (contain many zeros). We provide 3 visualisations: t-SNE, UMAP, and lvish. 

```
panc <- perform.umap(object = panc, 
                     assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                     reduction = 'pca', 
                     reduction.save = 'pca_umap', 
                     n.dim = 1:42, 
                     n_components = 3)
panc <- perform.umap(object = panc, 
                     assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                     reduction = 'dbmap', 
                     reduction.save = 'dbmap_umap', 
                     n_components = 3)
                     
panc <- perform.tsne(object = panc, 
                     assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                     reduction = 'pca', 
                     reduction.save = 'pca_tsne', 
                     n.dim = 1:42, n_components = 3)
panc <- perform.tsne(object = panc, 
                     assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), 
                     reduction = 'dbmap', 
                     reduction.save = 'dbmap_tsne', n_components = 3)
                     
panc <- perform.lvish(object = panc, 
                      assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                      reduction = 'pca', 
                      reduction.save = 'pca_lvish', 
                      n.dim = 1:42, n_components = 3)
panc <- perform.lvish(object = panc, 
                      assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                      reduction = 'dbmap', 
                      reduction.save = 'dbmap_lvish', n_components = 3)
```


Now, we can attempt to discover our cell populations, we provide 3 techniques that are known to function the best: seurat graph-based clustering (fast), kmeans clustering performed on t-SNE reduction (fast), and SC3 which iterates many kmeans clustering and perform hierachial clustering on the results (slow). 

Seurat requires our computational reductions as input, we have both PCA and dbMAP, thus it is appropriate to calculate clusters on both:
```
panc <- perform.seurat.cluster(object = panc, assay = c('SCT', 'SCRAN', 'TPM', 'SCANPY'), reduction = c('pca', 'dbmap'), 
                               assignment.df.name = c('pca_seurat', 'dbmap_seurat'), dims = 1:42)                                                           
```
kmeans/t-SNE requires our t-SNE projections, again we calculated t-SNE from our two computational reductions and will use both:
```
panc <- perform.reduction.kmeans(object = panc, assay = c('SCT','SCRAN','SCANPY','TPM'), 
                         reduction = c('dbmap_tsne', 'pca_tsne'), dims = list(NULL, NULL), k = 10:16, 
                         assignment.df.name = c('dbmap_tsne_kmeans', 'pca_tsne_kmeans'), method = 'kmeans')
```
Finally, SC3 is performed on whol count matrices or normalised matrices, however scaled is not recommended. 
```
panc <- perform.sc3.slot.cluster(object = panc, 
                                 assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                                 slot = 'normalised', 
                                 HVGs = TRUE, 
                                 assignment.df.name = 'counts_SC3', 
                                 ks = 10:16, n.core = 3)
```
Now we have calculated our clusters, we must assess how effective they were with our benchmarking function:

When the ground truth cell types are available we enable the use of 5 benchmarking metrics: Average Silhouette Width (ASW), Dunn Index (DI), Connectivity, Adjusted Rand Index (ARI), and Normalized Mutual Information (NMI). Without the ground truth we are only able to utilise ASW, DI, and connectivity. Essentially, we want to witness an elevated score in ASW, DI, ARI, and NMI; in contrast, we aim to see a depleted score for connectivity.
```
panc <- benchmark.clustering(object = panc, 
                     assay = c('SCT', 'SCRAN', 'SCANPY', 'TPM'), 
                     clustering = c('pca_seurat', 'dbmap_seurat', 'dbmap_tsne_kmeans', 'pca_tsne_kmeans', 'counts_SC3'), 
                     reduction = c('pca_umap', 'dbmap_umap', 'dbmap_umap', 'pca_umap', 'pca_umap'), 
                     components = 1:3, 
                     dist.method = 'euclidean', 
                     ground.truth = panc$celltype)
```

Congratulations, you have just run a series of successful IBRAP-based pipelines. Now, with this large amount of data we must ascertain which combinations worked best for our dataset, we have constructed and interactive application to aid users in this discovery. However first, we must save our data object in rds format: 

```
saveRDS(object = panc, file = '/path/to/folder/location/pancreas_data.rd', compress = TRUE)
```

### IBRAP application:

Firstly, we must upload our RDS file which can be browsed for under 'RDS file upload', once the upload is completed we must load the file into our application using the 'Activate' button. Now, we can select from our assays under the 'Select assay' selection menu, this can be changed at any time. Finally, we need to select a visualisation projection to view out data in, this can be selected in the dropdown menu under 'Select reduction technique'. Now, we are ready to visualise some data. 

#### TAB: clustering

We can find our cell assignment for our assay under 'Select cell assignment', for simplicity we have also added our metadata to this section incase we need to visualise this information. Under 'Generate a dimensionality column' we have 3 dropdown menus: 'Select cell assignment column' where we navigate between the columns of our selected cell assignment, '2D or 3D?' allow us to choose between visualising the reduction in 2 or 3 dimensions, and 'point size' which allows us to adjust the size of our cells. A legend is also displayed on the right-hand side of our plot which displays which colour corresponds to which cell assignment. 

![step1](/figures/step1.png)

To help us to understand which parameter functioned optimally, benchmarking metrics for the clustering is displayed underneath. These metrics can be used as a guide but may not always indicate the best result or the correct number of clusters. 

![step2](/figures/step2.png)

#### TAB: features

A fundamental part of scRNA-seq analyses is understanding the biology behind the results. Dividing our cell populations up according to their transcriptomic profiles is excellent, but we must understand which assignment belongs to which cell type. Therefore, you can visualise gene expression in our feature tab. 

In the first section of this tab we can produce multiple feature plots displaying the gene expression of specific markers. This enables us to identify canonical markers that aid in their biological identification. 

![step3](/figures/step3.png)

By constructing these leading methods into an ease-of-use interchangeable pipeline, we enable flexibility and uniquity for scRNA-seq analyses. As you can see, the feature we have selected are obvious exclusively expressed in high quantities in a particular demographic. We can begin to understand which cluster belongs to which cell type.

![step4](/figures/step4.png)

To help us get a more categorical understanding of expression we can utilise violin plots split by cluster assignments, this in combination with feature plots is a powerful tool for visualising gene expression. 

![step5](/figures/step5.png)

Finally, in the case where we wish to see the distribution of expression amongst a large number of cells we can use a heatmap. For easy interpretation of the values the expression matrix z-score has been calculated. A z-score is simple to interpret, a negative value means the expression is below the average of the whole population, whereas a positive value indicates a higher than average expression. 

With this visualisation application, we can begin to dissect our best pipeline combination and visually interpret the biology of the results. 

Our current repertoire of tools are outlined in the following table:

#### IBRAP repertoire

![tool_table](/figures/IBRAP_table.png)

#### IBRAP pipeline

![pipeline](/figures/IBRAP_schemata.png)

