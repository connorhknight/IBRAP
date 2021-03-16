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

panc <- filter_IBRAP(object = panc, RAW_total.features < max.features & RAW_percent.mt < 8)

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

Once your normalisation has been performed, the next step is to reduce our dataset dimension to a more manageable quantitiy, here we have adopted either PCA or diffusion-based Manifold Approximation and Projection (dbMAP). These methods differ slightly, where PCA requires our data that has been scaled, whereas dbMAP functions well with either counts or normalised counts, but works poorly with scaled data. 

Our current repertoire of tools are outlined in the following table:

#### IBRAP repertoire

![tool_table](/figures/IBRAP_table.png)

By constructing these leading methods into an ease-of-use interchangeable pipeline, we enable flexibility and uniquity for scRNA-seq analyses. 

#### IBRAP pipeline

![pipeline](/figures/IBRAP_schemata.png)

#### IBRAP Tutorial

