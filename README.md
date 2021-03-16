# IBRAP
 **I**ntegrated **B**enchmarking single-cell **R**NA-sequencing **A**utomated **P**ipeline endeavours to make optimal bioinformatic pipelines preconstructed into an ease-of-use tool, alongside providing gold-standard metrics to assess their effectiveness. 

## Tutorial

Quality Control & Pre-processing for droplet-based technology:

Since droplet-based technology is becoming the majorly used technology in scRNA-seq we have adopted 2 novel methods that aim to rectify their specific problems. WARNING: Do not use these methods if you are not using a drople-based method.  

Method one:

Droplets are designed to capture a singular cell. However, infrequently (a small number in a sample) 2 or more cells are captured by a single droplet and thus, do not represet a true cell; therefore requiring their omission. We have incorporated scrublet, a python-based module that identifies multiplets/doublets by simualting multiplets through combining singlets from our sample and finding similar profiles in our observations.

```
panc <- perform.scrublet(object = panc, 
                         assay = 'RAW', 
                         slot = 'counts',
                         split.by = 'sample', 
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

Our current repertoire of tools are outlined in the following table:

#### IBRAP repertoire

![tool_table](/figures/IBRAP_table.png)

By constructing these leading methods into an ease-of-use interchangeable pipeline, we enable flexibility and uniquity for scRNA-seq analyses. 

#### IBRAP pipeline

![pipeline](/figures/IBRAP_schemata.png)

#### IBRAP Tutorial

