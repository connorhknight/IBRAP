# IBRAP
 **I**ntegrated **B**enchmarking single-cell **R**NA-sequencing **A**utomated **P**ipeline endeavours to make optimal bioinformatic pipelines preconstructed into an ease-of-use tool, alongside providing gold-standard metrics to assess their effectiveness. 

## Tutorial

Quality Control & Pre-processing for droplet-based technology:

Since droplet-based technology is becoming the majorly used technology in scRNA-seq we have adopted 2 novel methods that aim to rectify their specific problems. WARNING: Do not use these methods if you are not using a drople-based method.  

Method one:

Droplets are designed to capture a singular cell. However, infrequently (a small number in a sample) 2 cells are captured by a single droplet and thus, do not represet a true cell; therefore requiring their omission. We have incorporated scrublet, a python-based module that identifies doublets by simualting their own doublets for our dataset and finding similar profiles in our observations.

```
panc <- perform.scrublet(object = panc, 
                         assay = 'RAW', 
                         slot = 'counts',
                         split.by = 'sample', 
                         expected_doublet_rate = 0.025)

```

Our current repertoire of tools are outlined in the following table:

#### IBRAP repertoire

![tool_table](/figures/IBRAP_table.png)

By constructing these leading methods into an ease-of-use interchangeable pipeline, we enable flexibility and uniquity for scRNA-seq analyses. 

#### IBRAP pipeline

![pipeline](/figures/IBRAP_schemata.png)

#### IBRAP Tutorial

