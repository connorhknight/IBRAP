Koh <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/Koh.SCE.rds')
Kumar.obj <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/Kumar.SCE.rds')

Zhengmix4eq <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/Zhengmix4eq.SCE.rds')
Zhengmix4uneq <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Zhengmix4uneq.SCE.rds')
Zhengmix8eq <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Zhengmix8eq.SCE.rds')

mixology10x3cl <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/mixology10x3cl.SCE.rds')
mixology10x5cl <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/mixology10x5cl.SCE.rds')

simMix1 <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/simMix1.SCE.rds')
simMix2 <- readRDS('/Users/knight05/Results/scRNA-seq/Benchmarking_datasets/simMix2.SCE.rds')

metadata <- readRDS(file = "/Users/knight05/Downloads/pancreas_v3_files/pancreas_metadata.rds")

reticulate::use_python('/Users/knight05/Library/r-miniconda/envs/r-reticulate/bin/python')
reticulate::import('scrublet', convert = FALSE)

# Read 10x

library(utils)
library(Matrix)

Read10X_output <- function(directory) {
  dir.files <- list.files(path = directory)
  if(!'matrix.mtx' %in% dir.files) {
    cat(crayon::cyan('Expected file: matrix.mtx\n'))
  }
  if(!'genes.tsv' %in% dir.files) {
    cat(crayon::cyan('Expected file: genes.tsv\n'))
  }
  if(!'barcodes.tsv' %in% dir.files) {
    cat(crayon::cyan('Expected file: barcodes.tsv\n'))
  }
  count.file <- paste0(directory, '/matrix.mtx')
  genes.file <- paste0(directory, '/genes.tsv')
  barcodes.files <- paste0(directory, '/barcodes.tsv')
  cat(crayon::cyan('Files read\n'))
  genes_ensembl <- read.table(file = genes.file, sep = '\t', header = FALSE)
  barcodes <- read.table(file = barcodes.files, sep = '\t', header = FALSE)
  mm <- Matrix::readMM(count.file)
  if(isUnique(genes_ensembl$V2) == FALSE) {
    cat(crayon::cyan('Non-unique features identified\n'))
    print(duplicated(genes_ensembl$V2))
    genes_ensembl$V2 <- make.unique(genes_ensembl$V2)
  }
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  cat(crayon::cyan('Success: Matrix concatenated\n'))
  mm <- as.matrix(mm)
  return(mm)
}

marrow_E <- Read10X_output(directory = '/Users/knight05/Raw_Data/Database_samples/healthy_references/BMMC_atlas/marrow_Ck')

# Metadata generator

library(SingleCellExperiment)

metadata.generator <- function(object, 
                               assay) {
  cat(crayon::cyan('Annotated with fundamental metadata for assay:', assay, '\n'))
  colData(object)[[paste0(assay, '_total.counts')]] <- colSums(assay(object, assay))
  colData(object)[[paste0(assay, '_total.features')]] <- colSums(assay(object, assay) != 0)
  cat(crayon::cyan('Completed\n'))
  return(object)
}

# SCE

create.sce.object <- function(counts, 
                              project.name, 
                              min.cells, 
                              metadata = NULL) {
  sce <- SingleCellExperiment(assays = list(counts = counts))
  cat(crayon::cyan('SCE object created\n'))
  colData(sce)$original.project <- project.name
  sce <- metadata.generator(object = sce, assay = 'counts')
  if(!is.null(metadata))  {
    colData(sce) <- cbind(colData(sce), metadata)
    cat(crayon::cyan('Custom metadata appended\n'))
  }
  sce <- sce[!duplicated(rownames(assay(sce, 'counts'))), ]
  cat(crayon::cyan('Completed\n'))
  return(sce)
}

# Load practice data

pancreas.data <- readRDS(file = "~/Raw_Data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS('~/Raw_Data/pancreas_v3_files/pancreas_metadata.rds')
pancreas.data <- as.matrix(pancreas.data)

panc <- create.sce.object(counts = pancreas.data, 
                              project.name = 'panc', 
                              min.cells = 3, metadata = metadata)

panc$celltype <- Koh$phenoid

#reticulate::import('scrublet')

# scrublet

perform.scrublet <- function(object, 
                             split.by, 
                             total_counts = NULL, 
                             sim_doublet_ratio = 2.0, 
                             n_neighbors = NULL, 
                             expected_doublet_rate = 0.075, 
                             stdev_doublet_rate = 0.02, 
                             random_state = 0L) {
  cat(crayon::cyan('Initialising scrublet\n'))
  scrublet <- reticulate::import('scrublet', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  if(isS4(object) == FALSE) {
    cat(crayon::cyan('Only an object of class S4 can be used\n'))
  } else {
    object
    raw_counts_list <- list()
    seperator <- unique(colData(object)[,split.by])
    for(l in seperator) {
      cat(crayon::cyan('###############################\n'))
      cat(crayon::cyan(paste0('scrublet analysing: ', l, '\n')))
      isolated <- object[,object[[split.by]]==l]
      isolated.2 <- assay(isolated, 'counts')
      raw_counts <- t(as.data.frame(as.matrix(isolated.2)))
      scrub1 <- scrublet$Scrublet(counts_matrix = reticulate::r_to_py(raw_counts))
      cat(crayon::cyan('scrublet object created\n'))
      res1 <- scrub1$scrub_doublets(min_counts = 1, 
                                    min_cells = 1, 
                                    min_gene_variability_pctl = 85, 
                                    verbose = TRUE)
      cat(crayon::cyan('doublets detected\n'))
      raw_counts <- t(as.data.frame(raw_counts))
      raw_counts <- raw_counts[,!reticulate::py_to_r(res1)[[2]] == TRUE]
      cat(crayon::cyan('matrix scrubbed\n'))
      raw_counts_list[[l]] <- raw_counts
    }
    raw_counts <- do.call('cbind', raw_counts_list)
    object <- object[,colnames(raw_counts)]
    assay(object, 'counts') <- raw_counts
    return(object)
    rm(obj, raw_counts, scrub1, res1, scrubbed)
  }
}

pancreas <- perform.scrublet(object = pancreas, 
                             split.by = 'tech', 
                             expected_doublet_rate = 0.025)

# decontamination

perform.decontX <- function(object,
                            z = NULL,
                            batch = NULL,
                            maxIter = 500,
                            delta = c(10, 10),
                            estimateDelta = TRUE,
                            convergence = 0.001,
                            iterLogLik = 10,
                            varGenes = 5000,
                            dbscanEps = 1,
                            seed = 12345,
                            logfile = NULL,
                            verbose = TRUE) {
  d <- celda::decontX(x = object,
                      z = z,
                      batch = colData(object)[,batch],
                      maxIter = maxIter,
                      delta = delta,
                      estimateDelta = estimateDelta,
                      convergence = convergence,
                      iterLogLik = iterLogLik,
                      varGenes = varGenes,
                      dbscanEps = dbscanEps,
                      seed = seed,
                      logfile = logfile,
                      verbose = verbose)
  cat(crayon::cyan('Decontamination comlpleted\n'))
  
  cat(crayon::cyan(paste0(formatC(sum(d$decontX_contamination)/length(d$decontX_contamination), 
                                  digits = 2), 
                          '% average contamination\n')))
  
  clean.matrix <- assay(d, 'decontXcounts')
  cat(crayon::cyan('Matrix isolated\n'))
  clean.matrix <- round(clean.matrix)
  zero.samples <- colSums(clean.matrix) > 0
  object <- object[,zero.samples]
  clean.matrix <- clean.matrix[,zero.samples]
  cat(crayon::cyan('converted to integer\n'))
  assay(object, 'decontXcounts') <- clean.matrix
  cat(crayon::cyan('Added matrix\n'))
  
  object <- metadata.generator(object = object, 
                               assay = 'decontXcounts')
  
  cat(crayon::cyan('Finished\n'))
  return(object)
}

pancreas <- perform.decontX(object = pancreas, 
                            batch = 'tech')

post_decontX <- pancreas

find_percentage_genes <- function(object, 
                                  pattern='^MT-', 
                                  which.assay='decontXcounts', 
                                  suffix = 'percent.mt') {
  cat(crayon::cyan('Calculating percentage\n'))
  temp <- colSums(assay(object, which.assay)[grep(pattern = pattern, x = rownames(assay(object, which.assay))),]) / colSums(assay(object, which.assay)) * 100
  cat(crayon::cyan('Percentage calculated\n'))
  temp <- as.data.frame(temp)
  colnames(temp) <- paste0(which.assay, '_', suffix)
  if(paste0(which.assay, '_', suffix) %in% colnames(colData(object))) {
    cat(crayon::cyan('Removing old metadata column\n'))
    colData(object) <- colData(object)[,colnames(colData(object)) != paste0(which.assay, '_', suffix)]
  }
  colData(object)[[paste0(which.assay, '_', suffix)]] <- temp
  cat(crayon::cyan('Appending new column\n'))
  return(object)
}

pancreas <- find_percentage_genes(object = pancreas, pattern = '^MT-', which.assay = 'counts', suffix = 'percent.mt')

filter_features_cells <- function(object, 
                                  use.assay, 
                                  max.features=NULL, 
                                  min.features=NULL, 
                                  max.mt=NULL, 
                                  min.cell=NULL, 
                                  omit.zero.features=FALSE) {
  if(omit.zero.features != NULL) {
    rowData(object)$total.cells <- rowSums(assay(object, use.assay) != 0)
    cat(crayon::cyan('Empty features removed\n'))
  } else if (omit.zero.features == NULL) {
    cat(crayon::cyan('Empty features not removed\n'))
  }
  if (max.features != NULL) {
    object <- object[rownames(subset(rowData(object), total.cells > as.numeric(min.cell))),]
    cat(crayon::cyan('Minimum cells removed\n'))
  } else if (omit.zero.features == NULL) {
    cat(crayon::cyan('Minimum cells not removed\n'))
  }
  if (min.features != NULL) {
    object <- object[,rownames(subset(colData(object), total.features < as.numeric(max.features)))]
    cat(crayon::cyan('Maximum features removed\n'))
  } else if (omit.zero.features == NULL) {
    cat(crayon::cyan('Minimum cells not removed\n'))
  }
  if (max.mt != NULL) {
    object <- object[,rownames(subset(colData(object), percent.mt < as.numeric(max.mt)))]
    cat(crayon::cyan('Maximum mitochondrial feature percentage removed\n'))
  } else if (max.mt == NULL) {
    cat(crayon::cyan('Maximum mitochondrial feature percentage not removed\n'))
  }
  cat(crayon::cyan('Filtration complete!\n'))
  return(object)
}

sd.value <- sd(pancreas$counts_total.features)
med.value <- median(pancreas$counts_total.features)
max.features <- (sd.value*3)+med.value

# pancreas <- filter_features_cells(object = pancreas, 
#                                   use.assay = 'decontXcounts', 
#                                   max.features = max.features, 
#                                   min.features = 200, 
#                                   max.mt = 10, 
#                                   min.cell = 3)



pancreas <- suppressWarnings(pancreas[,colData(pancreas)[['counts_total.features']] > 200 & < max.features & colData(pancreas)[['counts_percent.mt']] < 10])
pancreas <- suppressWarnings(pancreas[,colData(pancreas)[['counts_total.features']] ])
pancreas <- suppressWarnings(pancreas[,colData(pancreas)[['counts_percent.mt']] < 10])
pancreas <- pancreas[rowSums(assay(pancreas, 'counts') != 0) >= 3, ]

add.cell.cycle <- function(object, 
                           assay, 
                           transform, ...) {
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/data/Homo_sapiens.csv', header = TRUE, sep = ',')
  cat(crayon::cyan('Cell cycle genes loaded\n'))
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = assay(object, assay))
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Cell cycle scores identified\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- as.Seurat(object, counts = NULL, data = assay)
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::CellCycleScoring(object = seuobj, s.features = r[55:97,3], g2m.features = r[1:54,3], ...)
    cat(crayon::cyan('Data transformed\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, sum(length(colnames(seuobj@meta.data))-2):length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}

panc <- add.cell.cycle(object = panc, assay = 'counts', transform = TRUE)

add.feature.score <- function(object, 
                              assay, 
                              transform, 
                              features, 
                              column.name,
                              ...) {
  genes <- rownames(object)
  genes <- list(genes[genes %in% features])
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = assay(object, assay))
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = genes, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    colData(object)[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- as.Seurat(object, counts = NULL, data = assay)
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    colData(object)[[column.name]] <- seuobj@meta.data[, length(colnames(seuobj@meta.data))]
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}

pancreas <- add.feature.score(object = pancreas, 
                              assay = 'counts', 
                              transform = TRUE, 
                              features = c('BAG3', 'BLOC1S5-TXNDC5', 'CALU', 'DNAJB1', 'DUSP1', 'EGR1', 
                                           'FOS', 'FOSB', 'HIF1A', 'HSP90AA1', 'HSP90AB1', 'HSP90AB2P', 
                                           'HSP90AB3P', 'HSP90B1', 'HSPA1A', 'HSPA1B', 'HSPA6', 'HSPB1', 
                                           'HSPH1', 'IER2', 'JUN', 'JUNB', 'NFKBIA', 'NFKBIZ', 'RGS2', 
                                           'SLC2A3', 'SOCS3', 'UBC', 'ZFAND2A', 'ZFP36', 'ZFP36L1'), 
                              column.name = 'StressScore')

library(Seurat)
library(sctransform)

perform.sct.normalisation <- function(object, 
                                      which.assay, 
                                      prefix = 'sctransform',
                                      do.scale = TRUE,
                                      do.center = TRUE,
                                      min_cells = 3,
                                      ...) {
  cat(crayon::cyan('Converting to Seurat object\n'))
  seuratobj <- Seurat::CreateSeuratObject(counts = assay(object, which.assay), project = 'NA')
  cat(crayon::cyan('Initiating SCTransform\n'))
  seuratobj <- Seurat::SCTransform(object = seuratobj, do.scale = do.scale, do.center = do.center, min_cells = min_cells, ...)
  cat(crayon::cyan('SCTransform completed!\n'))
  genes <- intersect(rownames(object), rownames(as.matrix(seuratobj@assays$SCT@data)))
  object <- object[genes,]
  metadata(object)$HVGs <- seuratobj@assays$SCT@var.features
  assay(object, paste0(prefix, '_counts')) <- as.matrix(seuratobj@assays$SCT@counts)[genes,]
  assay(object, paste0(prefix, '_data')) <- as.matrix(seuratobj@assays$SCT@data)[genes,]
  reducedDim(object, paste0(prefix, '_scaled')) <- as.matrix(t(seuratobj@assays$SCT@scale.data))
  cat(crayon::cyan('Populated SCE\n'))
  return(object)
}

pancreas.tmp <- perform.sct.normalisation(object = pancreas, 
                                      which.assay = 'counts')

library(scran)
library(scater)

perform.scran.normalisation <- function(object, 
                                        new.assay = 'scran', 
                                        max.cluster.size = 1000, 
                                        scaling=NULL,
                                        do.log=TRUE, 
                                        center_size_factors=TRUE,
                                        split.by) {
  
  t <-unique(colData(object)[[split.by]])
  list.matrix <- list()
  
  for(o in t) {
    cat(crayon::cyan(paste0('analysing sample: ', o, '\n')))
    sub <- object[,colData(object)[[split.by]] == o]
    cat(crayon::cyan('initialisaing quickCluster\n'))
    clusters <- quickCluster(sub)
    cat(crayon::cyan('initialisaing computeSumFactors\n'))
    sce.scran <- computeSumFactors(sub, clusters=clusters, max.cluster.size=max.cluster.size, scaling=scaling)
    cat(crayon::cyan('initialisaing LogNormCounts\n'))
    log <- logNormCounts(x = sce.scran, log = do.log, center_size_factors=center_size_factors, exprs_values='counts')
    cat(crayon::cyan('Adding matrix\n'))
    y <- object
    altExp(y) <- NULL
    list.matrix[[o]] <- logcounts(log)
  }
  cat(crayon::cyan('Binding matrices\n'))
  f <- do.call('cbind', list.matrix)
  f <- as.matrix(as.data.frame(f))
  new <- assay(object, 'counts') 
  new <- new[rownames(f),]
  cat(crayon::cyan('New matrix created\n'))
  for(x in colnames(new)){
    new[,x] <- f[,x]
  }
  cat(crayon::cyan('Matrix populated\n'))
  empty <- object
  assay(empty, new.assay) <- new
  cat(crayon::cyan('Added to object\n'))
  return(empty)
}




perform.tpm.normalisation <- function(object, 
                                      which.assay, 
                                      new.assay = 'tpm', 
                                      log.transform) {
  
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/mart_export.csv', header = TRUE, sep = ',')
  r$Gene.length <- r$Gene.end..bp. - r$Gene.start..bp.
  
  subset <- r[r$Gene.name %in% rownames(object),]
  
  cat(crayon::cyan('Matrix subsetted\n'))
  
  rownames(subset) <- make.unique(names = as.character(subset$Gene.name), '.')
  
  cat(crayon::cyan('Rownames added\n'))
  
  meta <- rowData(object)[intersect(rownames(rowData(object)), rownames(subset)),]
  
  cat(crayon::cyan('Gene names interesected\n'))
  
  object <- object[intersect(rownames(rowData(object)), rownames(subset)),]
  
  ordered <- subset[match(rownames(rowData(object)), rownames(subset)),]
  
  cat(crayon::cyan('Matrices ordered\n'))
  
  rowData(object)$length <- ordered$Gene.length
  
  mat <- assay(object, which.assay)
  
  cat(crayon::cyan('Calculated counts/feature length\n'))
  
  calc <- sweep(mat, 1, as.numeric(rowData(object)$length), `/`)
  
  scale.factor <- colSums(calc)/1000000
  
  calc2 <- sweep(calc, 2, as.numeric(scale.factor), `/`)
  
  cat(crayon::cyan('Calculations completed\n'))
  
  if(log.transform == TRUE) {
    cat(crayon::cyan('log(x+1) transforming\n'))
    mat <- log(calc2+1)
    cat(crayon::cyan('Transformation completed\n'))
  }
  
  assay(object, new.assay) <- mat
  
  cat(crayon::cyan('Completed!\n'))
  
  return(object)
}

processed.scrublet.decontx.mt.filter.tpm <- perform.tpm(object = processed.scrublet.decontx.mt.filter, 
                                                        which.assay = 'decontXcounts', 
                                                        log.transform = TRUE)

# scobj <- sc$AnnData(X = t(assay(pancreas.test, 'decontXcounts')), obs = as.data.frame(colData(pancreas.test)))
# scobj$obs_names <- colnames(pancreas.test)
# scobj$var_names <- rownames(pancreas.test)
# sc$pp$normalize_total(scobj, target_sum = 1e6)
# sc$pp$log1p(scobj)
# sc$pp$highly_variable_genes(adata = scobj, min_mean=0.0125, max_mean=6, min_disp=0.5)
# varfeat <- rownames(scobj$var[scobj$var$highly_variable == TRUE,])
# scobj <- sc$AnnData(X = t(assay(pancreas.test, 'decontXcounts')[varfeat,]), obs = as.data.frame(colData(pancreas.test)))
# sc$pp$regress_out(adata = scobj, keys = c('decontXcounts_total.counts', 'decontXcounts_percent.mt'))
# sc$pp$scale(scobj)
# sc$tl$pca(scobj)
# sc$pp$neighbors(scobj)
# sc$tl$umap(scobj)
# sc$pl$umap(scobj, color = c('tech'))

perform.scanpy.normalisation <- function(object, assay, new.assay='scanpy', 
                                         target_sum = 1e6, exclude_highly_expressed = FALSE,  
                                         max_fraction = 0.05, key_added = 'scanpy_norm_factor') {
  sc <- reticulate::import('scanpy')
  scobj <- sc$AnnData(X = t(assay(object, assay)))
  scobj$obs_names <- as.factor(colnames(object))
  scobj$var_names <- as.factor(rownames(object))
  
  if(!is.null(target_sum) & !is.null(key_added)) {
    sc$pp$normalize_total(adata = scobj, target_sum = as.integer(target_sum), 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction), key_added = as.character(key_added))
  } else if (!is.null(target_sum)) {
    sc$pp$normalize_total(adata = scobj, target_sum = as.integer(target_sum), 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction))
  } else if (!is.null(key_added)) {
    sc$pp$normalize_total(adata = scobj, 
                          exclude_highly_expressed = as.logical(exclude_highly_expressed), 
                          max_fraction = as.integer(max_fraction))
  } else {
    sc$pp$normalize_total(adata = scobj)
  }
  
  sc$pp$log1p(scobj)
  
  mat <- t(scobj$X)
  rownames(mat) <- rownames(object)
  colnames(mat) <- colnames(object)
  
  assay(object, new.assay) <- mat
  return(object)
}

perform.scanpy.hvg <- function(object, assay, n_top_genes = 1000, 
                               max_mean = 6, min_mean = 0.0125, 
                               min_disp = 0.5, span = 0.3, 
                               n_bins = 20, flavor = 'seurat', batch_key = NULL) {
  
  sc <- reticulate::import('scanpy')
  scobj <- sc$AnnData(X = t(assay(object, assay)))
  scobj$obs_names <- as.factor(colnames(object))
  scobj$var_names <- as.factor(rownames(object))
  
  if (!is.null(n_top_genes) & !is.null(batch_key)) {
    
    sc$pp$highly_variable_genes(adata = scobj, 
                                n_top_genes = as.integer(n_top_genes), 
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor), 
                                batch_key = as.character(batch_key))
    
  } else if (!is.null(n_top_genes)) {
    
    sc$pp$highly_variable_genes(adata = scobj, 
                                n_top_genes = as.integer(n_top_genes), 
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor))
    
  } else if (!is.null(batch_key)) {
    
    sc$pp$highly_variable_genes(adata = scobj,
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor), 
                                batch_key = as.character(batch_key))
    
  } else {
    
    sc$pp$highly_variable_genes(adata = scobj, 
                                min_mean = as.integer(min_mean), 
                                max_mean = as.integer(max_mean), 
                                min_disp = as.integer(min_disp), 
                                span = as.integer(span),
                                n_bins = as.integer(n_bins), 
                                flavor = as.character(flavor))
    
  }
  
  vargenes <- rownames(scobj$var[scobj$var$highly_variable == TRUE,])
  cat(crayon::cyan(paste0('Discovered variable features: ', vargenes)))
  metadata(object)$HVGs <- vargenes
  
  return(object)
  
}

perform.scanpy.scale <- function(object, 
                                 assay,
                                 assay.name = 'scanpy_HVGs_scaled',
                                 var.to.regress, 
                                 n_jobs = NULL, 
                                 zero_center = TRUE, 
                                 max_value = NULL, 
                                 obsm = NULL) {
  
  sc <- reticulate::import('scanpy')
  scobj <- sc$AnnData(X = t(assay(object, assay))[,metadata(object)$HVGs])
  scobj$obs_names <- as.factor(colnames(object))
  scobj$var_names <- as.factor(metadata(object)$HVGs)
  
  if(as.data.frame(colData(object)) >= 1) {
    scobj$obs <- as.data.frame(colData(object))
  }
  
  if(!is.null(var.to.regress)) {
    if(!is.null(n_jobs)) {
      sc$pp$regress_out(adata = scobj, keys = var.to.regress, n_jobs = as.integer(n_jobs))
    } else {
      sc$pp$regress_out(adata = scobj, keys = var.to.regress)
    }
  }
  
  if(!is.null(max_value) & !is.null(obsm)) {
    sc$pp$scale(scobj, zero_center = as.logical(zero_center), max_value = as.integer(max_value), obsm = as.character(obsm))
  } else if(!is.null(max_value)) {
    sc$pp$scale(scobj, zero_center = as.logical(zero_center), max_value = as.integer(max_value))
  } else if(!is.null(obsm)) {
    sc$pp$scale(scobj, zero_center = as.logical(zero_center), obsm = as.character(obsm))
  } else {
    sc$pp$scale(scobj)
  }
  
  mat <- scobj$X
  
  colnames(mat) <- metadata(object)$HVGs
  rownames(mat) <- colnames(object)
  
  reducedDim(object, assay.name) <- mat
  
  return(object)
}

perform.bbknn <- function(object, 
                          reduction,
                          graph.name = 'bbknn',
                          column.correct,
                          n_pcs = NULL, 
                          trim = NULL, 
                          n_trees = 10,
                          use_faiss = TRUE,
                          set_op_mix_ratio = 1.0,
                          local_connectivity= 1) {
  
  sc <- reticulate::import('scanpy')
  print('.')
  scobj <- sc$AnnData(X = reducedDim(object, reduction))
  print('.')
  scobj$obs_names <- as.factor(colnames(object))
  print('.')
  scobj$var_names <- as.factor(colnames(reducedDim(object, reduction)))
  print('.')
  scobj$obsm$update(X_pca = reducedDim(object, reduction))
  print(scobj)
  
  if(as.data.frame(colData(object)) >= 1) {
    print('obs')
    print(as.data.frame(colData(object)))
    scobj$obs <- pd$DataFrame(data = as.data.frame(colData(object)))
  }
  
  if(is.null(n_pcs)) {
    print('npcs calculated')
    n_pcs <- as.integer(length(colnames(reducedDim(object, reduction))))
  }
  
  if(is.null(trim)) {
    print('initialising bbknn')
    sc$external$pp$bbknn(scobj, 
                         batch_key = as.character(column.correct), 
                         approx = as.logical(FALSE), 
                         n_pcs = n_pcs,
                         n_trees = as.integer(n_trees),
                         use_faiss = as.logical(use_faiss),
                         set_op_mix_ratio = set_op_mix_ratio,
                         local_connectivity = local_connectivity)
  } else if (!is.null(trim)) {
    sc$external$pp$bbknn(scobj, 
                         batch_key= as.character(column.correct), 
                         approx = as.logical(FALSE), 
                         n_pcs = n_pcs,
                         trim = as.integer(trim),
                         n_trees = as.integer(n_trees),
                         use_faiss = as.logical(use_faiss),
                         set_op_mix_ratio = set_op_mix_ratio,
                         local_connectivity = local_connectivity)
  }
  
  graph.list <- list()
  connectivities <- scobj$obsp[['connectivities']]
  colnames(connectivities) <- colnames(object)
  rownames(connectivities) <- colnames(object)
  distances <- scobj$obsp[['distances']]
  colnames(distances) <- colnames(object)
  rownames(distances) <- colnames(object)
  
  graph.list[['connectivities']] <- connectivities
  graph.list[['distances']] <- distances
  print(graph.list)
  metadata(object)[['graphs']][[graph.name]] <- graph.list
  
  return(object)
  
}

test <- perform.bbknn(object = pancreas.test, 
                      reduction = 'uncorrected_dbmap', 
                      graph.name = 'dbmap_bbknn', 
                      column.correct = 'tech')

perform.seurat.hvg <- function(object, 
                               assay, 
                               nfeatures=1500, 
                               feat.to.omit=NULL, 
                               ...) {
  if(typeof(object) != 'S4') {
    cat(crayon::cyan('Must be an S4 SCE object\n'))
    return(NULL)
  } else {
    tmp <- suppressWarnings(Seurat::as.Seurat(object, counts = NULL, data = assay))
    cat(crayon::cyan('SCE converted to seurat object\n'))
    f <- Seurat::FindVariableFeatures(object = tmp[!(rownames(tmp) %in% feat.to.omit)], nfeatures = nfeatures, ...)
    tmp <- f@assays$RNA@var.features[f@assays$RNA@var.features %in% feat.to.omit]
    cat(crayon::cyan('Variable features identified\n'))
    metadata(object)$HVGs <- tmp
    cat(crayon::cyan('HVGs added to object\n'))
  }
  return(object)
}

pancreas.test.1 <- perform.seurat.hvg(object = pancreas.test.1, assay = 'sctransform', nfeatures = 1500)

perform.scran.hvg <- function(object, 
                              assay.name, 
                              nfeatures = 1500, 
                              method = 'fish', 
                              unwanted.variance = NULL) {
  if(typeof(object) != 'S4') {
    cat(crayon::cyan('Must be an S4 SCE object\n'))
    return(NULL)
  } else {
    matrix <- assay(object, assay.name)
    cat(crayon::cyan('Matrix isolated\n'))
  }
  dec <- scran::modelGeneVar(x = matrix)
  top.hvgs <- scran::getTopHVGs(dec, n=nfeatures)
  cat(crayon::cyan('Variable features identified\n'))
  metadata(object)$HVGs <- top.hvgs
  cat(crayon::cyan('HVGs added to object\n'))
  return(object)
}

pancreas.scran <- perform.scran.hvg(object = pancreas, assay.name = 'sctransform')

perform.seurat.scale <- function(object, 
                                 use.assay, 
                                 unwanted.variance = NULL, 
                                 scale, 
                                 centre, 
                                 ...) {
  y <- object
  tmp <- suppressWarnings(Seurat::as.Seurat(y, counts = NULL, data = use.assay))
  cat(crayon::cyan('SCE converted to seurat object\n'))
  tmp <- Seurat::ScaleData(tmp, vars.to.regress = unwanted.variance, do.scale = scale, do.center = scale, ...)
  cat(crayon::cyan('Data scaled\n'))
  reducedDim(y, 'HVGs_scaled') <- t(as.matrix(tmp@assays$RNA@scale.data))
  cat(crayon::cyan('Scaled matrix attached to object\n'))
  return(y)
}

pancreas.test.1 <- perform.seurat.scale(object = pancreas.test.1, use.assay = 'sctransform', unwanted.variance = c('percent.mt'), scale = TRUE, centre = TRUE)

perform.pca <- function(object, 
                        reduced.dim = 'HVGs_scaled', 
                        reduction.save='pca', 
                        ...) {
  cat(crayon::cyan('Initialising PCA\n'))
  a <- PCAtools::pca(mat = t(reducedDim(object, reduced.dim)), center = FALSE, scale = FALSE, ...)
  reducedDim(object, reduction.save) <- as.matrix(a$rotated[,1:50])
  cat(crayon::cyan('PCA completed\n'))
  return(object)
}

pancreas <- perform.pca(object = pancreas, reduction.save)

plot.red.sd <- function(object, 
                        reduction, 
                        n.dim, 
                        cex.names = 0.6) {
  gg <- apply(X = reducedDim(object, as.character(reduction))[,n.dim], MARGIN = 2, FUN = sd)
  barplot(gg, cex.names = cex.names, las = 2, main = paste0(reduction, '_var'))
}

plot.red.sd(object = pancreas.scran, reduction = 'uncorrected_pca', n.dim = 1:25)

perform.dbmap <- function(object, 
                          assay, 
                          n_components = 100, 
                          n_neighbors = 15, 
                          reduction.save='dbmap') {
  scipy.sparse <- reticulate::import('scipy.sparse', convert = FALSE)
  print('.')
  dbmap <- reticulate::import('dbmap', convert = FALSE)
  print('.')
  cellnames <- colnames(assay)
  print('.')
  data <- scipy.sparse$csr_matrix(reticulate::r_to_py(t(assay)))
  print('.')
  diff <- dbmap$diffusion$Diffusor(n_components = as.integer(n_components), n_neighbors = as.integer(n_neighbors),
                                   transitions = as.logical(F),
                                   norm = as.logical(F), ann_dist = as.character('cosine'),
                                   n_jobs = as.integer(10), kernel_use = as.character('simple'))$fit(data)
  print('.')
  dbmap_components <- reticulate::py_to_r(diff$transform(data))
  print('.')
  res <- diff$return_dict()
  print(plot(reticulate::py_to_r(res$EigenValues)))
  print(barplot(reticulate::py_to_r(res$EigenValues)))
  print('.')
  rownames(dbmap_components) <- cellnames
  print('.')
  dim.names <- list()
  for(t in 1:length(colnames(dbmap_components))) {
    dim.names[[t]] <- paste0('dbmap_', t)
  }
  colnames(dbmap_components) <- unlist(dim.names)
  reducedDim(object, reduction.save) <- as.matrix(dbmap_components)
  return(object)
}

pancreas <- perform.dbmap(object = pancreas, assay = assay(pancreas, 'sctransform')[metadata(pancreas)$HVGs,], reduction.save = 'uncorrected_dbmap')

perform.umap <- function(object, 
                         reduction='pca',
                         graph = NULL,
                         reduction.save,
                         n.dim, 
                         n_components = 3, 
                         ...) {
  
  c <- uwot::umap(X = reducedDim(object, reduction)[,n.dim], n_components = n_components, verbose = TRUE, ...)
  dim.names <- list()
  for(l in 1:n_components) {
    dim.names[[l]] <- paste0('umap_', l)
  }
  colnames(c) <- unlist(dim.names)
  reducedDim(object, reduction.save) <- c
  
  return(object)
}

pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'scanorama_reduced_umap', 
                         reduction = 'scanorama_reduced', 
                         n.dim = 1:50)
pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'uncorrected_pca_umap', 
                         reduction = 'uncorrected_pca', 
                         n.dim = 1:12)

perform.tsne <- function(object, 
                         reduction, 
                         reduction.save, 
                         n.dim,
                         n_components = 3, 
                         ...) {
  if(isS4(object) == FALSE) {
    cat(crayon::cyan('Must be S4 class object\n'))
  }
  if(is.null(assay)) {
    cat(crayon::cyan(paste0('Please provide assay\n')))
  }
  cat(crayon::cyan('t-SNE reduction initialising\n'))
  c <- ProjectionBasedClustering::tSNE(DataOrDistances = as.matrix(reducedDim(object, reduction))[,n.dim], 
                                       OutputDimension = n_components, Iterations = 1000, ...)$ProjectedPoints
  cat(crayon::cyan('t-SNE reduction completed\n'))
  dim.names <- list()
  for(t in 1:n_components) {
    dim.names[[t]] <- paste0('tsne_', t)
  }
  colnames(c) <- unlist(dim.names)
  reducedDim(object, reduction.save) <- c
  cat(crayon::cyan('t-SNE data added\n'))
  return(object)
}

pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'uncorrected_dbmap_tsne', 
                         reduction = 'uncorrected_dbmap', 
                         n.dim = 1:104)
pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'uncorrected_pca_tsne', 
                         reduction = 'uncorrected_pca', 
                         n.dim = 1:12)

library(harmony)

perform.harmony <- function(object, group.by.vars, reduction = 'pca', dims.use = NULL,
                            theta = NULL, lambda = NULL, sigma = 0.1, nclust = NULL,
                            tau = 0, block.size = 0.05, max.iter.harmony = 10,
                            max.iter.cluster = 20, epsilon.cluster = 1e-05,
                            epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                            reference_values = NULL, reduction.save = "harmony", ...) {
  mat <- reducedDim(object, reduction)
  cat(crayon::cyan('Initialising harmony\n'))
  harm <- harmony::HarmonyMatrix(data_mat = mat[,dims.use], meta_data = colData(object), vars_use = group.by.vars, do_pca = FALSE, 
                                 theta = theta, lambda = lambda, sigma = sigma, nclust = nclust, tau = tau, 
                                 block.size = block.size, max.iter.harmony = max.iter.harmony, max.iter.cluster = max.iter.cluster, 
                                 epsilon.cluster = epsilon.cluster, epsilon.harmony = epsilon.harmony, plot_convergence = plot_convergence, 
                                 return_object = FALSE, verbose = verbose, reference_values = reference_values)
  reducedDim(object, reduction.save) <- harm
  cat(crayon::cyan('Harmony completed\n'))
  return(object)
}

pancreas <- perform.harmony(object = pancreas, group.by.vars = c('tech'), dims.use = 1:12, max.iter.harmony = 100, reduction = 'uncorrected_pca', theta = 1, reduction.save = 'pca_harmony')
pancreas <- perform.harmony(object = pancreas, group.by.vars = c('tech'), dims.use = 1:104, max.iter.harmony = 100, reduction = 'uncorrected_dbmap', theta = 1, reduction.save = 'dbmap_harmony')

perform.scanorama <- function(object, 
                              assay, 
                              split.by, 
                              reduced.dim, 
                              n.dims = 50, 
                              reduction.save='scanorama', 
                              assay.save='scanorama', 
                              batch_size = as.integer(5000), 
                              approx = TRUE, 
                              sigma = as.integer(15), 
                              alpha = as.integer(0.1), 
                              knn = as.integer(20), 
                              hvg = NULL) {
  cat(crayon::cyan('Initialising scanorama\n'))
  scanorama <- reticulate::import('scanorama', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  list.matrix <- list()
  column.names <- list()
  sep <- unique(colData(object)[,split.by])
  print('.')
  mat <- assay
  print('.')
  counter <- 1
  print('.')
  for(x in sep) {
    print('.')
    column.names[[counter]] <- colnames(mat[,colData(object)[,split.by] == x])
    print('.')
    list.matrix[[counter]] <- t(mat[,colData(object)[,split.by] == x])
    print('.')
    counter <- counter + 1
  }
  cat(crayon::cyan('Matrices isolated\n'))
  gene.list <- list()
  for(x in 1:length(sep)) {
    gene.list[[x]] <- rownames(mat[,colData(object)[,split.by] == x])
  }
  cat(crayon::cyan('Genes identified\n'))
  
  cat(crayon::cyan('Corrections starting\n'))
  integrated.corrected.data <- scanorama$correct(datasets_full = reticulate::r_to_py(list.matrix), 
                                                 genes_list = reticulate::r_to_py(gene.list), 
                                                 dimred = as.integer(n.dims), 
                                                 return_dimred=TRUE, 
                                                 return_dense=TRUE, 
                                                 verbose = TRUE, 
                                                 batch_size = as.integer(batch_size), 
                                                 approx = approx, 
                                                 sigma = as.integer(sigma), 
                                                 alpha = as.integer(alpha), 
                                                 knn = as.integer(knn), 
                                                 hvg = hvg)
  
  if(reduced.dim == TRUE) {
    dims <- list()
    cat(crayon::cyan('Isolating scanorama reduced dimensions\n'))
    dim.names <- list()
    for(c in 1:n.dims) {
      dim.names[[c]] <- paste0('scanorama_', c)
    }
    
    dim.names <- unlist(dim.names)
    
    for(x in 1:length(sep)) {
      transposed <- t(reticulate::py_to_r(integrated.corrected.data)[[1]][[x]])
      colnames(transposed) <- column.names[[x]]
      rownames(transposed) <- dim.names
      dims[[x]] <- transposed
    }
    cat(crayon::cyan('Combining samples\n'))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan('Samples concatenated\n'))
    reducedDim(object, reduction.save) <- t(combined)
    
    return(object)
  } else if (reduced.dim == FALSE) {
    dims <- list()
    cat(crayon::cyan('Isolating scanorama corrected gene matrix\n'))
    for(x in 1:length(sep)) {
      transposed <- t(reticulate::py_to_r(integrated.corrected.data)[[2]][[x]])
      print('.')
      colnames(transposed) <- column.names[[x]]
      print('.')
      rownames(transposed) <- reticulate::py_to_r(integrated.corrected.data)[[3]]
      print('.')
      dims[[x]] <- transposed
    }
    cat(crayon::cyan('Combining samples\n'))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan('Samples concatenated\n'))
    reducedDim(object, reduction.save) <- as.matrix(t(combined))
    cat(crayon::cyan('Scanorama completed\n'))
    return(object)
  }
  
}

pancreas <- perform.scanorama(object = pancreas, assay = t(reducedDim(pancreas, 'sctransform_scaled')), split.by = 'tech', reduced.dim = TRUE, n.dims = 50, reduction.save = 'scanorama_reduced')

perform.bbknn <- function(object, 
                          reduced.dim, 
                          dims, 
                          split.by, 
                          graph.save = 'bbknn') {
  anndata <- reticulate::import('anndata', convert = FALSE)
  sc <- reticulate::import('scanpy', convert = FALSE)
  bbknn <- reticulate::import('bbknn', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  pca <- reticulate::r_to_py(reducedDim(object, reduced.dim)[,dims])
  batch <- reticulate::r_to_py(object[[split.by]])
  adata <- anndata$AnnData(X = pca, obs = batch)
  cat(crayon::cyan('Anndata object created\n'))
  sc$tl$pca(adata)
  adata$obsm$X_pca <- pca
  cat(crayon::cyan('Correcting data\n'))
  bbknn$bbknn(adata = adata, batch_key = 0)
  cat(crayon::cyan('Data corrected\n'))
  
  mat  <- reticulate::py_to_r(adata$obsp[['connectivities']])
  rownames(mat) <- colnames(object)
  colnames(mat) <- colnames(object)
  
  metadata(object)[['graphs']][[graph.save]] <- as.Graph(as.matrix(mat))
  cat(crayon::cyan('BBKNN completed\n'))
  return(object)
}

pancreas <- perform.bbknn(object = pancreas, reduced.dim = 'uncorrected_pca', dims = 1:12, split.by = 'tech', dim.save = 'pca_bbknn')
pancreas <- perform.bbknn(object = pancreas.test, reduced.dim = 'uncorrected_dbmap', dims = 1:length(colnames(reducedDim(pancreas.test, 'uncorrected_dbmap'))), split.by = 'tech', graph.save = 'dbmap_bbknn')

pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'pca_harmony_tsne', 
                         reduction = 'pca_harmony', 
                         n.dim = 1:12)
pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'dbmap_harmony_tsne', 
                         reduction = 'dbmap_harmony', 
                         n.dim = 1:104)
pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'pca_bbknn_tsne', 
                         reduction = 'pca_bbknn', 
                         n.dim = 1:12)
pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'dbmap_bbknn_tsne', 
                         reduction = 'dbmap_bbknn', 
                         n.dim = 1:104)
pancreas <- perform.tsne(object = pancreas, 
                         reduction.save = 'scanorama_reduced_tsne', 
                         reduction = 'scanorama_reduced', 
                         n.dim = 1:50)

pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'pca_harmony_umap', 
                         reduction = 'pca_harmony', 
                         n.dim = 1:12)
pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'dbmap_harmony_umap', 
                         reduction = 'dbmap_harmony', 
                         n.dim = 1:104)
pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'pca_bbknn_umap', 
                         reduction = 'pca_bbknn', 
                         n.dim = 1:12)
pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'dbmap_bbknn_umap', 
                         reduction = 'dbmap_bbknn', 
                         n.dim = 1:104)
pancreas <- perform.umap(object = pancreas, 
                         reduction.save = 'scanorama_reduced_umap', 
                         reduction = 'scanorama_reduced', 
                         n.dim = 1:50)

perform.seurat.cluster <- function(object, 
                                   reduction='pca', 
                                   data_frame_name,
                                   res=c(0.1,0.2,0.3,0.4,0.5,
                                         0.6,0.7,0.8,0.9,1,
                                         1.1,1.2,1.3,1.4,1.5), 
                                   dims=1:34,
                                   prune.SNN=0, 
                                   nn.method='annoy', 
                                   annoy.metric='euclidean', 
                                   nn.eps=0.0, 
                                   ...) {
  
  tmp <- suppressWarnings(Seurat::as.Seurat(x = object, counts='counts', data=NULL))
  cat(crayon::cyan('Converted SCE to Seurat object\n'))
  tmp <- Seurat::FindNeighbors(object = tmp, reduction = reduction, verbose = TRUE, dims = dims, compute.SNN = TRUE, prune.SNN = prune.SNN,
                               nn.method = nn.method, annoy.metric = annoy.metric, nn.eps = nn.eps)
  cat(crayon::cyan('Neighbours identified\n'))
  tmp <- Seurat::FindClusters(object = tmp, resolution = res, ...)
  cat(crayon::cyan('Clusters identified\n'))
  temp <- object
  orig.names <- colnames(colData(object))
  new.names <- colnames(tmp@meta.data)
  sep.names <- new.names[!(new.names %in% orig.names)]
  sep.names <-sep.names[1:length(sep.names)-1]
  new.clusters <- tmp@meta.data[,sep.names]
  z <- list()
  for(t in res) {
    z[length(z)+1] <- paste0('Seurat_res_', t)
  }
  
  colnames(new.clusters) <- unlist(z)
  metadata(temp)[['clustering']][[data_frame_name]] <- new.clusters
  cat(crayon::cyan('Seurat clusters added\n'))
  return(temp)
}

# pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'uncorrected_dbmap', dims = 1:104, data_frame_name = 'seurat:uncorrected_dbmap')
# pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'uncorrected_pca', dims = 1:12, data_frame_name = 'seurat:uncorrected_pca')
# pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'dbmap_bbknn', dims = 1:104, data_frame_name = 'seurat:dbmap_bbknn')
# pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'pca_bbknn', dims = 1:12, data_frame_name = 'seurat:pca_bbknn')
pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'scanorama_reduced', dims = 1:50, data_frame_name = 'seurat:scanorama_reduced')
# pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'pca_harmony', dims = 1:12, data_frame_name = 'seurat:pca_harmony')
# pancreas <- perform.seurat.cluster(object = pancreas, reduction = 'dbmap_harmony', dims = 1:104, data_frame_name = 'seurat:dbmap_harmony')


library(SC3)

perform.sc3.cluster <- function(object, 
                                reduction, 
                                dims, 
                                data_frame_name,
                                ks, 
                                n.core=3) {
  cat(crayon::cyan('Initialising SC3 clustering\n'))
  temp.2 <- object
  temp.2 <- SingleCellExperiment(list('logcounts' = t(reducedDim(object, reduction))[dims,]))
  rowData(temp.2)$feature_symbol <- rownames(temp.2)
  temp.2 <- temp.2[!duplicated(rowData(temp.2)$feature_symbol), ]
  temp.2 <- sc3_prepare(temp.2, gene_filter = FALSE, n_cores = n.core)
  temp.2 <- sc3_calc_dists(temp.2)
  temp.2 <- sc3_calc_transfs(temp.2)
  temp.2 <- sc3_kmeans(temp.2, ks = ks)
  temp.2 <- sc3_calc_consens(temp.2)
  orig.names <- colnames(colData(object))
  new.names <- colnames(colData(temp.2))
  sep.names <- new.names[!(new.names %in% orig.names)]
  new.clusters <- colData(temp.2)[,sep.names]
  metadata(object)[['clustering']][[data_frame_name]] <- as.data.frame(new.clusters)
  cat(crayon::cyan('SC3 clustering completed\n'))
  return(object)
}

processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust.sc3clust <- perform.sc3(object = processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust, reduction = 'scanorama', dims = 1:50, ks = 11:15, n.core = 3)

perform.tsne.kmeans <- function(object, 
                                reduction=NULL, 
                                k=NULL,
                                data_frame_name,
                                method='kmeans',
                                ...) {
  
  clusters <- data.frame(kmeans=numeric(length(colnames(object))))
  rownames(clusters) <- colnames(object)
  
  if(is.null(k)) {
    cat(crayon::cyan('Specify number of clusters\n'))
  }
  if(is.null(reduction)) {
    cat(crayon::cyan('Provide assay\n'))
  }
  if(method == 'pam') {
    for(i in k) {
      clusters[,paste0('pam_clustering_K_', i)] <- cluster::pam(x = reducedDim(object, reduction), k = i, ...)$clustering
    }
  }
  if(method == 'kmeans') {
    for(i in k) {
      clusters[[paste0('kmeans_clustering_K_', i)]] <- kmeans(x = reducedDim(object, reduction), centers = i, ...)$cluster
    }
  } else {
    cat(crayon::cyan('Please specify method: pam or kmeans\n'))
  }
  print(clusters)
  print(data_frame_name)
  clusters <- clusters[,2:length(colnames(clusters))]
  metadata(object)[['clustering']][[data_frame_name]] <- clusters
  return(object)
}

pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'uncorrected_dbmap_tsne', k=10:16, data_frame_name = 'tsne_kmeans:uncorrected_dbmap')
pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'uncorrected_pca_tsne', k=10:16, data_frame_name = 'tsne_kmeans:uncorrected_pca')
pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'dbmap_bbknn_tsne', k=10:16, data_frame_name = 'tsne_kmeans:dbmap_bbknn')
pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'pca_bbknn_tsne', k=10:16, data_frame_name = 'tsne_kmeans:pca_bbknn')
pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'scanorama_reduced_tsne', k=10:16, data_frame_name = 'tsne_kmeans:scanorama_reduced')
pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'pca_harmony_tsne', k=10:16, data_frame_name = 'tsne_kmeans:pca_harmony')
pancreas <- perform.tsne.kmeans(object = pancreas, reduction = 'dbmap_harmony_tsne', k=10:16, data_frame_name = 'tsne_kmeans:dbmap_harmony')


library(SingleCellExperiment)

benchmark.clustering <- function(object, 
                                 components, 
                                 reduction, 
                                 dist.method='euclidean',
                                 ground.truth=NULL) {
  
  for(x in names(metadata(object)[['clustering']])) {
    print('.')
    all.clusters <- metadata(object)[['clustering']][[x]]
    print('.')
    print(paste0(unlist(strsplit(x = x, split = ':'))[2], '_', reduction))
    dims <- reducedDim(object, reduction)[,components]
    print('.')
    dist.matrix <- dist(x = dims, method = dist.method)
    print('.')
    sil.results <- data.frame(average_silhoeutte=NA)
    print('.')
    for (v in colnames(all.clusters)[1:length(colnames(all.clusters))]) {
      print(paste0('Calculating silhouette for ', v))
      tmp <- cluster::silhouette(x = as.numeric(x = as.factor(x = all.clusters[,v])), dist = dist.matrix)
      average <- sum(tmp[,3])/length(tmp[,3])
      sil.results[v,] <- average
    }
    sil.results <- sil.results[complete.cases(sil.results),]
    max.AS <- max(sil.results)
    print(max.AS)
    
    dunn.results <- data.frame(dunn.index=NA)
    for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
      print(paste0('Calculating dunn index for ', p))
      dunn.results[p,] <- clValid::dunn(distance = dist.matrix, clusters = as.numeric(x = as.factor(x = all.clusters[,p])))
    }
    dunn.results <- dunn.results[complete.cases(dunn.results),]
    max.dunn <- max(dunn.results)
    print(max.dunn)
    
    conn.results <- data.frame(connectivity=NA)
    for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
      print(paste0('Calculating connectivity for ', p))
      conn.results[p,] <- clValid::connectivity(distance = dist.matrix, clusters = all.clusters[,p])
    }
    conn.results <- conn.results[complete.cases(conn.results),]
    max.conn <- max(conn.results)
    print(max.conn)
    
    if(!is.null(ground.truth)) {
      ARI.results <- data.frame(ARI=NA)
      for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
        print(paste0('Calculating ARI for ', p))
        ARI.results[p,] <- mclust::adjustedRandIndex(x = all.clusters[,p], y = ground.truth)
      }
      ARI.results <- ARI.results[complete.cases(ARI.results),]
      max.ARI <- max(ARI.results)
      print(max.ARI)
      NMI.results <- data.frame(NMI=NA)
      for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]) {
        print(paste0('Calculating NMI for ', p))
        NMI.results[p,] <- aricode::AMI(c1 = all.clusters[,p], c2 = ground.truth)
      }
      NMI.results <- NMI.results[complete.cases(NMI.results),]
      results <- cbind(sil.results, dunn.results, conn.results, ARI.results, NMI.results)
      rownames(results) <- colnames(all.clusters)
      colnames(results) <- c(paste0(x, '_sil.results'), paste0(x, '_dunn.results'), paste0(x, '_conn.results'), paste0(x, '_ARI.results'), paste0(x, '_NMI.results'))
      metadata(object)[['benchmarking_clustering']][[as.character(x)]] <- results
    } else {
      results <- cbind(sil.results, dunn.results, conn.results)
      rownames(results) <- colnames(all.clusters)
      colnames(results) <- c(paste0(x, '_sil.results'), paste0(x, '_dunn.results'), paste0(x, '_conn.results'))
      metadata(object)[['benchmarking_clustering']][[as.character(x)]] <- results
    }
  }
  return(object)
}

pancreas.test <- benchmark.clustering(object = pancreas, 
                                      components = 1:3, 
                                      reduction = 'umap', 
                                      dist.method = 'euclidean', 
                                      ground.truth = pancreas$celltype)

library(plotly)
library(ggplot2)
library(egg)

panc <- readRDS("~/Results/scRNA-seq/IBRAP_development/Example_processed_results/RShiny_test.RData")

plot.reduced.dim <- function(object, 
                             reduction='', 
                             pt.size=5, 
                             metadata.access='clustering',
                             sub.access='metadata',
                             group.by, 
                             dimensions) {
  
  if(is.null(metadata(object)[['clustering']][['metadata']])){
    metadata(object)[['clustering']][['metadata']] <- as.data.frame(colData(object))
  }
  
  print('plot_cluster_dr_started')
  results <- as.data.frame(reducedDim(object, as.character(reduction)))[,1:3]
  print('.')
  results[,'variable'] <- metadata(object)[[metadata.access]][[sub.access]][[group.by]]
  rownames(results) <- colnames(object)
  print('.')
  if(dimensions == 3) {
    print('3')
    print(plotly::plot_ly(data = results,
                          x = as.formula(paste0('~', colnames(results)[1])), 
                          y = as.formula(paste0('~', colnames(results)[2])),
                          z = as.formula(paste0('~', colnames(results)[3])), 
                          color = as.formula(paste0('~',colnames(results)[4])),  
                          mode = "markers", 
                          marker = list(size = pt.size, width=0.5), 
                          text=as.formula(paste0('~',colnames(results)[4])), 
                          hoverinfo="text", plot_bgcolor = 'black'))
    
  } else if (dimensions == 2) {
    print('2')
    print(plotly::plot_ly(data = as.data.frame(results), 
                          x = as.formula(paste0('~', colnames(results)[1])), 
                          y = as.formula(paste0('~', colnames(results)[2])), 
                          color = as.formula(paste0('~',colnames(results)[4])),  
                          mode = "markers", 
                          marker = list(size = pt.size, width=0.5), 
                          text=as.formula(paste0('~',colnames(results)[4])), 
                          hoverinfo="text", plot_bgcolor = 'black'))
  }
  
}

plot.reduced.dim(object = pancreas, reduction = 'scanorama_reduced_umap',
                metadata.access = 'clustering', sub.access = 'seurat:scanorama_reduced',
                group.by = 'Seurat_res_0.2', dimensions = 2)

plot.features <- function(object, 
                          reduction='', 
                          pt.size=10, 
                          assay, 
                          feature,
                          dimensions) {
  results <- as.data.frame(reducedDim(object, reduction))[,1:3]
  print(dimensions)
  iso <- assay(object, assay)[feature,]
  print('.')
  results[,feature] <- iso
  print('.')
  colnames(results)[4] <- 'feature'
  print('.')
  if(dimensions == 3){
    print('3')
    p <- plotly::plot_ly(data = results[order(results$feature),], 
                         x = as.formula(paste0('~', colnames(results)[1])), 
                         y = as.formula(paste0('~', colnames(results)[2])),
                         z = as.formula(paste0('~', colnames(results)[3])), 
                         color = as.formula(paste0('~',colnames(results)[4])), 
                         mode = "markers", colors = RColorBrewer::brewer.pal(n = 9, name = 'Blues')[3:9],
                         marker = list(size = 5, width=5), 
                         text=as.formula(paste0('~',colnames(results)[4])), 
                         hoverinfo="text", plot_bgcolor = 'black')
    p <- p %>% layout(title=feature)
    print(p)
  } 
  if (dimensions == 2) {
    print('2')
    p <- plotly::plot_ly(data = results[order(results$feature),], 
                         x = as.formula(paste0('~', colnames(results)[1])), 
                         y = as.formula(paste0('~', colnames(results)[2])),
                         color = as.formula(paste0('~',colnames(results)[4])),  
                         mode = "markers",
                         colors = RColorBrewer::brewer.pal(n = 9, name = 'Blues')[3:9],
                         marker = list(size = 5, width=5), 
                         text=as.formula(paste0('~',colnames(results)[4])), 
                         hoverinfo="text", plot_bgcolor = 'black')
    p <- p %>% layout(title=feature)
    print(p)
  } 
}

#plot.features(object = panc, reduction = 'scanorama_reduced_umap', 
#              pt.size = 1, assay = 'decontXcounts', feature = "A1BG-AS1", 
#              dimensions = 2)

plot.features.multiple <- function(object, 
                                   assay, 
                                   reduction, 
                                   lab.key, 
                                   features) {
  
  plot.list <- list()
  
  for(x in features) {
    print('1')
    results <- as.data.frame(reducedDim(object, reduction))[,1:2]
    print('2')
    print(assay)
    iso <- assay(object, assay)[x,]
    print('3')
    colnames(results) <- c('red_1', 'red_2')
    print('4')
    results[,x] <- iso
    print('5')
    colnames(results)[3] <- 'feature'
    print('6')
    plot.list[[x]] <- ggplot(data = results[order(results$feature),], aes(x = red_1, y = red_2)) + 
      geom_point(aes(color=feature)) + 
      scale_color_gradient(low = '#C6DBEF', high = '#08306B') + 
      theme_bw() + labs(title=x, x=paste0(lab.key,'_1'), y=paste0(lab.key,'_2')) + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10))
    print('7')
  }
  
  is.even <- function(x) {
    if(as.integer(x %% 2) == 0) {
      print(TRUE)
    } else {
      print(FALSE)
    }
  }
  
  if(!is.even(length(plot.list))) {
    plot.list[length(plot.list)+1] <- plot.list[1] + geom_blank()
  }
  
  if(length(plot.list) > 1){
    print('multi')
    if(length(plot.list) <3) {
      print('1x2')
      do.call('ggarrange', c(plots = plot.list, ncol = 1, nrow = 2))
    } else if(length(plot.list) <5) {
      print('2x2')
      do.call('ggarrange', c(plots = plot.list, ncol = 2, nrow = 2))
    } else if(length(plot.list) <7) {
      print('3x2')
      do.call('ggarrange', c(plots = plot.list, ncol = 3, nrow = 2))
    }
  } else {
    print('single')
    plot.list[1]
  }
}

plot.features.multiple(object = panc, assay = 'decontXcounts', 
                       reduction = 'scanorama_reduced_umap', 
                       lab.key = 'umap', features = features[1:3])

plot.barplot <- function(object, 
                         x.value, 
                         y.value) {
  object[['var']] <- y.value 
  object[['group']] <- x.value
  print(dittoSeq::dittoBarPlot(object = object, var = 'var', group.by = 'group') + labs(title='') + ggplot2::scale_fill_hue())
}

plot.barplot(object = panc, y.value = metadata(panc)$clustering[['seurat:pca_harmony']][['Seurat_res_0.2']], x.value = panc$tech)

plot.benchmarking <- function(object, 
                              clust.method, 
                              ARI){
  print('.')
  clust.bench <- metadata(object)[['benchmarking_clustering']][[as.character(clust.method)]]
  print('.')
  clust.bench <- as.data.frame(clust.bench)
  print(clust.bench)
  clust.bench[,'cluster_index'] <- rownames(clust.bench)
  if(ARI == TRUE) {
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'ARI', 'NMI', 'cluster_index')
  } else {
    labels <- c('ASW', 'Dunn_index', 'Connectivity', 'cluster_index')
  }
  
  colnames(clust.bench) <- labels
  print(clust.bench)
  list.plot <- list()
  
  for(o in 1:sum(length(labels)-2)) {
    
    label <- labels[as.numeric(o)]
    print(label)
    fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(label), group = 1)) +
      geom_point() +
      geom_line() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    list.plot[[as.numeric(o)]] <- fig
    
  }
  
  last.label <- labels[as.numeric(sum(length(labels)-1))]
  
  last.fig <- fig <- ggplot(clust.bench, aes_string(x = 'cluster_index', y = as.character(last.label), group = 1)) +
    geom_point() +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  list.plot[[as.numeric(sum(length(labels)-1))]] <- last.fig
  print('/')
  do.call('ggarrange', c(plots = list.plot, ncol = 5))
}

plot.heatmap <- function(object, 
                         assay, 
                         features, 
                         ...) {
  
  ass <- object[features,]
  ass <- assay(ass, assay)
  cat(crayon::cyan('Isolated assay\n'))
  
  logged <- log(ass+1)
  means <- apply(X = logged, MARGIN = 1, FUN = mean)
  standev <- apply(X = logged, MARGIN = 1, FUN = sd)
  z_scores <- logged - means / standev 
  cat(crayon::cyan('z-scores calculated\n'))
  
  p <- gplots::heatmap.2(x = z_scores, 
                         key.title = 'colour key', 
                         key.xlab = 'z-score', 
                         key.ylab = '', xlab = 'cells', 
                         ylab = 'genes', trace = 'none', margins = c(7,7))
  p
}

plot.heatmap(object = pancreas, assay = 'counts', features = features)

plot.vln <- function(object, 
                     assay,
                     features, 
                     group.by, 
                     title = NULL, 
                     xlab = 'group', 
                     ylab = 'expression') {
  
  if(!length(group.by) == length(colnames(object))) {
    cat(crayon::cyan('group.by variable must match the length of column names in SCE.'))
    return(NULL)
  }
  
  plot.list <- list()
  
  for(x in features) {
    ass <- t(assay(object[x,], assay))
    df <- data.frame(groups = colnames(object))
    rownames(df) <- df$barcodes
    df <- cbind(df, ass)
    df[,'groups'] <- group.by
    colnames(df) <- c('group', 'feature')
    
    p <- ggplot(data = df, aes(x = group, y = feature, color = group)) + 
      geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + 
      theme_light() + labs(title = x, y = 'expression', x = 'group') + 
      theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 20))
    
    plot.list[[x]] <- p
  }

  if(length(plot.list) > 1) {
    do.call('ggarrange', c(plots = plot.list, ncol = 1))
  } else {
    plot.list[1]
  }
}

plot.vln(object = panc, 
         assay = 'sctransform', 
         features = c("A1BG-AS1", "A1BG", "A1CF"),
         group.by = metadata(panc)[['clustering']][['seurat:scanorama_reduced']][['Seurat_res_0.2']])

perform.seurat.diffexp <- function(object, counts.assay, norm.assay, scaled.assay, test, cluster.method, cluster.column, ...) {
  seuobj <- Seurat::as.Seurat(x = object, counts = counts.assay, data = norm.assay)
  seuratobj$clusters <- metadata(pancreas)[['clustering']][['seurat:scanorama_reduced']][['Seurat_res_0.2']]
  Seurat::Idents(seuratobj) <- 'clusters'
  results <- Seurat::FindAllMarkers(object = seuobj, test.use = test)
  return(results)
}

tester <- perform.seurat.diffexp(object = panc, 
                                 counts.assay = 'sctransform_counts', 
                                 norm.assay = 'sctransform_data', 
                                 test = 'MAST', 
                                 cluster.method = 'seurat:scanorama_reduced', 
                                 cluster.column = 'Seurat_res_0.2')





