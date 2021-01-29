Koh <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Koh.SCE.rds')
Kumar <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Kumar.SCE.rds')
Zhengmix4eq <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Zhengmix4eq.SCE.rds')
Zhengmix4uneq <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Zhengmix4uneq.SCE.rds')
Zhengmix8eq <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/Zhengmix8eq.SCE.rds')
mixology10x3cl <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/mixology10x3cl.SCE.rds')
mixology10x5cl <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/mixology10x5cl.SCE.rds')
simMix1 <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/simMix1.SCE.rds')
simMix2 <- readRDS('~/Results/scRNA-seq/IBRAP/IBRAP/datasets/simMix2.SCE.rds')
pancreas.data <- readRDS(file = "~/Raw_Data/pancreas_v3_files/pancreas_expression_matrix.rds")
metadata <- readRDS('~/Raw_Data/pancreas_v3_files/pancreas_metadata.rds')
pancreas.data <- as.matrix(pancreas.data)
metadata <- readRDS(file = "/Users/knight05/Downloads/pancreas_v3_files/pancreas_metadata.rds")

reticulate:py_config()

library(utils)
library(Matrix)

# Read 10x

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

library(SingleCellExperiment)

# Metadata generator

metadata.generator <- function(object, 
                               assay, 
                               prefix) {
  cat(crayon::cyan('Annotated with fundamental metadata for assay: ', 'assay', '\n'))
  colData(object)[[paste0(prefix, '_total.counts')]] <- colSums(assay(object, assay))
  colData(object)$total.features <- colSums(assay(object, assay) != 0)
  cat(crayon::cyan('Completed\n'))
}

# SCE

create_sce_object <- function(counts, 
                              project.name, 
                              min.cells, 
                              metadata) {
  sce <- SingleCellExperiment(assays = list(counts = counts))
  cat(crayon::cyan('SCE object created\n'))
  colData(sce)$original.project <- project.name
  sce <- metadata.generator(object = sce, assay = 'counts', prefix = 'raw_')
  if(!is.null(metadata))  {
    colData(sce) <- cbind(colData(sce), metadata)
    cat(crayon::cyan('Custom metadata appended\n'))
  }
  sce <- sce[!duplicated(rownames(assay(sce, 'counts'))), ]
  cat(crayon::cyan('Completed\n'))
  return(sce)
}

pancreas <- create_sce_object(counts = pancreas.data, 
                              project.name = 'pancreas', 
                              min.cells = 3, 
                              metadata = metadata)

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
                     assay = 'decontXcounts', 
                     prefix = 'decontaminated_')
  
  cat(crayon::cyan('Finished\n'))
  return(object)
}

pancreas <- perform.decontX(object = pancreas, 
                            batch = 'tech')



find_percentage_genes <- function(object, 
                                  pattern='^MT-', 
                                  which.assay='counts', 
                                  prefix) {
  temp <- data.frame(tmp=colSums(assay(object, which.assay)[grep(pattern = pattern, x = rownames(assay(object, which.assay))),]) / colData(object)$total.counts * 100)
  colnames(temp) <- c(column.name)
  if(column.name %in% colnames(colData(object))) {
    colData(object) <- colData(object)[,colnames(colData(object)) != column.name]
  }
  colData(object) <- cbind(colData(object), temp)
  return(object)
}

pancreas <- find_percentage_genes(object = pancreas, pattern = '^MT-', which.assay = 'decontXcounts', column.name = 'percent.mt')

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

sd.value <- sd(pancreas$total.features)
med.value <- median(pancreas$total.features)
max.features <- (sd.value*3)+med.value

pancreas <- filter_features_cells(object = pancreas, 
                                  use.assay = 'decontXcounts', 
                                  max.features = max.features, 
                                  min.features = 200, 
                                  max.mt = 10, 
                                  min.cell = 3)

add.cell.cycle <- function(object, 
                           assay, 
                           transform, ...) {
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/Homo_sapiens.csv', header = TRUE, sep = ',')
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

test <- add.cell.cycle(object = pancreas.scran, assay = 'decontXcounts', transform = TRUE)

add.feature.score <- function(object, 
                              assay, 
                              transform, 
                              features, 
                              ...) {
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP_development/IBRAP/Homo_sapiens.csv', header = TRUE, sep = ',')
  if(transform == TRUE) {
    seuobj <- Seurat::CreateSeuratObject(counts = assay(object, assay))
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::NormalizeData(object = seuobj)
    cat(crayon::cyan('Data transformed\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  } else {
    seuobj <- as.Seurat(object, counts = NULL, data = assay)
    cat(crayon::cyan('Converted to Seurat object\n'))
    seuobj <- Seurat::AddModuleScore(object = seuobj, features = features, ...)
    cat(crayon::cyan('Seurat gene score calculated\n'))
    colData(object) <- cbind(colData(object), seuobj@meta.data[, length(colnames(seuobj@meta.data))])
    cat(crayon::cyan('New metadata added\n'))
  }
  return(object)
}

library(Seurat)
library(sctransform)

perform.sct.normalisation <- function(object, 
                                      split.by=NULL, 
                                      which.assay, 
                                      new.assay = 'sctransform', 
                                      ...) {
  if(is.null(split.by)) {
    seuratobj <- CreateSeuratObject(counts = assay(object, which.assay), project = 'NA')
    seuratobj <- SCTransform(object = seuratobj, ...)
    genes <- intersect(rownames(object), rownames(as.matrix(seuratobj@assays$SCT@data)))
    object <- object[genes,]
    assay(object, new.assay) <- as.matrix(seuratobj@assays$SCT@data)[genes,]
    return(object)
  } else {
    t <-unique(colData(object)[[split.by]])
    list.matrix <- list()
    cat(crayon::cyan('loading datasets\n'))
    for(o in t) {
      cat(crayon::cyan(paste0('calculating: ', o, '\n')))
      sub <- object[,colData(object)[[split.by]] == o]
      seuratobj <- CreateSeuratObject(counts = assay(sub, which.assay), project = 'NA')
      seuratobj <- SCTransform(object = seuratobj, ...)
      list.matrix[[o]] <- seuratobj@assays$SCT@data
    }
    cat(crayon::cyan('Aligning features: stage 1\n'))
    genes.length <- list()
    for(o in t) {
      genes.length[[o]] <- length(rownames(list.matrix[[o]]))
    }
    
    cat(crayon::cyan('Aligning features: stage 2\n'))
    highest.rows <- names(which.max(rank(x = unlist(genes.length))))
    gene.list <- rownames(list.matrix[[highest.rows]])
    
    for(o in t[!t %in% highest.rows]){
      gene.list <- intersect(gene.list, rownames(list.matrix[[o]]))
    }
    
    cat(crayon::cyan('Aligning features: stage 3\n'))
    
    for(o in t) {
      list.matrix[[o]] <- list.matrix[[o]][gene.list,]
    }
    
    cat.mat <- do.call('cbind', list.matrix)
    cat(crayon::cyan('Aligning features: complete\n'))
    empty <- object
    common_rows <- intersect(rownames(cat.mat), rownames(empty))
    empty <- empty[common_rows,]
    cat.mat <- cat.mat[common_rows,]
    
    cat(crayon::cyan('Matrix added to S4 object\n'))
    
    assay(empty, new.assay) <- cat.mat
    cat(crayon::cyan('Complete\n'))
    return(empty)
  }
}

pancreas.test.1 <- perform.sct.normalisation(object = pancreas, 
                                      which.assay = 'decontXcounts', 
                                      split.by = NULL)

library(scran)
library(scater)

perform.scran.normalisation <- function(object, 
                                        new.assay = 'scran', 
                                        max.cluster.size = 1000, 
                                        scaling=NULL,
                                        do.log=TRUE, 
                                        center_size_factors=TRUE) {
  
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
  new <- assay(object, 'decontXcounts') 
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

perform.tpm <- function(object, 
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
                                 unwanted.variance, 
                                 scale, 
                                 centre, 
                                 ...) {
  y <- object
  tmp <- suppressWarnings(Seurat::as.Seurat(y[metadata(y)$HVGs,], counts = NULL, data = use.assay))
  cat(crayon::cyan('SCE converted to seurat object\n'))
  tmp <- Seurat::ScaleData(tmp, vars.to.regress = unwanted.variance, do.scale = scale, do.center = scale, ...)
  cat(crayon::cyan('Data scaled\n'))
  reducedDim(y, 'scaled') <- t(as.matrix(tmp@assays$RNA@scale.data))
  cat(crayon::cyan('Scaled matrix attached to object\n'))
  return(y)
}

pancreas.test.1 <- perform.seurat.scale(object = pancreas.test.1, use.assay = 'sctransform', unwanted.variance = c('percent.mt'), scale = TRUE, centre = TRUE)

perform.pca <- function(object, 
                        assay, 
                        reduction.save='pca', 
                        ...) {
  cat(crayon::cyan('Initialising PCA\n'))
  a <- PCAtools::pca(mat = t(reducedDim(object, 'scaled')), center = FALSE, scale = FALSE, ...)
  reducedDim(object, reduction.save) <- as.matrix(a$rotated[,1:50])
  cat(crayon::cyan('PCA completed\n'))
  return(object)
}

pancreas.seurat <- perform.pca(object = pancreas.seurat, assay = reducedDim(pancreas.scran, 'scanorama'))
pancreas.scran <- perform.pca(object = pancreas.scran, scaled_data = TRUE, use.assay = NULL)

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

pancreas.scran <- perform.dbmap(object = pancreas.scran, assay = assay(pancreas.scran)[metadata(pancreas.scran)$HVGs,], reduction.save = 'uncorrected_dbmap')

perform.umap <- function(object, 
                         reduction='pca', 
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

pancreas.scran <- perform.umap(object = pancreas.scran, 
                       reduction.save = 'scanorama_reduced_umap', 
                       reduction = 'scanorama_reduced', 
                       n.dim = 1:50)
pancreas.scran <- perform.umap(object = pancreas.scran, 
                               reduction.save = 'uncorrected_pca_umap', 
                               reduction = 'uncorrected_pca', 
                               n.dim = 1:12)
pancreas.scran <- perform.umap(object = pancreas.scran, 
                               reduction.save = 'scanorama_pca_umap', 
                               reduction = 'scanorama_pca', 
                               n.dim = 1:12)
pancreas.scran <- perform.umap(object = pancreas.scran, 
                               reduction.save = 'scanorama_dbmap_umap', 
                               reduction = 'scanorama_dbmap', 
                               n.dim = 1:104)
pancreas.scran <- perform.umap(object = pancreas.scran, 
                               reduction.save = 'uncorrected_dbmap_umap', 
                               reduction = 'uncorrected_dbmap', 
                               n.dim = 1:104)

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

pancreas.scran <- perform.tsne(object = pancreas.scran, 
                               reduction.save = 'scanorama_reduced_umap', 
                               reduction = 'scanorama_reduced', 
                               n.dim = 1:50)
pancreas.scran <- perform.tsne(object = pancreas.scran, 
                               reduction.save = 'uncorrected_pca_umap', 
                               reduction = 'uncorrected_pca', 
                               n.dim = 1:12)
pancreas.scran <- perform.tsne(object = pancreas.scran, 
                               reduction.save = 'scanorama_pca_umap', 
                               reduction = 'scanorama_pca', 
                               n.dim = 1:12)
pancreas.scran <- perform.tsne(object = pancreas.scran, 
                               reduction.save = 'scanorama_dbmap_umap', 
                               reduction = 'scanorama_dbmap', 
                               n.dim = 1:104)
pancreas.scran <- perform.tsne(object = pancreas.scran, 
                               reduction.save = 'uncorrected_dbmap_umap', 
                               reduction = 'uncorrected_dbmap', 
                               n.dim = 1:104)

library(harmony)

perform.harmony <- function(object, group.by.vars, reduction = 'pca', dims.use = NULL,
                            theta = NULL, lambda = NULL, sigma = 0.1, nclust = NULL,
                            tau = 0, block.size = 0.05, max.iter.harmony = 10,
                            max.iter.cluster = 20, epsilon.cluster = 1e-05,
                            epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                            reference_values = NULL, reduction.save = "harmony", ...) {
  mat <- reducedDim(object, reduction)
  cat(crayon::cyan('Initialising harmony\n'))
  harm <- HarmonyMatrix(data_mat = mat, meta_data = colData(object), vars_use = group.by.vars, do_pca = FALSE, 
                        theta = theta, lambda = lambda, sigma = sigma, nclust = nclust, tau = tau, 
                        block.size = block.size, max.iter.harmony = max.iter.harmony, max.iter.cluster = max.iter.cluster, 
                        epsilon.cluster = epsilon.cluster, epsilon.harmony = epsilon.harmony, plot_convergence = plot_convergence, 
                        return_object = FALSE, verbose = verbose, reference_values = reference_values)
  reducedDim(object, reduction.save) <- harm
  cat(crayon::cyan('Harmony completed\n'))
  return(object)
}

DLBCL <- perform.harmony(object = DLBCL, group.by.vars = c('original.project'), dims.use = 1:18, reduction = 'pca', theta = 1)

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
    reducedDim(object, reduction.save) <- t(combined)
    cat(crayon::cyan('Scanorama completed\n'))
    return(object)
  }
  
}

pancreas.scran <- perform.scanorama(object = pancreas.scran, assay = t(reducedDim(pancreas.scran, 'scaled')), split.by = 'tech', reduced.dim = TRUE, n.dims = 50, reduction.save = 'scanorama_reduced')

pancreas.scran <- perform.pca(object = pancreas.scran, assay = t(reducedDim(pancreas.scran, 'scanorama_matrix')),reduction.save = 'scanorama_pca')

perform.bbknn <- function(object, 
                          reduced.dim, 
                          dims, 
                          split.by, 
                          dim.save = 'bbknn') {
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
  reducedDim(object, as.character(dim.save)) <- reticulate::py_to_r(adata$X)
  cat(crayon::cyan('BBKNN completed\n'))
  return(object)
}

processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.harmony.bbknn <- perform.bbknn(object = processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.harmony, reduced.dim = 'pca', dims = 1:18, split.by = 'tech', dim.save = 'bbknn')

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

pancreas.scran <- perform.seurat.cluster(object = pancreas.scran, reduction = 'scanorama_reduced', dims = 1:50, data_frame_name = 'scanorama_reduced')
pancreas.scran <- perform.seurat.cluster(object = pancreas.scran, reduction = 'uncorrected_dbmap', dims = 1:104, data_frame_name = 'uncorrected_dbmap')
pancreas.scran <- perform.seurat.cluster(object = pancreas.scran, reduction = 'scanorama_dbmap', dims = 1:104, data_frame_name = 'scanorama_dbmap')
pancreas.scran <- perform.seurat.cluster(object = pancreas.scran, reduction = 'uncorrected_pca', dims = 1:12, data_frame_name = 'uncorrected_pca')
pancreas.scran <- perform.seurat.cluster(object = pancreas.scran, reduction = 'scanorama_pca', dims = 1:12, data_frame_name = 'scanorama_pca')

library(SC3)

perform.sc3 <- function(object, 
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
  metadata(object)[['clustering']][[data_frame_name,]] <- as.data.frame(new.clusters)
  cat(crayon::cyan('SC3 clustering completed\n'))
  return(object)
}

processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust.sc3clust <- perform.sc3(object = processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust, reduction = 'scanorama', dims = 1:50, ks = 11:15, n.core = 3)

perform.tsne.kmeans <- function(object, 
                                reduction=NULL, 
                                k=NULL,
                                data_frame_name,
                                ...) {
  if(is.null(k)) {
    cat(crayon::cyan('Specify number of clusters\n'))
  }
  if(is.null(tsne.assay)) {
    cat(crayon::cyan('Provide assay\n'))
  }
  if(method == 'pam') {
    for(i in k) {
      object[[paste0('pam_clustering_K_', i)]] <- cluster::pam(x = reducedDim(object, reduction), k = i, ...)$clustering
    }
  }
  if(method == 'kmeans') {
    for(i in k) {
      clusters[[paste0('kmeans_clustering_K_', i)]] <- kmeans(x = reducedDim(object, reduction), ...)$cluster
    }
  } else {
    cat(crayon::cyan('Please specify method: pam or kmeans\n'))
  }
  object[['clustering']][[data_frame_name]] <- clusters
  return(object)
}

library(SingleCellExperiment)

benchmark.clustering <- function(object, 
                                    components, 
                                    reduction, 
                                    dist.method='euclidean',
                                    ground.truth=NULL) {

  for(x in names(metadata(object)[['clustering']])) {
    all.clusters <- metadata(object)[['clustering']][[x]]
    dims <- reducedDim(object, paste0(x, '_', reduction))[,components]
    dist.matrix <- dist(x = dims, method = dist.method)
    sil.results <- data.frame(average_silhoeutte=NA)
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

pancreas.scran <- benchmark.clustering(object = pancreas.scran, 
                                       components = 1:3, 
                                       reduction = 'umap', 
                                       dist.method = 'euclidean', 
                                       ground.truth = pancreas.scran$celltype)

library(plotly)
library(ggplot2)

plot.reduced.dim <- function(object, 
                            reduction='', 
                            pt.size=5, 
                            metadata.access='clustering',
                            sub.access='metadata',
                            group.by, 
                            dimensions) {
  if(is.null(metadata(object)[['clustering']][['metadata']])){
    metadata(object) <- colData(object)
  }
  l <- metadata(object)[[metadata.access]][[sub.access]][[group.by]]
  print('plot_cluster_dr_started')
  reduction <- as.data.frame(reducedDim(object, as.character(reduction)))
  print('.')
  barcodes <- colnames(object)
  print('.')
  if(dimensions == 3) {
    print('.')
    x.val <- reduction[,1]
    print('.')
    y.val <- reduction[,2]
    print('.')
    z.val <- reduction[,3]
    print('.')
    results <- cbind(l, x.val, y.val, z.val)
    print('.')
    results <- as.data.frame(results)
    print('.')
    results[,1] <- as.factor(results[,1])
    print('.')
    rownames(results) <- barcodes
    print('.')
    if(sub.access != 'metadata') {
      print('.')
      colnames(results) <- c('variable', colnames(reduction)[1:3])
    } else {
      print('.')
      colnames(results) <- c('clusters', colnames(reduction)[1:3])
    }
    print('.')
    print(plotly::plot_ly(data = results, 
                    x = as.formula(paste0('~', colnames(results)[2])), 
                    y = as.formula(paste0('~', colnames(results)[3])),
                    z = as.formula(paste0('~', colnames(results)[4])), 
                    color = as.formula(paste0('~',colnames(results)[1])),  
                    mode = "markers", 
                    marker = list(size = pt.size, width=0.5), 
                    text=as.formula(paste0('~',colnames(results)[1])), 
                    hoverinfo="text", plot_bgcolor = 'black'))
                    
  } else if (dimensions == 2) {
    print('.')
    x.val <- reduction[,1]
    print('.')
    y.val <- reduction[,2]
    print('.')
    results <- cbind(l, x.val, y.val)
    print('.')
    rownames(results) <- barcodes
    print('.')
    results <- as.data.frame(results)
    print('.')
    results[,1] <- as.factor(results[,1])
    print('.')
    if(sub.access != 'metadata') {
      print('.')
      colnames(results) <- c('variable', colnames(reduction)[1:2])
    } else {
      print('.')
      colnames(results) <- c('clusters', colnames(reduction)[1:2])
    }
    print('.')
    print(plotly::plot_ly(data = as.data.frame(results), 
                          x = as.formula(paste0('~', colnames(results)[2])), 
                          y = as.formula(paste0('~', colnames(results)[3])), 
                          color = as.formula(paste0('~',colnames(results)[1])),  
                          mode = "markers", 
                          marker = list(size = pt.size, width=0.5), 
                          text=as.formula(paste0('~',colnames(results)[1])), 
                                          hoverinfo="text", plot_bgcolor = 'black'))
  }

}

plot.reduced.dim(object = x, reduction = 'scanorama_reduced_umap', 
                 metadata.access = 'clustering', sub.access = 'scanorama_reduced',  
                 group.by = 'Seurat_res_0.1', dimensions = 2)

plot.features <- function(object, 
                          reduction='', 
                          pt.size=10, 
                          assay, 
                          feature,
                          dimensions) {
  results <- as.data.frame(reducedDim(object, reduction))[,1:3]
  print('.')
  iso <- assay(object, assay)[feature,]
  print('.')
  results[,feature] <- iso
  print('.')
  colnames(results)[4] <- feature
  print('.')
  if(dimensions == 3){
    print('3')
    print(plotly::plot_ly(data = results, 
                          x = as.formula(paste0('~', colnames(results)[1])), 
                          y = as.formula(paste0('~', colnames(results)[2])),
                          z = as.formula(paste0('~', colnames(results)[3])), 
                          color = as.formula(paste0('~',colnames(results)[4])), 
                          mode = "markers", colors = RColorBrewer::brewer.pal(n = 9, name = 'Blues')[3:9],
                          marker = list(size = 5, width=5), 
                          text=as.formula(paste0('~',colnames(results)[4])), 
                          hoverinfo="text", plot_bgcolor = 'black'))
    print(p)
  } 
  if (dimensions == 2) {
    print('2')
    p <- plotly::plot_ly(data = results, 
                          x = as.formula(paste0('~', colnames(results)[1])), 
                          y = as.formula(paste0('~', colnames(results)[2])),
                          color = as.formula(paste0('~',colnames(results)[4])),  
                          mode = "markers",
                          colors = RColorBrewer::brewer.pal(n = 9, name = 'Blues')[3:9],
                          marker = list(size = 5, width=5), 
                          text=as.formula(paste0('~',colnames(results)[4])), 
                          hoverinfo="text", plot_bgcolor = 'black')
    print(p)

  } else {
    cat(crayon::cyan('Dimensions must be either 2 or 3'))
    return(NULL)
  }
}

plot.features(object = x, reduction = 'scanorama_reduced_umap', 
              pt.size = 1, assay = 'decontXcounts', feature = 'AACSP1', 
              dimensions = 3)

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










