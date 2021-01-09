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
#mode(pancreas.data) <- 'integer'
metadata <- readRDS(file = "/Users/knight05/Downloads/pancreas_v3_files/pancreas_metadata.rds")

#Sys.setenv(RETICULATE_PYTHON = "/Users/knight05/Library/r-miniconda/envs/r-reticulate/bin/python")
#reticulate::py_config()

library(utils)
library(Matrix)

# Read 10x

Read10X_output <- function(directory) {
  count.file <- paste0(directory, '/matrix.mtx')
  genes.file <- paste0(directory, '/genes.tsv')
  barcodes.files <- paste0(directory, '/barcodes.tsv')
  cat(crayon::cyan('Files read\n'))
  genes_ensembl <- read.table(file = genes.file, sep = '\t', header = FALSE)
  barcodes <- read.table(file = barcodes.files, sep = '\t', header = FALSE)
  mm <- Matrix::readMM(count.file)
  rownames(mm) <- genes_ensembl$V2
  colnames(mm) <- barcodes$V1
  cat(crayon::cyan('Success: Matrix concatenated\n'))
  mm <- as.matrix(mm)
  return(mm)
}

R3104 <- Read10X_output(directory = '~/Raw_Data/DLBCL_analysis_Jess/new/R3104')
marrow_A <- Read10X_output(directory = '/Users/knight05/Raw_Data/Database_samples/healthy_references/BMMC_atlas/marrow_A')

library(SingleCellExperiment)

# SCE

create_sce_object <- function(counts, project.name, min.cells, metadata) {
  sce <- SingleCellExperiment(assays = list(counts = counts))
  cat(crayon::cyan('SCE object created\n'))
  colData(sce)$original.project <- project.name
  colData(sce)$total.counts <- colSums(assay(sce, 'counts'))
  colData(sce)$total.features <- colSums(assay(sce, 'counts') != 0)
  cat(crayon::cyan('Annotated with fundamental metadata\n'))
  if(!is.null(metadata))  {
    colData(sce) <- cbind(colData(sce), metadata)
    cat(crayon::cyan('Custom metadata appended\n'))
  }
  sce <- sce[!duplicated(rownames(assay(sce, 'counts'))), ]
  cat(crayon::cyan('Completed\n'))
  return(sce)
}

DLBCL1 <- create_sce_object(as.matrix(DLBCL1@assays$RNA@counts), project.name = 'DLBCL1', metadata = NULL)
R3104 <- create_sce_object(s$decontXcounts, project.name = 'DLBCL2', metadata = NULL)
genes <- intersect(rownames(DLBCL1), rownames(DLBCL2))
DLBCL.tmp <- cbind(DLBCL1[genes,], DLBCL2[genes,])

# scrublet

perform.scrublet <- function(object, split.by, total_counts = NULL, sim_doublet_ratio = 2.0, n_neighbors = NULL, expected_doublet_rate = 0.075, stdev_doublet_rate = 0.02, random_state = 0L) {
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
      res1 <- scrub1$scrub_doublets(min_counts = 1, min_cells = 1, min_gene_variability_pctl = 85, verbose = TRUE)
      cat(crayon::cyan('doublets detected\n'))
      raw_counts <- t(as.data.frame(raw_counts))
      raw_counts <- raw_counts[,!reticulate::py_to_r(res1)[[2]] == TRUE]
      cat(crayon::cyan('matrix scrubbed\n'))
      raw_counts_list[[l]] <- raw_counts
    }
    raw_counts <- do.call('cbind', raw_counts_list)
    storage.mode(raw_counts) <- 'integer'
    object <- object[,colnames(raw_counts)]
    assay(object, 'counts') <- raw_counts
    return(object)
    rm(obj, raw_counts, scrub1, res1, scrubbed)
  }
}

DLBCL.tmp <- perform.scrublet(object = DLBCL.tmp, split.by = 'original.project', expected_doublet_rate = 0.025)

library(celda)

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
  cat(crayon::cyan(paste0(formatC(sum(d$decontX_contamination)/length(d$decontX_contamination), digits = 2), '% average contamination\n')))
  clean.matrix <- assay(d, 'decontXcounts')
  cat(crayon::cyan('Matrix isolated\n'))
  mode(clean.matrix) <- 'integer'
  zero.samples <- colSums(clean.matrix) > 0
  object <- object[,zero.samples]
  clean.matrix <- clean.matrix[,zero.samples]
  cat(crayon::cyan('converted to integer\n'))
  assay(object, 'decontXcounts') <- clean.matrix
  cat(crayon::cyan('Added matrix\n'))
  cat(crayon::cyan('Finished\n'))
  return(object)
}

DLBCL <- perform.decontX(object = DLBCL, batch = 'original.project')

find_percentage_genes <- function(object, pattern='^MT', column.name) {
  temp <- data.frame(tmp=colSums(assay(object, 'counts')[grep(as.character(pattern),rownames(assay(object, 'counts'))),]) / colData(object)$total.counts * 100)
  colnames(temp) <- c(column.name)
  colData(object) <- cbind(colData(object), temp)
  return(object)
}

s <- find_percentage_genes(object = s, pattern = '^mt-', column.name = 'percent.mt')

filter_features_cells <- function(object, use.assay, max.features, min.features, max.mt, min.cell) {
  rowData(object)$total.cells <- rowSums(assay(object, use.assay) != 0)
  cat(crayon::cyan('Empty features removed\n'))
  object <- object[rownames(subset(rowData(object), total.cells > as.numeric(min.cell))),]
  cat(crayon::cyan('Minimum cells removed\n'))
  object <- object[,rownames(subset(colData(object), total.counts > as.numeric(min.features)))]
  cat(crayon::cyan('Minimum features removed\n'))
  object <- object[,rownames(subset(colData(object), total.features < as.numeric(max.features)))]
  cat(crayon::cyan('Maximum features removed\n'))
  object <- object[,rownames(subset(colData(object), percent.mt < as.numeric(max.mt)))]
  cat(crayon::cyan('Maximum Mitochondrial feature percentage removed\n'))
  cat(crayon::cyan('Complete!\n'))
  return(object)
}

sd.value <- sd(s$total.features)
med.value <- median(R3104$total.features)
max.features <- (sd.value*3)+med.value

s <- filter_features_cells(object = s, use.assay = 'decontXcounts', max.features = max.features, min.features = 200, max.mt = 10, min.cell = 3)



library(Seurat)
library(sctransform)

#processed.scrublet.decontx.mt.filter$original.project <- processed.scrublet.decontx.mt.filter$tech
#processed.scrublet.decontx.mt.filter <- processed.scrublet.decontx.mt.filter[,processed.scrublet.decontx.mt.filter$tech != 'smartseq2']

perform.sct.normalisation <- function(object, split.by, which.assay, new.assay = 'sctransform', ...) {
  if(is.null(split.by)) {
    seuratobj <- CreateSeuratObject(counts = assay(object, which.assay), project = 'NA')
    seuratobj <- SCTransform(object = seuratobj, do.scale = FALSE, do.center = FALSE, return.only.var.genes = FALSE)
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
      seuratobj <- SCTransform(object = seuratobj, do.scale = FALSE, do.center = FALSE, return.only.var.genes = FALSE)
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

DLBCL2 <- perform.sct.normalisation(object = DLBCL2, which.assay = 'counts', split.by = NULL)

library(scran)
library(scater)

perform.scran.normalisation <- function(object, new.assay = 'scran', max.cluster.size = 1000, scaling=NULL,
                                        do.log=TRUE, center_size_factors=TRUE) {
  
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

perform.scran.normalisation <- function(object, new.assay = 'scran', max.cluster.size = 1000, scaling=NULL,
                                        do.log=TRUE, center_size_factors=TRUE) {
  
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

zheng <- perform.scran.normalisation(object = zheng, )

perform.tpm <- function(object, which.assay, new.assay = 'tpm', log.transform) {
  
  r <- read.csv('/Users/knight05/Results/scRNA-seq/IBRAP/IBRAP/mart_export.csv', header = TRUE, sep = ',')
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
                                                        which.assay = 'decontXcounts', log.transform = TRUE)

perform.seurat.hvg <- function(object, assay, nfeatures=1500, feat.to.omit=NULL, ...) {
  if(typeof(object) != 'S4') {
    cat(crayon::cyan('Must be an S4 SCE object\n'))
    return(NULL)
  } else {
    tmp <- suppressWarnings(Seurat::as.Seurat(object, counts = NULL, data = assay))
    cat(crayon::cyan('SCE converted to seurat object\n'))
    f <- Seurat::FindVariableFeatures(object = tmp[!(rownames(tmp) %in% feat.to.omit)], nfeatures = nfeatures, ...)
    tmp <- f@assays$RNA@var.features
    cat(crayon::cyan('Variable features identified\n'))
    metadata(object)$HVGs <- tmp
    cat(crayon::cyan('HVGs added to object\n'))
  }
  return(object)
}

DLBCL2 <- perform.seurat.hvg(object = DLBCL2, assay = 'sctransform', nfeatures = 1500)


perform.scran.hvg <- function(object, assay.name, nfeatures = 1500, method = 'fish', unwanted.variance = NULL) {
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

processed.scrublet.decontx.mt.filter.tpm.scHVGs <- perform.scran.hvg(object = processed.scrublet.decontx.mt.filter.tpm, assay.name = 'tpm')

perform.seurat.scale <- function(object, use.assay, unwanted.variance, scale, centre, ...) {
  y <- object
  tmp <- suppressWarnings(Seurat::as.Seurat(y[metadata(y)$HVGs,], counts = NULL, data = use.assay))
  cat(crayon::cyan('SCE converted to seurat object\n'))
  tmp <- Seurat::ScaleData(tmp, vars.to.regress = unwanted.variance, do.scale = scale, do.center = scale, ...)
  cat(crayon::cyan('Data scaled\n'))
  reducedDim(y, 'scaled') <- t(as.matrix(tmp@assays$RNA@scale.data))
  cat(crayon::cyan('Scaled matrix attached to object\n'))
  return(y)
}

DLBCL2 <- perform.seurat.scale(object = DLBCL2, use.assay = 'sctransform', unwanted.variance = c('percent.mt'), scale = TRUE, centre = TRUE)

perform.pca <- function(object, scaled_data, use.assay, reduction.save='pca', ...) {
  if(scaled_data == TRUE) {
    cat(crayon::cyan('Initialising PCA\n'))
    a <- PCAtools::pca(mat = t(reducedDim(object, 'scaled')), center = FALSE, scale = FALSE, ...)
    reducedDim(object, 'pca') <- as.matrix(a$rotated[,1:50])
    cat(crayon::cyan('PCA completed\n'))
    return(object)
  } else if(is.na(use.assay) == TRUE) {
    cat(crayon::cyan('Please specify an assay\n'))
    return(NULL) 
  } else if(is.na(use.assay) == FALSE & scaled_data == FALSE) {
    cat(crayon::cyan('Initialising PCA\n'))
    a <- PCAtools::pca(mat = assay(object, use.assay), center = FALSE, scale = FALSE, ...)
    reducedDim(object, reduction.save) <- as.matrix(a$rotated[,1:50])
    cat(crayon::cyan('PCA completed\n'))
    return(object)
  }
}

DLBCL2 <- perform.pca(object = DLBCL2, scaled_data = TRUE, use.assay = NULL)

plot.red.sd <- function(object, reduction, n.dim, cex.names = 0.6) {
  gg <- apply(X = reducedDim(object, as.character(reduction)), MARGIN = 2, FUN = sd)
  barplot(gg, cex.names = cex.names, las = 2, main = paste0(reduction, '_var'))
}

plot.red.sd(object = DLBCL2, reduction = 'pca', n.dim = 1:25)

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
      colnames(transposed) <- column.names[[x]]
      rownames(transposed) <- integrated.corrected.data[[3]]
      dims[[x]] <- transposed
    }
    cat(crayon::cyan('Combining samples\n'))
    combined <- do.call('cbind', dims)
    cat(crayon::cyan('Samples concatenated\n'))
    assay(object, reduction.save) <- t(combined)
    cat(crayon::cyan('Scanorama completed\n'))
    return(object)
  }
  
}

DLBCL <- perform.scanorama(object = DLBCL, assay = t(reducedDim(DLBCL, 'scaled')), split.by = 'original.project', reduced.dim = TRUE, n.dims = 50, reduction.save = 'scanorama')

perform.bbknn <- function(object, reduced.dim, dims, split.by, dim.save = 'bbknn') {
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

perform.seurat.cluster <- function(object, reduction='pca', res=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5), dims=1:34,
                                   prune.SNN=0, nn.method='annoy', annoy.metric='euclidean', nn.eps=0.0, ...) {
  
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
  colData(temp) <- cbind(colData(temp), new.clusters)
  cat(crayon::cyan('Seurat clusters added\n'))
  return(temp)
  
}

DLBCL2 <- perform.seurat.cluster(object = DLBCL2, reduction = 'pca', dims = 1:23)

library(SC3)

perform.sc3 <- function(object, reduction, dims, ks, n.core=3) {
  cat(crayon::cyan('Initialising SC3 clustering\n'))
  temp.2 <- object
  temp.2 <- SingleCellExperiment(list('logcounts' = t(reducedDim(object, reduction))[dims,]))
  #logcounts(temp.2) <- assay(temp.2, 'data')
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
  metadata(object)$clustering$sc3 <- as.data.frame(new.clusters)
  cat(crayon::cyan('SC3 clustering completed\n'))
  return(object)
}

processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust.sc3clust <- perform.sc3(object = processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust, reduction = 'scanorama', dims = 1:50, ks = 11:15, n.core = 3)

perform.umap <- function(object, reduction='pca', reduction.save, n.dim, n_components = 3, ...) {
  c <- uwot::umap(X = reducedDim(object, reduction)[,n.dim], n_components = n_components, verbose = TRUE, ...)
  dim.names <- list()
  for(l in 1:n_components) {
    dim.names[[l]] <- paste0('umap_', l)
  }
  colnames(c) <- unlist(dim.names)
  reducedDim(object, reduction.save) <- c
  return(object)
}

DLBCL2 <- perform.umap(object = DLBCL2, reduction.save = 'umap', reduction = 'pca', n.dim = 1:23)

perform.tsne <- function(object, n.dims, reduction.save, n_components = 3, reduction='pca', ...) {
  if(isS4(object) == FALSE) {
    cat(crayon::cyan('Must be S4 class object\n'))
  }
  if(is.null(reducedDim(object, reduction))) {
    cat(crayon::cyan(paste0(reduction, ' does not exist\n')))
  }
  if(max(n.dims) > length(colnames(reducedDim(object, reduction)))) {
    cat(crayon::cyan('Reduction dimensions exceed actual amount available: ', length(colanmes(reducedDim(object, reduction))), '\n'))
  }
  cat(crayon::cyan('t-SNE reduction initialising\n'))
  c <- ProjectionBasedClustering::tSNE(DataOrDistances = reducedDim(object, reduction)[,n.dims], 
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

DLBCL <- perform.tsne(object = DLBCL, reduction.save = 'scanorama_tsne', n.dims = 1:50, reduction = 'harmony')

perform.dbmap <- function(object, reduction, n.dims, n_components = 3, dim.save='dbmap') {
  scipy.sparse <- reticulate::import('scipy.sparse', convert = FALSE)
  print('.')
  dbmap <- reticulate::import('dbmap', convert = FALSE)
  print('.')
  reduced.dim <- reducedDim(object, reduction)[,n.dims]
  print('.')
  cellnames <- rownames(reduced.dim)
  print('.')
  data <- scipy.sparse$csr_matrix(reticulate::r_to_py(reduced.dim))
  print('.')
  diff <- dbmap$diffusion$Diffusor(n_components = as.integer(n_components))$fit(data)
  print('.')
  res <- diff$transform(data)
  print('.')
  dbmap.embeddings <- reticulate::py_to_r(res)
  print('.')
  names(dbmap.embeddings) <- cellnames
  print('.')
  dim.names <- list()
  for(t in 1:n_components) {
    dim.names[[t]] <- paste0('dbmap_', t)
  }
  colnames(dbmap.embeddings) <- unlist(dim.names)
  reducedDim(object, dim.save) <- dbmap.embeddings
  return(object)
}

processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust.sc3clust.umap.tsne.dbmap <- perform.dbmap(object = processed.scrublet.decontx.mt.filter.tpm.hvgs.scale.pca.scanorama.seuratclust.sc3clust.umap.tsne.dbmap, reduction = 'pca', n.dims = 1:50)

library(cluster)
library(clValid)
library(mclust)

benchmarking_clustering <- function(object, components, reduction, dist.method='euclidean',
                                    cluster.columns, ARI=FALSE, ground.truth) {
  all.clusters <- colData(object)[,cluster.columns]
  print(all.clusters)
  dims <- reducedDim(object, reduction)[,components]
  dist.matrix <- dist(x = dims, method = dist.method)
  #all.clusters <- all.clusters[, colSums(all.clusters != 0) > 0]
  #all.clusters <- all.clusters[!is.na(names(all.clusters))]
  print(all.clusters)
  sil.results <- data.frame(average_silhoeutte=NA)
  for (v in colnames(all.clusters)[1:length(colnames(all.clusters))]) {
    print(paste0('Calculating silhouette for ', v))
    tmp <- silhouette(x = as.numeric(x = as.factor(x = all.clusters[,v])), dist = dist.matrix)
    average <- sum(tmp[,3])/length(tmp[,3])
    sil.results[v,] <- average
  }
  sil.results <- sil.results[complete.cases(sil.results),]
  max.AS <- max(sil.results)
  print(max.AS)
    
  dunn.results <- data.frame(dunn.index=NA)
  for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
    print(paste0('Calculating dunn index for ', p))
    dunn.results[p,] <- dunn(distance = dist.matrix, clusters = as.numeric(x = as.factor(x = all.clusters[,p])))
  }
  dunn.results <- dunn.results[complete.cases(dunn.results),]
  max.dunn <- max(dunn.results)
  
  conn.results <- data.frame(connectivity=NA)
  for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
    print(paste0('Calculating connectivity for ', p))
    conn.results[p,] <- connectivity(distance = dist.matrix, clusters = all.clusters[,p])
  }
  conn.results <- conn.results[complete.cases(conn.results),]
  max.conn <- max(conn.results)
  
  if(ARI == TRUE) {
    ARI.results <- data.frame(ARI=NA)
    for (p in colnames(all.clusters)[1:length(colnames(all.clusters))]){
      print(paste0('Calculating ARI for ', p))
      ARI.results[p,] <- adjustedRandIndex(x = all.clusters[,p], y = ground.truth)
    }
    ARI.results <- ARI.results[complete.cases(ARI.results),]
    max.ARI <- max(ARI.results)
    results <- cbind(sil.results, dunn.results, conn.results, ARI.results)
    rownames(results) <- colnames(all.clusters)
    colnames(results) <- c(paste0(g, '_sil.results'), paste0(g, '_dunn.results'), paste0(g, '_conn.results'), paste0(g, '_ARI.results'))
    metadata(object)[['benchmarking_clustering']][[as.character(g)]] <- results
  } else {
    results <- cbind(sil.results, dunn.results, conn.results)
    rownames(results) <- colnames(all.clusters)
    colnames(results) <- c(paste0(g, '_sil.results'), paste0(g, '_dunn.results'), paste0(g, '_conn.results'))
    metadata(object)[['benchmarking_clustering']][[as.character(g)]] <- results
  }
  return(object)
}

DLBCL <- benchmarking_clustering(object = DLBCL, components = 1:3, reduction = 'umap', dist.method = 'euclidean', cluster.columns = c('Seurat_res_0.3', 'Seurat_res_0.4', 'Seurat_res_0.5', 'Seurat_res_0.6', 'Seurat_res_0.7', 'Seurat_res_0.8', 'Seurat_res_0.9', 'Seurat_res_1'), ARI = FALSE, ground.truth = NULL)

library(plotly)
library(ggplot2)

plot_cluster_dr <- function(object, reduction='', pt.size=0.5, group.by, interactive=FALSE) {
  print('plot_cluster_dr_started')
  iso <- colData(object)[,group.by]
  print('.')
  reduction <- as.data.frame(reducedDim(object, as.character(reduction)))
  print('.')
  meta <- colData(object)
  print('.')
  meta.colnames <- colnames(meta)
  print('.')
  barcodes <- colnames(object)
  print('.')
  x.val <- reduction[,1]
  print('.')
  y.val <- reduction[,2]
  print('.')
  results <- cbind(barcodes, meta, iso, x.val, y.val)
  print('.')
  colnames(results) <- c('barcode', meta.colnames, 'key', 'x_value', 'y_value')
  print('.')
  results <- as.data.frame(results)
  print(results)
  plot <- ggplot(data = as.data.frame(results), aes(ident = `barcode`, x = `x_value`,
                                                                      y = `y_value`, color = `key`)) +
    geom_point(shape = 16, size = pt.size) +
    scale_color_hue() +
    guides(fill=FALSE, alpha=FALSE, size=FALSE) +
    theme_minimal() +
    #xlab(paste0(reduction, '_1')) + ylab(paste0(reduction, '_2')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                   plot.title = element_text(hjust = 0.5))
  if(interactive == TRUE) {
    z <- plotly::ggplotly(plot)
    return(z)
  } else {
    return(plot)
  }
  print('plot_cluster_dr_completed')
}
plot_cluster_dr(object = DLBCL2, reduction = 'umap', group.by = 'Seurat_res_0.5', interactive = TRUE)





