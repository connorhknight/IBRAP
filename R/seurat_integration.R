#' @name perform.seurat.cca
#' @aliases perform.seurat.cca
#' 
#' @title Performs Seurat Integration
#'
#' @description Performs Seurat integration on the supplied assay names. Results are saved under integration_reductions 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param batch Character. Which column in the metadata defines the batches
#' @param nfeatures Numerical. How many features should be found as integration anchors. Default = 3000
#' @param l2.norm Logical. Perform L2 normalization on the CCA cell embeddings after dimensional reduction. Default = TRUE
#' @param k.anchor Numerical. How many neighbors (k) to use when picking anchors. Default = 5
#' @param k.filter Numerical. How many neighbors (k) to use when filtering anchors. Default = 200
#' @param k.score Numerical. How many neighbors (k) to use when scoring anchors. Default = 30
#' @param max.features. Numerical. The maximum number of features to use when specifying the neighborhood search space in the anchor filtering. Default = 200
#' @param nn.method Character. Method for nearest neighbor finding. Options include: rann, annoy. Default = annoy
#' @param n.tree Numerical. More trees gives higher precision when using annoy approximate nearest neighbor search. Default = 50
#' @param anchor.eps Numerical. Error bound on the neighbor finding algorithm (from RANN/Annoy) when finding integration genes. 
#' @param features Character. Vector of features to use when computing the PCA to determine the weights. Only set if you want a different set from those used in the anchor finding process. Default = NULL
#' @param integrate.dims Numerical. Number of dimensions to use in the anchor weighting procedure. Default = 1:30
#' @param k.weight Numerical. Number of neighbors to consider when weighting anchors. Default = 100
#' @param sd.weight Numerical. Controls the bandwidth of the Gaussian kernel for weighting. Default = 1
#' @param sample.tree Character. Specify the order of integration. If NULL, will compute automatically. Default = NULL
#' @param integrate.eps Numerical. Error bound on the neighbor finding algorithm (from RANN)
#' @param save.plot Boolean. Should the automatically genewrated plot be saved? Default = TRUE
#' 
#' @return Produces a new 'methods' assay containing normalised, scaled and HVGs.
#' 
#' @examples perform.seurat.cca <- function(object = object, 
#'                                          assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                                          normalisation.method = c('perform.sct', 'perform.scran', 'perform.scanpy'), 
#'                                          batch = original.project, 
#'                                          reduction.name.suffix=NULL,
#'                                          nfeatures = 3000,
#'                                          reduction = 'cca',
#'                                          anchors.dims = 1:30,
#'                                          l2.norm = T,
#'                                          k.anchor = 5,
#'                                          k.filter = 200,
#'                                          k.score = 30, 
#'                                          max.features = 200, 
#'                                          nn.method = 'annoy', 
#'                                          n.trees = 50, 
#'                                          anchor.eps = 0, 
#'                                          features = NULL,
#'                                          integrate.dims = 1:30,
#'                                          k.weight = 100,
#'                                          sd.weight = 1, 
#'                                          sample.tree = NULL, 
#'                                          integrate.eps = 0)
#'
#' @export

perform.seurat.integration <- function(object, 
                                       assay, 
                                       normalisation.method, 
                                       batch, 
                                       reduction.name.suffix=NULL,
                                       nfeatures = 2000,
                                       anchors.dims = 1:30,
                                       l2.norm = T,
                                       k.anchor = 5,
                                       k.filter = 200,
                                       k.score = 30, 
                                       max.features = 200, 
                                       nn.method = 'annoy', 
                                       n.trees = 50, 
                                       anchor.eps = 0, 
                                       features = NULL,
                                       integrate.dims = 1:30,
                                       k.weight = 100,
                                       sd.weight = 1, 
                                       sample.tree = NULL, 
                                       integrate.eps = 0,
                                       print.variance = TRUE,
                                       verbose=FALSE,
                                       seed = 1234,
                                       ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string(s) \n')
    
  } else if(is.character(assay)) {
    
    for(x in assay) {
      
      if(!x %in% names(object@methods)) {
        
        stop(paste0(x, ' is not contained within object@methods \n'))
        
      }
      
    }
    
  }
  
  if(!is.character(normalisation.method)) {
    
    stop('normalisation.method must be character string(s) \n')
    
  } else if(is.character(normalisation.method)) {
    
    for(x in normalisation.method) {
      
      if(!x %in% c('perform.sct', 'perform.scran', 'perform.scanpy', 'perform.tpm')) {
        
        if(isFALSE(x = x)) {
          
          stop(paste0(x, ' function does not exist in the global environment. Does it exist? \n'))
          
        }
        
      }
      
    }
    
  }
  
  if(!is.character(batch)) {
    
    stop('batch must be character string(s) \n')
    
  } else if (is.character(batch)) {
    
    for(x in batch) {
      
      if(!x %in% names(object@sample_metadata)) {
        
        stop(paste0(x, ' is not contained within object@sample_metadata \n'))
        
      }
      
    }
    
  }
  
  if(!is.null(reduction.name.suffix)) {
    
    if(!is.character(reduction.name.suffix)) {
      
      stop('reduction.name.suffix must be a character string \n')
      
    }
    
  }
  
  if(!is.numeric(nfeatures)) {
    
    stop('nfeatures must be cnumerical \n')
    
  }
  
  if(!is.logical(l2.norm)) {
    
    stop('l2.norm must be logical. TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(anchors.dims)) {
    
    stop('anchors.dims must be numerical \n')
    
  }
  
  if(!is.numeric(k.anchor)) {
    
    stop('k.anchor must be numerical \n')
    
  }
  
  if(!is.numeric(k.filter)) {
    
    stop('k.filter must be numerical \n')
    
  }
  
  if(!is.numeric(k.score)) {
    
    stop('k.score must be numerical \n')
    
  }
  
  if(!is.numeric(max.features)) {
    
    stop('max.features must be numerical \n')
    
  }
  
  if(!is.character(nn.method)) {
    
    stop('nn.method must be numerical \n')
    
  }
  
  if(!is.numeric(anchor.eps)) {
    
    stop('anchor.eps must be numerical \n')
    
  }
  
  if(!is.null(features)) {
    
    if(!is.character(features)) {
      
      stop('features must be character string(s)')
      
    } 
    
  }
  
  if(!is.numeric(integrate.dims)) {
    
    stop('integrate.dims must be numerical \n')
    
  }
  
  if(!is.numeric(k.weight)) {
    
    stop('k.weight must be numerical \n')
    
  }
  
  if(!is.numeric(sd.weight)) {
    
    stop('sd.weight must be numerical \n')
    
  }
  
  if(!is.null(sample.tree)) {
    
    if(!is.character(sample.tree)) {
      
      stop('must be character stirng \n')
      
    }
    
  }
  
  if(!is.numeric(integrate.eps)) {
    
    stop('integrate.eps must be numerical. TRUE/FALSE \n')
    
  }
  
  if(!is.logical(print.variance)) {
    
    stop('print.variance must be boolean. TRUE/FALSE \n')
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose must be boolean. TRUE/FALSE \n')
    
  }
  
  if(length(batch) > 1) {
    
    temp <- function(x) {
      
      return(paste(x, collapse = '_'))
      
    }
    
    df <- object@sample_metadata[,batch]
    object2 <- object
    object2@sample_metadata$batch <- apply(X = df, MARGIN = 1, FUN = temp)
    
    batch <- 'batch'
    
  } else {
    
    object2 <- object
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed must be numerical \n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  if(!'integration_method' %in% colnames(object@pipelines)) {
    
    tmp <- tibble::add_column(.data = object@pipelines, integration_method=NA, integration_time=NA)
    
  } else {
    
    tmp <- object@pipelines
    
  }
  
  count <- 1
  
  for(a in assay) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': initialising seurat integration for assay: ', a, '\n')))
      
    }
    
    if(normalisation.method[[count]] != 'perform.sct') {
      
      if(normalisation.method[[count]] == 'perform.scran') {
        
        normalise <- function(...) {
          
          perform.scran(...)
          
        }
        
      } else if (normalisation.method[[count]] == 'perform.scanpy') {
        
        normalise <- function(...) {
          
          perform.scanpy(...)
          
        }
        
      } else if (normalisation.method[[count]] == 'perform.tpm') {
        
        normalise <- function(...) {
          
          IBRAP::perform.tpm(...)
          
        }      
        
      } else {
        
        stop('normalisation.method must be either: perform.sct, perform.scran, perfrom.scanpy or perform.tpm \n')
        
      }
      
      obj.list <- list()
      
      for(x in unique(object2@sample_metadata[,batch])) {
        print('.')
        obj.list[[x]] <- suppressMessages(IBRAP::createIBRAPobject(counts = object@methods$RAW@counts[,object2@sample_metadata[,batch]==x], 
                                                                   original.project = x, add.suffix = F, verbose = F, 
                                                                   min.cells = 0, min.features = 0))

      }

      obj.list <- lapply(X = obj.list, FUN = 'normalise')

      new.list <- list()
      
      for(x in names(obj.list)) {
        
        new.list[[x]] <- suppressWarnings(Seurat::CreateSeuratObject(counts = obj.list[[x]]@methods[[2]]@normalised, project = x))
        
        new.list[[x]]@assays$RNA@var.features <- obj.list[[x]]@methods[[2]]@highly.variable.genes
        
      }
      
      features.list <- suppressWarnings(Seurat::SelectIntegrationFeatures(object.list = new.list, nfeatures = nfeatures, verbose = verbose, ... = ))
      
      anchors <- suppressWarnings(Seurat::FindIntegrationAnchors(object.list = new.list, 
                                                                 anchor.features = features.list, 
                                                                 reduction = reduction,
                                                                 scale = T, 
                                                                 l2.norm = l2.norm, 
                                                                 dims = anchors.dims, 
                                                                 k.anchor = k.anchor, 
                                                                 k.filter = k.filter, 
                                                                 k.score = k.score, 
                                                                 max.features = max.features,
                                                                 nn.method = nn.method, 
                                                                 n.trees = n.trees, 
                                                                 eps = anchor.eps, 
                                                                 verbose = verbose))

      to_integrate <- Reduce(intersect, lapply(anchors@object.list, rownames))
      
      combined <- suppressWarnings(Seurat::IntegrateData(anchorset = anchors, 
                                                         features = features,
                                                         features.to.integrate = to_integrate, 
                                                         dims = integrate.dims,
                                                         k.weight = k.weight,
                                                         sd.weight = sd.weight, 
                                                         sample.tree = sample.tree, 
                                                         preserve.order = T, 
                                                         eps = integrate.eps, 
                                                         verbose = verbose))
      
      Seurat::DefaultAssay(combined) <- 'integrated'
      
      combined <- suppressWarnings(Seurat::ScaleData(object = combined))
      
      tmp.obj <- IBRAP::createIBRAPobject(counts = combined@assays$integrated@data, 
                                          original.project = 'tmp', 
                                          add.suffix = F)
      
      tmp.obj@methods[['RAW']]@norm.scaled <- as.matrix(combined@assays$integrated@scale.data)
      
      tmp.obj@methods[['RAW']]@highly.variable.genes <- combined@assays$integrated@var.features
      
      tmp.obj <- IBRAP::perform.pca(object = tmp.obj, assay = 'RAW', slot = 'norm.scaled', print.variance = T)
      
      if(is.null(reduction.name.suffix)) {
        
        object@methods[[a]]@integration_reductions[[paste0('seurat', reduction.name.suffix[[count]])]] <- tmp.obj@methods[[1]]@computational_reductions[[1]]
        
      } else if(!is.null(reduction.name.suffix)) {
        
        if(!is.character(reduction.name.suffix)) {
          
          stop('reduction.name.suffix must be character string \n')
          
        } else {
          
          if('_' %in% unlist(strsplit(x = reduction.name.suffix, split = ''))) {
            
            if(isTRUE(verbose)) {
              
              cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in reduction.name.suffix, replacing with - \n')))
              
            }
            
            reduction.name.suffix <- sub(pattern = '_', replacement = '-', x = reduction.name.suffix)
            
          }
          
          object@methods[[a]]@integration_reductions[[paste0('seurat', reduction.name.suffix[[count]])]] <- tmp.obj@methods[[1]]@computational_reductions[[1]]
          
        }
        
      }
      
    } else if(normalisation.method[[count]] == 'perform.sct') {
      
      normalise <- function(...) {
        
        Seurat::SCTransform(...)
        
      }
      
      obj.list <- list()
      
      for(x in unique(object2@sample_metadata[,batch])) {
        
        obj.list[[x]] <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods$RAW@counts[,object2@sample_metadata[,batch]==x]))
        
      }
      
      obj.list <- lapply(X = obj.list, FUN = 'normalise')
      
      features <- suppressWarnings(Seurat::SelectIntegrationFeatures(object.list = obj.list, nfeatures = nfeatures, verbose = verbose))
      
      obj.list <- suppressWarnings(Seurat::PrepSCTIntegration(object.list = obj.list,anchor.features = features, verbose = verbose))
      
      anchors <- suppressWarnings(Seurat::FindIntegrationAnchors(object.list = obj.list, 
                                                                 normalization.method = "SCT", 
                                                                 anchor.features = features, 
                                                                 reduction = 'cca',
                                                                 scale = F, 
                                                                 l2.norm = l2.norm, 
                                                                 dims = anchors.dims, 
                                                                 k.anchor = k.anchor, 
                                                                 k.filter = k.filter, 
                                                                 k.score = k.score, 
                                                                 max.features = max.features,
                                                                 nn.method = nn.method, 
                                                                 n.trees = n.trees, 
                                                                 eps = anchor.eps, 
                                                                 verbose = verbose))
      
      to_integrate <- Reduce(intersect, lapply(anchors@object.list, rownames))
      
      combined <- suppressWarnings(Seurat::IntegrateData(anchorset = anchors, 
                                                         normalization.method = "SCT", 
                                                         features = features,
                                                         features.to.integrate = to_integrate, 
                                                         dims = integrate.dims,
                                                         k.weight = k.weight,
                                                         sd.weight = sd.weight, 
                                                         sample.tree = sample.tree, 
                                                         preserve.order = T, 
                                                         eps = integrate.eps, 
                                                         verbose = verbose))
      
      tmp.obj <- IBRAP::createIBRAPobject(counts = combined@assays$integrated@data, 
                                          original.project = 'tmp', 
                                          add.suffix = F)
      
      tmp.obj@methods[['RAW']]@norm.scaled <- as.matrix(combined@assays$integrated@scale.data)
      
      tmp.obj@methods[['RAW']]@highly.variable.genes <- combined@assays$integrated@var.features
      
      tmp.obj <- IBRAP::perform.pca(object = tmp.obj, assay = 'RAW', slot = 'norm.scaled', print.variance = print.variance)
      
      if(is.null(reduction.name.suffix)) {
        
        object@methods[[a]]@integration_reductions[[paste0('seurat', reduction.name.suffix[[count]])]] <- tmp.obj@methods[[1]]@computational_reductions[[1]]
        
      } else if(!is.null(reduction.name.suffix)) {
        
        if(!is.character(reduction.name.suffix)) {
          
          stop('reduction.name.suffix must be character string \n')
          
        } else {
          
          if('_' %in% unlist(strsplit(x = reduction.name.suffix, split = ''))) {
            
            if(isTRUE(verbose)) {
              
              cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in reduction.name.suffix, replacing with - \n')))
              
            }

            reduction.name.suffix <- sub(pattern = '_', replacement = '-', x = reduction.name.suffix)
            
          }
          
          object@methods[[a]]@integration_reductions[[paste0('seurat', reduction.name.suffix[[count]])]] <- tmp.obj@methods[[1]]@computational_reductions[[1]]
          
        }
        
      }
      
    }
    
    count <- count + 1
    
    end_time <- Sys.time()
    
    function_time <- end_time - start_time
    
    if(!'integration_method' %in% colnames(object@pipelines)) {
      
      tmp[which(x = tmp$normalisation_method==p),'integration_method'] <- paste0('CCA', reduction.save.suffix)
      
      tmp[which(x = tmp$normalisation_method==p),'integration_time'] <- as.difftime(function_time, units = 'secs')
      
    }
    
    if('integration_method' %in% colnames(object@pipelines)) {
      
      if(paste0('CCA', reduction.save.suffix) %in% tmp$integration_method) {
        
        tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('CCA', reduction.save.suffix)),] <- c(tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('CCA', reduction.save.suffix)),c('normalisation_method','normalisation_time')], paste0('CCA', reduction.save.suffix), as.difftime(function_time, units = 'secs'))  
        
      }
      
      if(!paste0('CCA', reduction.save.suffix) %in% object@pipelines$integration_method) {
        
        df <- tmp[which(tmp$normalisation_method==p),]
        
        df <- df[!duplicated(df$normalisation_method),]
        
        df[,'integration_method'] <- paste0('CCA', reduction.save.suffix)
        
        df[,'integration_time'] <- function_time
        
        tmp <- rbind(tmp, df)
        
      }
      
    }
    
  }
  
  if(!'integration_method' %in% colnames(object@pipelines)) {
    
    tmp$integration_time <- as.difftime(tim = tmp$integration_time, units = 'secs')
    
    rownames(tmp) <- 1:nrow(tmp)
    
    object@pipelines <- tmp
    
  } else if ('integration_method' %in% colnames(object@pipelines)) {
    
    tmp$integration_time <- as.difftime(tim = tmp$integration_time, units = 'secs')
    
    rownames(tmp) <- 1:nrow(tmp)
    
    object@pipelines <- tmp
    
  }
  
  return(object)
  
}


