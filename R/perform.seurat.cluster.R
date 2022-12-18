#' @name perform.graph.cluster
#' @aliases perform.graph.cluster
#' 
#' @title Performs graph-based clustering
#'
#' @description Performs graph-based clustering on previously generated neighbouhood graphs. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param neighbours Character. String indicating which neighbourhood graphs should be used. 
#' @param algorithm Numerical. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default = 1 Default = NULL
#' @param cluster.df.name Character. What to call the df contained in clusters. Default = 'seurat
#' @param res Numerical vector. Which resolution to run the clusterign algorithm at, a smaller and larger value identified less and more clusters, respectively. Default = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
#' @param verbose Logical. Should system information be printed. Default = FALSE
#' @param seed Numeric. What should the seed be set as. Default = 1234
#' @param ... arguments to be passed to Seurat::FindClusters
#' 
#' @return Cluster assignments using the list of resolutions provided contained within cluster_assignments under cluster.df.name
#'
#' @examples 
#' 
#' object <- perform.nn.v1(object = object, assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                         reduction = c('pca_harmony','scanorama'), 
#'                         dims = list(0,0), generate.diffmap = T)
#' 
#' object <- perform.nn.v1(object = object, assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                         reduction = c('pca_bbknn_bbknn:diffmap','pca_harmony_nn.v1:diffmap', 'scanorama_nn.v1:diffmap'), 
#'                         dims = list(0,0,0))
#' 
#' object <- perform.nn.v2(object = object, assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                        reduction = c('pca_harmony','scanorama','pca_bbknn_bbknn:diffmap',
#'                                      'pca_harmony_nn.v1:diffmap', 'scanorama_nn.v1:diffmap'), 
#'                        dims = list(0,0,0,0,0))
#'                        
#' object <- perform.graph.cluster(object = object, assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                                 neighbours = c("pca_bbknn_bbknn",
#'                                                "pca_harmony_nn.v1",
#'                                                "scanorama_nn.v1",
#'                                                "pca_bbknn_bbknn:diffmap_nn.v1",
#'                                                "pca_harmony_nn.v1:diffmap_nn.v1",
#'                                                "scanorama_nn.v1:diffmap_nn.v1",
#'                                                "pca_harmony_nn.v2",
#'                                                "scanorama_nn.v2",
#'                                                "pca_bbknn_bbknn:diffmap_nn.v2",
#'                                                "pca_harmony_nn.v1:diffmap_nn.v2",
#'                                                "scanorama_nn.v1:diffmap_nn.v2" ), 
#'                                 algorithm = 1)
#'
#' @export

perform.graph.cluster <- function(object, 
                                  assay,
                                  neighbours, 
                                  algorithm=1,
                                  cluster.df.name.suffix = '',
                                  res=c(0.1,0.2,0.3,0.4,0.5,
                                        0.6,0.7,0.8,0.9,1,
                                        1.1,1.2,1.3,1.4,1.5), 
                                  verbose=FALSE,
                                  seed=1234,
                                  ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string(s)')
    
  } else if (is.character(assay)) {
    
    for(x in assay) {
      
      if(!x %in% names(object@methods)) {
        
        stop(paste0('assay: ', x, ' is not contained within object@methods \n'))
        
      }
      
    }
    
  }
  
  for(p in assay) {
    
    if(!is.character(neighbours)) {
      
      stop('neighbours must be character string(s)')
      
    } else if (is.character(neighbours)) {
      
      for(x in neighbours) {
        
        if(!x %in% names(object@methods[[p]]@neighbours)) {
          
          stop(paste0('neighbours: ', x, ' is not contained within assay: ', p,  ' \n'))
          
        }
        
      }
      
    }
    
  }
  
  if(!is.numeric(algorithm)) {
    
    stop('method must either be numeric\n')
    
  }
  
  if(!algorithm %in% c(1,2,3,4)) {
    
    stop('method must either: 1 (original Louvain), 2 (Louvain with multilevel refinement), 3 (SLM) or 4 (Leiden)\n')
    
  }
  
  if(algorithm==1) {
    
    algo.name <- 'LOUVAIN'
    
  }
  
  if(algorithm==2) {
    
    algo.name <- 'LOUVAINMLR'
    
  }
  
  if(algorithm==3) {
    
    algo.name <- 'SLM'
    
  }
  
  if(algorithm==4) {
    
    algo.name <- 'LEIDEN'
    
  }
  
  if('_' %in% unlist(strsplit(x = cluster.df.name.suffix, split = ''))) {
    
    if(isTRUE(verbose)) {
      
      cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in cluster.df.name.suffix, replacing with - \n')))
      
    }
    
    cluster.df.name.suffix <- sub(pattern = '_', replacement = '-', x = cluster.df.name.suffix)
    
  }
  
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical \n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  
  # if(!'clustering_method' %in% colnames(object@pipelines)) {
  #   
  #   tmp <- tibble::add_column(.data = object@pipelines, clustering_method=NA, clustering_time=NA)
  #   
  # } else {
  #   
  #   tmp <- object@pipelines
  #   
  # }
    
  for(p in assay) {
    
    for(t in neighbours) {
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': calculating Seurat clusters for assay: ', p, ' using graph: ', t, '\n')))
        
      }
      
      seuobj <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@methods[[p]]@counts))
      
      seuobj@graphs[['neighbourhood_graph']] <- suppressWarnings(Seurat::as.Graph(x = object@methods[[p]]@neighbours[[t]]$connectivities))
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': graph added\n')))
        
      }
      
      for(x in res) {
        
        start_time <- Sys.time()
        
        if(isTRUE(verbose)) {
          
          seuobj <- suppressWarnings(Seurat::FindClusters(object = seuobj, resolution = x, graph.name = 'neighbourhood_graph', 
                                                          algorithm = algorithm, verbose = TRUE, ...))
          
        } else if(isFALSE(verbose)) {
          
          seuobj <- suppressWarnings(Seurat::FindClusters(object = seuobj, resolution = x, graph.name = 'neighbourhood_graph', 
                                                          algorithm = algorithm, verbose = FALSE, ...))
          
        }
        
        # end_time <- Sys.time()
        # 
        # function_time <- end_time - start_time
        # 
        # if('integration_method' %in% colnames(tmp) & 'normalisation_method' %in% colnames(tmp)) {
        #   
        #   if(nrow(tmp[which(tmp$normalisation_method==p & tmp$integration_method==unique(tmp$integration_method[which(tmp$integration_method %in% unlist(strsplit(x = toupper(t), split = '_')))]) & tmp$clustering_method==paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix)),])!=0) {
        #     
        #     tmp[which(tmp$normalisation_method==p & tmp$integration_method==unique(tmp$integration_method[which(tmp$integration_method %in% unlist(strsplit(x = toupper(t), split = '_')))]) & tmp$clustering_method==paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix)),'clustering_method'] <- paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix) 
        #     
        #     tmp[which(tmp$normalisation_method==p & tmp$integration_method==unique(tmp$integration_method[which(tmp$integration_method %in% unlist(strsplit(x = toupper(t), split = '_')))]) & tmp$clustering_method==paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix)),'clustering_time'] <- as.difftime(function_time, units = 'secs')
        #     
        #   } else {
        #     
        #     if(any(tmp$integration_method %in% unlist(strsplit(x = toupper(t), split = '_'))==TRUE)){
        #       
        #       int_meth <- unique(tmp$integration_method[which(tmp$integration_method %in% unlist(strsplit(x = toupper(t), split = '_')))])
        #       
        #       roww <- which(tmp$normalisation_method==p & tmp$integration_method==int_meth)
        #       
        #       tmp_row <- tmp[which(tmp$normalisation_method==p & tmp$integration_method==int_meth),]
        #       
        #       if(NA %in% tmp_row[,'clustering_method']) {
        #         
        #         tmp <- tmp[-roww,]
        #         
        #       }
        #       
        #       tmp_row <- tmp_row[1,]
        #       
        #       tmp_row[,'clustering_method'] <- paste0(algo.name, '_resolution(', x, ')', cluster.df.name.suffix)
        #       
        #       tmp_row[,'clustering_time'] <- as.difftime(function_time, units = 'secs')
        #       
        #       tmp <- rbind(tmp, tmp_row)
        #       
        #     }
        #     
        #     if(all(tmp$integration_method %in% unlist(strsplit(x = toupper(t), split = '_'))==FALSE)) {
        #       
        #       tmp_row <- tmp[which(x = tmp$normalisation_method==p),]
        #       
        #       tmp_row <- tmp_row[1,]
        #       
        #       tmp_row[,c('integration_method','integration_time')] <- NA
        #       
        #       tmp_row[,'clustering_method'] <- paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix)
        #       
        #       tmp_row[,'clustering_time'] <- as.difftime(function_time, units = 'secs')
        #       
        #       tmp <- rbind(tmp, tmp_row)
        #       
        #     } 
        #     
        #   }
        #   
        # }  else if ('normalisation_method' %in% colnames(tmp) & !'integration_method' %in% colnames(tmp)) {
        #   
        #   tmp_sub <- tmp[tmp$normalisation_method==p,]
        #   
        #   if(paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix) %in% tmp_sub$clustering_method) {
        #     
        #     tmp[which(tmp$normalisation_method==p & tmp$clustering_method==paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix)),'clustering_method'] <- paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix) 
        #     
        #     tmp[which(tmp$normalisation_method==p & tmp$clustering_method==paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix)),'clustering_time'] <- as.difftime(function_time, units = 'secs')
        #     
        #   } else if(!paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix) %in% tmp_sub$clustering_method) {
        #     
        #     tmp_sub[1,'clustering_method'] <- paste0(algo.name, '_(resolution=', x, ')', cluster.df.name.suffix) 
        #     
        #     tmp_sub[1,'clustering_time'] <- as.difftime(function_time, units = 'secs')
        #     
        #     tmp <- rbind(tmp, tmp_sub)
        #     
        #     tmp <- tmp[!duplicated(tmp),]
        #     
        #   }
        #   
        # }
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': clusters identified\n')))
          
        }
        
      }
      
      new.clusters <- seuobj@meta.data[,colnames(seuobj@meta.data)[grepl(pattern = 'res', x = colnames(seuobj@meta.data))]]
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': seurat clusters added\n')))
        
      }
      
      if(!is.data.frame(new.clusters)) {
        
        new.clusters <- as.data.frame(new.clusters)
        
        rownames(new.clusters) <- rownames(seuobj@meta.data)
        
        colnames(new.clusters) <- paste0('res_', res)
        
      } else {
        
        colnames(new.clusters) <- paste0('res_', res)
        
      }
      
      object@methods[[p]]@cluster_assignments[[paste0(t, ':', algo.name, cluster.df.name.suffix)]] <- new.clusters
      
    }
    
  }
  
  # for(p in assay) {
  #   
  #   if(any(is.na(tmp[which(tmp$normalisation_method==p),'clustering_method']))) {
  #     
  #     tmp <- tmp[-which(tmp$normalisation_method==p & is.na(tmp$clustering_method)),]
  #     
  #   }
  #   
  # }
  # 
  # tmp$clustering_time <- as.difftime(tim = tmp$clustering_time, units = 'secs')
  # 
  # if('integration_method' %in% colnames(tmp) & 'normalisation_method' %in% colnames(tmp)) {
  #   
  #   tmp$clustering_time <- tmp$normalisation_time + tmp$integration_time + tmp$clustering_time
  #   
  # } else if ('normalisation_method' %in% colnames(tmp) & !'integration_method' %in% colnames(tmp)) {
  #   
  #   tmp$clustering_time <- tmp$normalisation_time + tmp$clustering_time
  #   
  # }
  # 
  # rownames(tmp) <- 1:nrow(tmp)
  # 
  # object@pipelines <- tmp
  # 
  return(object)
  
}
