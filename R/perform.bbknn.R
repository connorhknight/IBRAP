#' @name perform.bbknn
#' @aliases perform.bbknn
#' 
#' @title Performs BBKNN integration
#'
#' @description Performs BBKNN integration on defined method-assays and reductions contained within. This is performed on reductions. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param reduction Character. String defining the name of the reduction to provide for BBKNN. Default = NULL
#' @param graph.name.suffix Character. Should a suffix be added to the end of bbknn as the graph name, i.e. parameter changes?
#' @param batch Character. Column name in metadata indicating batch. Can be multiple.
#' @param approx Character. Employs annoy's approximate neighbour finding. Useful for large datasets but may increase correction. 
#' @param metric. Character. Which distance metric to use when approx is TRUE, options: 'angular', 'euclidean', 'manhattan' or 'hamming'. Default = 'euclidean'
#' @param neighbors_within_batch Numerical. How many neighbours to report per batch. Default = 3
#' @param n_pcs Numerical. Range of principal components to use. Default = NULL
#' @param trim Numerical. Trims the n of neighbours per cell to this value. Helps with population independence. Default = NULL
#' @param annoy_n_trees Numerical. Number of trees to generate in annoy forest. More trees provides higher precision at the cost of increased resource demand and run time. Default = 10
#' @param use_faiss Boolean. Uses faiss package to compute nearest neighbour, this improves run time at the cost of precision. Default = TRUE
#' @param set_op_mix_ratio Numerical. UMAP connectivity parameter between 0 and 1. controls the blen d between a connectivity matrix formed exclusively from mutual nearest neighbour pairs (0) and a union of all observed neighbour relationships with the mutual pairs emphasised (1). Default = 1.0
#' @param local_connectivity Numerical. How many nearest neighbours of each cell are assumed to be fully connected. Default = 1
#' @param generate.diffmap Boolean. Should diffusion maps be generated from the neighourhood graphs, these will be stored in computational_reductions and can be used for umap generation and further neighbourhood generation. Default = TRUE
#' @param n_comps Numerical. How many components should be generated for the diffusion maps. Default = 15
#' @param diffmap.name.sufix Character. Should a suffix be added to the end of bbknn:diffmap as the reduction name, i.e. parameter changes?
#' @param verbose Logical Should function messages be printed?
#' @param seed Numerical What seed should be set. Default = 1234
#' 
#' @return BBKNN connectivity graph contained in graphs in the indicated method-assays
#' 
#' @examples 
#' 
#' object <- perform.bbknn(object = object, 
#'                         assay = c('SCT', 'SCANPY', 'SCRAN'), 
#'                         reduction = c('pca'),
#'                         batch = 'tech')
#'                         
#' @export

perform.bbknn <- function(object,
                          assay,
                          reduction,
                          graph.name.suffix = '',
                          batch,
                          approx = FALSE,
                          metric = 'euclidean',
                          neighbors_within_batch = 3,
                          n_pcs = NULL,
                          trim = NULL,
                          annoy_n_trees = 10,
                          use_faiss = TRUE,
                          set_op_mix_ratio = 1.0,
                          local_connectivity= 1,
                          generate.diffmap = FALSE,
                          n_comps = 15,
                          diffmap.name.suffix='',
                          verbose = FALSE,
                          seed=1234) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP\n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('reduction: ', x, 'does not exist\n'))
      
    }
    
  }
  
  if(!is.character(reduction)) {
    
    stop('reduction must be character string\n')
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        stop(paste0('reduction: ', x, '  does not exist\n'))
        
      }
      
    }
    
  }
  
  if(is.null(n_pcs)) {
    
    n_pcs <- list()
    
    for(x in 1:length(reduction))  {
      
      n_pcs[[x]] <- 0
      
    }
    
  }
  
  if(!is.character(graph.name.suffix)) {
    
    stop('graph.name.suffix must be character string\n')
    
  }
  
  if(!is.character(batch)) {
    
    stop('batch must be character string\n')
    
  }
  
  for(r in batch) {
    
    if(!r %in% colnames(object@sample_metadata)) {
      
      stop(crayon::cyan(paste0(r, ' in batch does not exist\n')))
      
    }
    
  }
  
  if(!is.logical(approx)) {
    
    stop('approx must be logical: TRUE/FALSE \n')
    
  }
  
  if(!is.character(metric)) {
    
    stop('metric must be character \n')
    
  }
  
  if(!is.numeric(neighbors_within_batch)) {
    
    stop('neighbors_within_batch must be numerical \n')
    
  }
  
  if(!is.null(trim)) {
    
    if(!is.numeric(trim)) {
      
      stop('trim must be character string\n')
      
    }
    
  }
  
  if(!is.numeric(annoy_n_trees)) {
    
    stop('annoy_n_trees must be numerical\n')
    
  }
  
  if(!is.logical(use_faiss)) {
    
    stop('use_faiss must be logical: TRUE/FALSE\n')
    
  }
  
  if(!is.numeric(set_op_mix_ratio)) {
    
    stop('set_op_mix_ratio must be numerical\n')
    
  }
  
  if(!is.numeric(local_connectivity)) {
    
    stop('local_connectivity must be numerical\n')
    
  }
  
  if(!is.logical(generate.diffmap)) {
    
    stop('generate.diffmap should be boolean. either TRUE/FALSE\n')
    
  }
  
  if(isTRUE(generate.diffmap)) {
    
    if(!is.numeric(n_comps)) {
      
      stop('n_comps must be numerical \n')
      
    }
    
    if(!is.character(diffmap.name.suffix)) {
      
      stop('diffmap.name.suffix \n')
      
    }
    
  }
  
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(seed)) {
    
    stop('seed should be numerical\n')
    
  }
  
  sc <- reticulate::import('scanpy')
  pd <- reticulate::import('pandas')
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  reticulate::py_set_seed(seed, disable_hash_randomization = TRUE)
  
  if(!'integration_method' %in% colnames(object@pipelines)) {
    
    tmp <- tibble::add_column(.data = object@pipelines, integration_method=NA, integration_time=NA)
    
  } else {
    
    tmp <- object@pipelines
    
  }
  
  for(p in assay) {

    start_time <- Sys.time()
    
    count <- 1
    
    for(r in reduction) {
      
      dims <- n_pcs[[count]]
    
      scobj <- sc$AnnData(X = object@methods[[p]]@computational_reductions[[reduction[count]]])
      
      scobj$obs_names <- as.factor(colnames(object))
      
      scobj$var_names <- as.factor(colnames(object@methods[[p]]@computational_reductions[[reduction[count]]]))
      
      scobj$obsm$update(X_pca = object@methods[[p]]@computational_reductions[[reduction[count]]])
      
      if(length(colnames(as.data.frame(object@sample_metadata))) >= 1) {
        
        if(length(batch) > 1) {
          
          temp <- function(x) {
            
            return(paste(x, collapse = '_'))
            
          }
          
          df <- object@sample_metadata[,batch]
          df2 <- apply(X = df, MARGIN = 1, FUN = temp)
          df <- cbind(df, batch = df2)
          
          batch <- 'batch'
          
          scobj$obs <- pd$DataFrame(data = as.data.frame(df))  
          
        } else {
          
          scobj$obs <- pd$DataFrame(data = as.data.frame(object@sample_metadata))  
          
        }
        
      }
      
      if(dims == 0) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': npcs calculated\n')))
        
        dims <- ncol(object@methods[[p]]@computational_reductions[[reduction[count]]])
        
      }
      
      if(is.null(trim)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': initialising BBKNN for assay: ', p,  ', reduction: ', r, '\n')))
        
        sc$external$pp$bbknn(scobj,
                             batch_key = reticulate::r_to_py(as.character(batch)),
                             approx = as.logical(FALSE),
                             metric = as.character(metric),
                             neighbors_within_batch = as.integer(neighbors_within_batch),
                             n_pcs = dims,
                             annoy_n_trees = as.integer(annoy_n_trees),
                             use_faiss = as.logical(use_faiss),
                             set_op_mix_ratio = set_op_mix_ratio,
                             local_connectivity = local_connectivity)
        
        cat(crayon::cyan(paste0(Sys.time(), ': BBKNN complete \n')))
        
      } else if (!is.null(trim)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': initialising BBKNN for assay: ', p,  ', reduction: ', r, '\n')))
        
        sc$external$pp$bbknn(scobj,
                             batch_key= reticulate::r_to_py(as.character(batch)),
                             approx = as.logical(FALSE),
                             metric = as.character(metric),
                             neighbors_within_batch = as.integer(neighbors_within_batch),
                             n_pcs = dims,
                             trim = as.integer(trim),
                             annoy_n_trees = as.integer(annoy_n_trees),
                             use_faiss = as.logical(use_faiss),
                             set_op_mix_ratio = set_op_mix_ratio,
                             local_connectivity = local_connectivity)
        
        cat(crayon::cyan(paste0(Sys.time(), ': BBKNN complete \n')))
        
      }
      
      if(isTRUE(generate.diffmap)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': calcualting diffusion map\n')))
        
        sc$tl$diffmap(adata = scobj, n_comps = as.integer(n_comps))
        
        cat(crayon::cyan(paste0(Sys.time(), ': diffusion map calculated\n')))
        
        diffmap <- as.matrix(scobj$obsm[['X_diffmap']])
        
        DC_names <- list()
        
        counter <- 1
        
        for(x in 1:ncol(diffmap)) {
          
          DC_names[[counter]] <- paste0('DC_', counter)
          
          counter <- counter + 1
          
        }
        
        DC_names <- unlist(DC_names)
        
        colnames(diffmap) <- DC_names
        rownames(diffmap) <- colnames(object)
        
      }
      
      graph.list <- list()
      connectivities <- scobj$obsp[['connectivities']]
      colnames(connectivities) <- colnames(object)
      rownames(connectivities) <- colnames(object)
      distances <- scobj$obsp[['distances']]
      colnames(distances) <- colnames(object)
      rownames(distances) <- colnames(object)
      connectivities <- as(object = as.matrix(connectivities), Class = 'dgCMatrix')
      distances <- as(object = as.matrix(distances), Class = 'dgCMatrix')
      
      graph.list[['connectivities']] <- connectivities
      graph.list[['distances']] <- distances
      
      if('_' %in% unlist(x = strsplit(graph.name.suffix, split = ''))) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in graph.name.suffix, replacing with -\n')))
        
        graph.name.suffix <- sub(pattern = '_', replacement = '-', x = diffmap.name.suffix)
        
      }
      
      object@methods[[p]]@neighbours[[paste0(r, '_bbknn_bbknn', graph.name.suffix[count])]] <- graph.list
      
      cat(crayon::cyan(paste0(Sys.time(), ': bbknn results added to IBRAP object\n')))
      
      if(isTRUE(generate.diffmap)) {
        
        if('_' %in% (strsplit(x = diffmap.name.suffix, split = ''))) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in diffmap.name.suffix, replacing with -\n')))
          
          diffmap.name.suffix <- sub(pattern = '_', replacement = '-', x = diffmap.name.suffix)
          
        }
        
        object@methods[[p]]@computational_reductions[[paste0(r, '_bbknn_bbknn:diffmap', diffmap.name.suffix[count])]] <- diffmap
        
      }

      count <- count + 1
      
      end_time <- Sys.time()
      
      function_time <- end_time - start_time
      
      if(!'integration_method' %in% colnames(object@pipelines)) {
        
        tmp[which(x = tmp$normalisation_method==p),'integration_method'] <- paste0('BBKNN', graph.name.suffix)
        
        tmp[which(x = tmp$normalisation_method==p),'integration_time'] <- as.difftime(function_time, units = 'secs')
        
      }
      
      if('integration_method' %in% colnames(object@pipelines)) {
        
        if(paste0('BBKNN', graph.name.suffix) %in% tmp$integration_method) {
          
          tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('BBKNN', graph.name.suffix)),] <- c(tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('BBKNN', graph.name.suffix)),c('normalisation_method','normalisation_time')], paste0('BBKNN', graph.name.suffix), as.difftime(function_time, units = 'secs'))  
          
        }
        
        if(!paste0('BBKNN', graph.name.suffix) %in% object@pipelines$integration_method) {

          df <- tmp[which(tmp$normalisation_method==p),]

          df <- df[!duplicated(df$normalisation_method),]

          df[,'integration_method'] <- paste0('BBKNN', graph.name.suffix)

          df[,'integration_time'] <- function_time

          tmp <- rbind(tmp, df)
          
        }
        
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
