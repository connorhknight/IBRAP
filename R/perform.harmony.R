#' @name perform.harmony
#' @aliases perform.harmony
#' 
#' @title Performs Harmony integration
#'
#' @description Performs Harmony integration on defined method-assays and reductions contained within. This is performed on reductions. 
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param batch Character. A string of the column nmae that contains variables to regress. 
#' @param reduction Character. String defining the name of the reduction to provide for HARMONY. Default = 'pca'
#' @param reduction.save.suffix Character. What should be appended to the end of harmony as the reduction name
#' @param dims.use Numerical. Number of dimensions of the provided reduction to input into harmony, NULL equates to all dimensions. Default = NULL
#' @param print.harmony.plot Boolean. Should the automatically generated plot be printed? Default = FALSE
#' @param verbose Logical Should function messages be printed?
#' @param seed Numeric. What should the seed be set as. Default = 1234
#' @param ... Arguments to be passed to harmony::HarmonyMatrix
#' 
#' @return PCA reductions contained within the computational_reduction list in the defined assays
#' 
#' @examples
#' 
#' object <- perform.harmony(object = object, 
#'                           assay = c('SCRAN', 'SCT', 'SCANPY'), 
#'                           batch = 'original.project', 
#'                           reduction = c('pca'),  
#'                           max.iter.harmony = 100,
#'                           dims.use = list(NULL))
#'
#' @export

perform.harmony <- function(object, 
                            assay, 
                            batch, 
                            reduction = 'pca', 
                            reduction.save.suffix = '',
                            dims.use = NULL, 
                            print.harmony.plot = FALSE,
                            seed=1234,
                            verbose = FALSE,
                            ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string\n')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('reduction: ', x, 'does not exist\n'))
      
    }
    
  }
  
  if(!is.character(batch)) {
    
    stop('batch must be character string(s)\n')
    
  }
  
  if(!batch %in% names(object@sample_metadata)) {
    
    stop('batch does not exist\n')
    
  }
  
  if(!is.character(reduction)) {
    
    stop('reduction must be character string(s)\n')
    
  }
  
  for(x in reduction) {
    
    for(i in assay) {
      
      if(!x %in% names(c(object@methods[[i]]@computational_reductions, 
                         object@methods[[i]]@visualisation_reductions, 
                         object@methods[[i]]@integration_reductions))) {
        
        stop(paste0('reduction: ', x, ' does not exist\n'))
        
      }
      
    }
    
  }
  
  if(!is.character(reduction.save.suffix)) {
    
    stop('reduction.save.suffix must be character string(s)\n')
    
  }
  
  if(is.list(dims.use)) {
    
    for(x in dims.use) {
      
      if(!is.null(x)) {
        
        if(!is.numeric(x)) {
          
          stop('dims must either be numeric or NULL\n')
          
        } 
        
      } else if (!is.numeric(x)) {
        
        if(!is.null(x)) {
          
          stop('dims must either be numeric or NULL\n')
          
        }
        
      }
      
    }
    
  }
  
  if(is.null(dims.use)) {
    
    dims.use <- list()
    
    for(x in 1:length(reduction))  {
      
      dims.use[[x]] <- 0
      
    }
    
  }
  
  if(!is.logical(verbose)) {
    
    stop('verbose should be logical, TRUE/FALSE \n')
    
  } 
  
  if(!is.numeric(seed)) {
    
    stop('seed must be a numerical value \n')
    
  }
  
  set.seed(seed = seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  
  # if(!'integration_method' %in% colnames(object@pipelines)) {
  #   
  #   tmp <- tibble::add_column(.data = object@pipelines, integration_method=NA, integration_time=NA)
  #   
  # } else {
  #   
  #   tmp <- object@pipelines
  #   
  # }
  
  for(p in assay) {
    
    start_time <- Sys.time()
    
    reduction.list <- list()
    red.names <- c(names(object@methods[[p]]@computational_reductions), 
                   names(object@methods[[p]]@integration_reductions),
                   names(object@methods[[p]]@visualisation_reductions))
    
    for(i in red.names) {
      
      if(i %in% names(object@methods[[p]]@computational_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@computational_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@integration_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@integration_reductions[[i]]
        
      }
      
      if(i %in% names(object@methods[[p]]@visualisation_reductions)) {
        
        reduction.list[[i]] <- object@methods[[p]]@visualisation_reductions[[i]]
        
      }
      
    }
    
    for(r in reduction) {
      
      if(!r %in% names(reduction.list)) {
        
        stop('reductions could not be found\n')
        
      }
      
    }
    
    count <- 1
    
    for(g in reduction) {
      
      red <- reduction.list[[g]]
      
      dims <- dims.use[[count]]
      
      if(dims == 0) {
        
        dims <- ncol(red)
        
      }
      
      if('_' %in% unlist(strsplit(x = reduction.save.suffix, split = ''))) {
        
        if(isTRUE(verbose)) {
          
          cat(crayon::cyan(paste0(Sys.time(), ': _ cannot be used in reduction.save.suffix, replacing with - \n')))
          
        }
        
        reduction.save.suffix <- sub(pattern = '_', replacement = '-', x = reduction.save.suffix)
        
      }
      
      if(isTRUE(verbose)) {
        
        cat(crayon::cyan(paste0(Sys.time(), ': initialising harmony for assay: ', p, ', reduction: ', g, '\n')))
        
      }
      
      if(isTRUE(print.harmony.plot)) {
        
        if(isTRUE(verbose)) {
          
          harm <- harmony::HarmonyMatrix(data_mat = red[,1:dims], meta_data = object@sample_metadata, vars_use = batch, do_pca = FALSE, verbose = TRUE, plot_convergence = TRUE, ...)
          
        }
        
        if(isFALSE(verbose)) {
          
          harm <- harmony::HarmonyMatrix(data_mat = red[,1:dims], meta_data = object@sample_metadata, vars_use = batch, do_pca = FALSE, verbose = FALSE, plot_convergence = TRUE, ...)
          
        }
        
      }
      
      if(isFALSE(print.harmony.plot)) {
        
        if(isTRUE(verbose)) {
          
          harm <- harmony::HarmonyMatrix(data_mat = red[,1:dims], meta_data = object@sample_metadata, vars_use = batch, do_pca = FALSE, verbose = TRUE, plot_convergence = FALSE, ...)
          
        }
        
        if(isFALSE(verbose)) {
          
          harm <- harmony::HarmonyMatrix(data_mat = red[,1:dims], meta_data = object@sample_metadata, vars_use = batch, do_pca = FALSE, verbose = FALSE, plot_convergence = FALSE, ...)
          
        }
        
      }
      
      object@methods[[p]]@integration_reductions[[paste0(r, '_HARMONY', reduction.save.suffix)]] <- harm
      
      count <- count + 1
      
      end_time <- Sys.time()
      
      function_time <- end_time - start_time
      
      # if(!'integration_method' %in% colnames(object@pipelines)) {
      #   
      #   tmp[which(x = tmp$normalisation_method==p),'integration_method'] <- paste0('HARMONY', reduction.save.suffix)
      #   
      #   tmp[which(x = tmp$normalisation_method==p),'integration_time'] <- as.difftime(function_time, units = 'secs')
      #   
      # }
      # 
      # if('integration_method' %in% colnames(object@pipelines)) {
      #   
      #   if(paste0('HARMONY', reduction.save.suffix) %in% tmp$integration_method) {
      #     
      #     tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('HARMONY', reduction.save.suffix)),] <- c(tmp[which(tmp$normalisation_method==p & tmp$integration_method==paste0('HARMONY', reduction.save.suffix)),c('normalisation_method','normalisation_time')], paste0('HARMONY', reduction.save.suffix), as.difftime(function_time, units = 'secs'))  
      #     
      #   }
      #   
      #   if(!paste0('HARMONY', reduction.save.suffix) %in% object@pipelines$integration_method) {
      #     
      #     df <- tmp[which(tmp$normalisation_method==p),]
      #     
      #     df <- df[!duplicated(df$normalisation_method),]
      #     
      #     df[,'integration_method'] <- paste0('HARMONY', reduction.save.suffix)
      #     
      #     df[,'integration_time'] <- function_time
      #     
      #     tmp <- rbind(tmp, df)
      #     
      #   }
      #   
      # }
      
    }
    
  }
  
  if(isTRUE(verbose)) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': harmony completed\n')))
    
  }
  
  # tmp$integration_time <- as.difftime(tim = tmp$integration_time, units = 'secs')
  # 
  # rownames(tmp) <- 1:nrow(tmp)
  # 
  # object@pipelines <- tmp
  
  
  return(object)
}