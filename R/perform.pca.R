#' @name perform.pca
#' @aliases perform.pca
#' 
#' @title Performs PCA reduction
#'
#' @description Performs PCA reduction on defined method-assays. Data should be HVG subset, normalised and scaled (in the norm.scaled assay)
#' 
#' @param object IBRAP S4 class object
#' @param assay Character. String containing indicating which assay to use
#' @param slot Character. String indicating which slot within the assay should be sourced
#' @param n.pcs Numerical. How many principal components should be produced. Default = 50
#' @param reduction.save Character. What should this reduction be saved as in computation_reduction. Default = 'pca'
#' @param ... Arguments to be passed to PCAtools::pca
#' 
#' @return PCA reductions contained within the computational_reduction list in the defined assays
#' 
#' @examples 
#' 
#' object <- perform.pca(object = object, 
#'                       assay = c('SCT', 'SCRAN', 'SCANPY'), 
#'                       n.pcs = 50, 
#'                       reduction.save = 'pca')
#'
#' @export

perform.pca <- function(object, 
                        assay,
                        slot='norm.scaled',
                        n.pcs=50,
                        reduction.save='pca', 
                        save.plot = TRUE, 
                        ...) {
  
  if(!is(object = object, class2 = 'IBRAP')) {
    
    stop('object must be of class IBRAP')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string')
    
  }
  
  for(x in assay) {
    
    if(!x %in% names(object@methods)) {
      
      stop(paste0('assay: ', x, ' does not exist\n'))
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be character string')
    
  }
  
  if(!is.numeric(n.pcs)) {
    
    stop('n.pcs must be numerical')
    
  }
  
  if(!is.character(reduction.save)) {
    
    stop('reduction.save must be numerical')
    
  }
  
  if(!is.logical(save.plot)) {
    
    stop('save.plot must be boolean. TRUE/FALSE \n')
    
  }
  
  ggarrange.tmp <- function(...) {
    
    egg::ggarrange(...)
    
  }
  
  list.of.figs <- list()
  
  for(t in assay) {
    
    if(is.null(object@methods[[t]]@highly.variable.genes)) {
      
      stop(paste0('no variable features have been identified for assay: ', t))
      
    }
    
    mat <- as.matrix(object@methods[[t]][[slot]][rownames(object@methods[[t]][[slot]]) %in% object@methods[[t]]@highly.variable.genes,])
    
    cat(crayon::cyan(paste0(Sys.time(), ': initialising PCA for assay:', t, '\n')))
    
    ass <- strsplit(x = names(object@methods)[which(names(object@methods)==t)], split = '_')[[1]][1]
    
    if(ass %in% c('SCT','SCRAN','TPM')) {
      
      a <- PCAtools::pca(mat = mat, center = F, scale = F, ...)
      eig <- a$sdev^2/sum(a$sdev^2)
      eig <- as.data.frame(eig*100)
      temp <- as.character(colnames(a$rotated))
      eig[,2] <- factor(x = temp, levels = unique(temp))
      colnames(eig) <- c ('Variance', 'PCs')
      
      p <- ggplot2::ggplot(data = eig[1:n.pcs,], mapping = ggplot2::aes(x = PCs, y = Variance)) + 
        ggplot2::geom_point() + 
        egg::theme_article() + 
        ggplot2::ylab('Explained Variance (%)') +
        ggplot2::ggtitle(t) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = ggplot2::element_text(hjust = 0.5))
      
      cat(crayon::cyan(paste0(Sys.time(), ': PCA completed\n')))
      
      object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,1:n.pcs])
      
    } else if (ass == 'SCANPY') {
      
      sc <- reticulate::import('scanpy')
      
      scobj <- sc$AnnData(X = t(as.matrix(object@methods[[t]][['norm.scaled']])))
      scobj$obs_names <- as.factor(colnames(object@methods[[t]][['norm.scaled']]))
      scobj$var_names <- as.factor(rownames(object@methods[[t]][['norm.scaled']]))
      
      sc$tl$pca(data = scobj, n_comps = as.integer(n.pcs), use_highly_variable = as.logical(F))
      
      tmp <- scobj$obsm[['X_pca']]
      rownames(tmp) <- colnames(object)
      
      pc.names <- list()
      count <- 1
      
      for(x in 1:n.pcs) {
        
        pc.names[[count]] <- paste0('PC',count)
        
        count <- count + 1
        
      }
      
    colnames(tmp) <- unlist(pc.names)
    rownames(tmp) <- colnames(object@methods[[t]][['norm.scaled']])
    
    eig <- as.data.frame(apply(X = tmp, MARGIN = 2, FUN = sd)^2/sum(apply(X = tmp, MARGIN = 2, FUN = sd)^2)*100)
    eig[,2] <- factor(x = colnames(tmp), levels = unique(colnames(tmp)))
    colnames(eig) <- c('Variance', 'PCs')
    
    p <- ggplot2::ggplot(data = eig[1:n.pcs,], mapping = ggplot2::aes(x = PCs, y = Variance)) + 
      ggplot2::geom_point() + 
      egg::theme_article() + 
      ggplot2::ylab('Explained Variance (%)') +
      ggplot2::ggtitle(t) + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = ggplot2::element_text(hjust = 0.5))
    
    cat(crayon::cyan(paste0(Sys.time(), ': PCA completed\n')))
    
    object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(tmp)
      
    } else {
      
      a <- PCAtools::pca(mat = mat, center = F, scale = F, ...)
      eig <- a$sdev^2/sum(a$sdev^2)
      eig <- as.data.frame(eig*100)
      temp <- as.character(colnames(a$rotated))
      eig[,2] <- factor(x = temp, levels = unique(temp))
      colnames(eig) <- c('Variance', 'PCs')
      
      p <- ggplot2::ggplot(data = eig[1:n.pcs,], mapping = ggplot2::aes(x = PCs, y = Variance)) + 
        ggplot2::geom_point() + 
        egg::theme_article() + 
        ggplot2::ylab('Explained Variance (%)') +
        ggplot2::ggtitle(t) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = ggplot2::element_text(hjust = 0.5))
      
      cat(crayon::cyan(paste0(Sys.time(), ': PCA completed\n')))
      
      object@methods[[t]]@computational_reductions[[reduction.save]] <- as.matrix(a$rotated[,1:n.pcs])
      
    }
    
    if(isTRUE(save.plot)) {
      
      pdf(file = paste0('PCA_', t, '.pdf'), onefile = TRUE)
      
    }
    
    print(p)
    
    if(isTRUE(save.plot)) {
      
      dev.off()
      
    }

  }
  
  return(object)
  
}