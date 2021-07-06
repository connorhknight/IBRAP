#' @name plot.dot.plot
#' @aliases plot.dot.plot
#' 
#' @title Plots a dot plot of gene expression
#'
#' @description Creates a dot plot of gene expression
#'
#' @param object An object of IBRAP class
#' @param assay Character. Which assay to access
#' @param slot Character. Which expression matrix slot should be used. Default = 'normalised'
#' @param clust.method Character. Which cluster method should be used, metadata accesses objects metadata 
#' @param column Character. Which column to access in the defined clust.method
#' @param features Character. Which features to plot
#'
#' @export plot.dot.plot

plot.dot.plot <- function(object, assay, slot='normalised', clust.method, column, features) {
  
  if(!is(object, 'IBRAP')) {
    
    stop('object must be of class IBRAP \n')
    
  }
  
  if(!is.character(assay)) {
    
    stop('assay must be character string \n')
    
  } else if (is.character(assay)) {
    
    if(!assay %in% names(object@methods)) {
      
      stop('assay is not contained within object@methods \n')
      
    }
    
  }
  
  if(!is.character(slot)) {
    
    stop('slot must be character string \n')
    
  } else if(is.character(slot)) {
    
    if(!slot %in% c('counts','normalised','norm.scaled')) {
      
      stop('slot must be either counts, normalised or norm.scaled \n')
      
    }
    
  }
  
  if(!is.character(clust.method)) {
    
    stop('clust.method must be character string \n')
    
  } 
  
  if(!is.character(column)) {
    
    stop('column must be character string \n')
    
  }
  
  if(clust.method != 'metadata') {
    
    if(!clust.method %in% names(object@methods[[assay]]@cluster_assignments)) {
      
      stop('clust.method not contained within cluster assignments \n')
      
    } else {
      
      if(!column %in% colnames(object@methods[[assay]]@cluster_assignments[[clust.method]][,column])) {
        
        stop('column not contained within clust.method \n')
        
      }
      
    }
    
  }
  
  for(x in features) {
    
    if(!x %in% rownames(object@methods[[assay]][[slot]])) {
      
      stop(paste0(x ' not contained within expression matrix \n'))
      
    }
    
  }
  
  expr <- object@methods[[assay]][[slot]]
  
  if(clust.method != 'metadata') {
    
    assignment <- object@methods[[assay]]@cluster_assignments[[clust.method]]
    
  } else if (clust.method == 'metadata') {
    
    assignment <- object@sample_metadata
    
  }
  
  expr <- as.matrix(expr[,match(x = colnames(expr), table = rownames(assignment))])
  features <- features[features %in% rownames(expr)]
  expr <- t(expr[features,])
  old.names <- colnames(expr)
  result <- cbind(expr, assignment[,column])
  colnames(result) <- c(old.names, 'variable')
  result <- as.data.frame(result)
  df <- result[1:length(unique(result$variable)),]
  df[,ncol(df)] <- unique(result$variable)
  df2 <- df
  
  count <- 1
  
  calc_nonzero <- function(x) {
    
    sum(sum(x != 0)/length(x))*100
    
  }
  
  for(d in df[,ncol(df)]) {
    
    sub_result <- result[result$variable == d,]
    
    apply(X = sub_result[,1:ncol(sub_result)-1], MARGIN = 2, FUN = mean)
    
    df[count,1:sum(ncol(df)-1)] <- apply(X = sub_result[,1:ncol(sub_result)-1], MARGIN = 2, FUN = mean)
    
    df2[count,1:sum(ncol(df)-1)] <- apply(X = sub_result[,1:ncol(sub_result)-1], MARGIN = 2, FUN = calc_nonzero)
    
    count <- count + 1
    
  }
  
  tmp <- data.frame(expression = numeric(), cluster = character(), percent_expressing = numeric(), gene = character())
  
  count <- 0
  count2 <- ncol(df)-1
  
  for(x in seq_along(1:nrow(df))) {
    
    tmp[sum(count + 1):count2,'gene expression'] <- as.numeric(t(df[x,1:sum(ncol(df)-1)])[,1])
    
    tmp[sum(count + 1):count2,'cluster assignments'] <- as.character(x)
    
    tmp[sum(count + 1):count2,'genes'] <- colnames(df)[1:sum(ncol(df)-1)]
    
    tmp[sum(count + 1):count2,'percentage of cells (%)'] <- as.numeric(t(df2[x,1:sum(ncol(df2)-1)])[,1])
    
    count <- count + ncol(df)-1
    count2 <- count2 + ncol(df)-1
    
  }
  
  ggplot2::ggplot(data = tmp, mapping = ggplot2::aes(x = `cluster assignments`, y = genes, col = `gene expression`, size = `percentage of cells (%)`)) + 
    ggplot2::geom_point() + 
    ggplot2::theme_classic() + 
    ggplot2::scale_color_gradient(low = 'gray93', high = 'mediumorchid4')
  
}
