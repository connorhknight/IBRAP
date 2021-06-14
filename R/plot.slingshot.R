#' @name plot.slingshot
#' @aliases plot.slingshot
#' 
#' @title Plots slingshot results
#'
#' @description Plots the results of the slingshot analysis using ggplots2
#' 
#' @param result A SlingshotDataSet class results object
#' @param clusters Vector. A vector of clusters to use in ggplot, if NULL the clusters used to generate results will be used. Default = NULL
#' @param lineages Boolean. Should the direct lineages be applied or curving, TRUE = curves, FALSE = lineages. Default = FALSE
#' @param which.lineage Numerical. If NULL, then all lineages are plotted, otherwise a numerical value is supplied indicating which lineage to use. Default = NULL
#' @param pt_size Numerical. What size should the cell points be. Default = 0.1
#' @param line_size Numerical. What size should the lineage lines be. Default = 0.1
#' @param title Character. Name the graph. Default = NULL
#' @param lab.clusters Boolean. Whether the clusters within a lineage should be labelled, note this only works for lineages and not curves. Default = TRUE
#' @param ... arguments to be passed to ggplot::geom_point
#' 
#' @return A ggplot of the reduced cellular embbedings and trajectories. 
#'
#' @export

plot.slingshot <- function(result, 
                           clusters = NULL, 
                           lineages = FALSE,
                           which.lineage = NULL,
                           pt_size = 0.1, 
                           line_size = 0.1,
                           title = NULL,
                           lab.clusters = TRUE,
                           ...) {
  
  if(isFALSE(is(object = result, 'SlingshotDataSet'))) {
    
    cat(crayon::cyan('result must be a SlingshotDataSet \n'))
    return(NULL)
    print('.')
  }
  
  if(!is.null(clusters)) {
    
    clusters <- as.character(clusters)
    
    if(length(clusters) != nrow(test@clusterLabels)) {
      
      cat(crayon::cyan('supplied clusters must be the same length as number of cells contained in results \n'))
      return(NULL)
      
    }
    
    if(!is.character(clusters)) {
      
      cat(crayon::cyan('unable to convert supplied clusters to character string \n'))
      return(NULL)
      
    }
    
  }
  
  if(!is.logical(lineages)) {
    
    cat(crayon::cyan('lineages should be boolean, TRUE/FALSE \n'))
    return(NULL)
    
  }
  
  if(!is.null(which.lineage)) {
    
    if(!is.numeric(which.lineage)) {
      
      cat(crayon::cyan('which.lineage should be numerical \n'))
      return(NULL)
      
    }
    
  }
  
  if(!is.numeric(pt_size)) {
    
    cat(crayon::cyan('pt_size msut be numerical \n'))
    return(NULL)
    
  }
  
  if(!is.numeric(line_size)) {
    
    cat(crayon::cyan('line_size msut be numerical \n'))
    return(NULL)
    
  }
  
  if(!is.null(title)) {
    
    if(!is.character(title)) {
      
      cat(crayon::cyan('title must be character string \n'))
      return(NULL)
      
    }
    
  }
  
  if(!is.logical(lab.clusters)) {
    
    
    
  }
  
  red <- slingshot::reducedDim(result)[,1:2]
  
  if(is.null(clusters)) {
    
    clusters <- slingshot::slingClusterLabels(result)
    clusters <- apply(clusters, 1, which.max)
    
  }
  
  df <- data.frame(dim1 = red[,1], dim2 = red[,2], cluster = as.character(clusters))
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = dim1, y = dim2, col = cluster)) + 
    ggplot2::geom_point(size = pt_size, ...) + 
    ggplot2::theme(legend.title.align=0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = colorspace::qualitative_hcl(n = length(unique(df[,3])), palette = 'Dark 3'))
  
  clust_centres <- data.frame(clusters = unique(df$cluster))
  
  centre_1 <- list()
  centre_2 <- list()
  
  count <- 1
  
  for(x in clust_centres$clusters) {
    
    centre_1[[count]] <- mean(df[df[,'cluster'] == x,][,'dim1'])
    centre_2[[count]] <- mean(df[df[,'cluster'] == x,][,'dim2'])
    
    count <- count + 1
    
  }
  
  clust_centres[,'dim1'] <- unlist(centre_1)
  clust_centres[,'dim2'] <- unlist(centre_2)
  
  if(is.null(which.lineage)) {
    
    if (isFALSE(lineages)) {
      
      for(i in seq_along(slingshot::slingCurves(result))) {
        
        curve_i <- slingshot::slingCurves(result)[[i]]
        curve_i <- curve_i$s[curve_i$ord, ][,1:2]
        colnames(curve_i) <- c("dim1", "dim2")
        p <- p + ggplot2::geom_path(data = as.data.frame(curve_i), col = "black", size = line_size)
        
      }
      
      if(!is.null(title)) {
        
        p <- p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
      }
      
    } else if (isTRUE(lineages)) {
      
      for(i in seq_along(slingshot::slingLineages(result))) {
        
        lineage_i <- slingshot::slingLineages(result)[[i]]
        coord <- clust_centres[clust_centres[lineage_i,'clusters'],]
        p <- p + ggplot2::geom_line(data = coord, mapping = ggplot2::aes(group = 1), color = 'black')
        
      }
      
      if(!is.null(title)) {
        
        p <- p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
      }
      
      if(isTRUE(lab.clusters)) {
        
        p <- p + ggrepel::geom_label_repel(mapping = 
                                             ggplot2::aes(x = dim1, y = dim2, label = clusters), 
                                           data = clust_centres, color = 'black', force = 100)
        
      }
      
    }
    
  } else if (!is.null(which.lineage)) {
    
    if (isFALSE(lineages)) {
      
      curve_i <- slingshot::slingCurves(result)[[which.lineage]]
      
      curve_i <- curve_i$s[curve_i$ord, ][,1:2]
      
      colnames(curve_i) <- c("dim1", "dim2")
      
      p <- p + ggplot2::geom_path(data = as.data.frame(curve_i), col = "black", size = line_size)
      
      if(!is.null(title)) {
        
        p <- p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
      }
      
    } else if (isTRUE(lineages)) {
      
      lineage_i <- slingshot::slingLineages(result)[[which.lineage]]
      coord <- clust_centres[clust_centres[lineage_i,'clusters'],]
      p <- p + ggplot2::geom_line(data = coord, mapping = ggplot2::aes(group = 1), color = 'black')
      
      if(isTRUE(lab.clusters)) {
        
        p <- p + ggrepel::geom_label_repel(mapping = 
                                                   ggplot2::aes(x = dim1, y = dim2, label = clusters), 
                                                 data = clust_centres[clust_centres[lineage_i,'clusters'],], color = 'black', force = 100)
        
      }
      
      if(!is.null(title)) {
        
        p <- p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        
      }
      
    }
    
  }
  
  print(p)
  
}
