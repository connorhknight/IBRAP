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
#' @param pt_size Numerical. What size should the cell points be. Default = 0.1
#' @param line_size Numerical. What size should the lineage lines be. Default = 0.1
#' @param relevant Should only pseudotime relevant cells be plotted. Default = TRUE
#' @param title Character. Name the graph. Default = NULL
#' @param lab.clusters Boolean. Whether the clusters within a lineage should be labelled, note this only works for lineages and not curves. Default = TRUE
#' @param ... arguments to be passed to ggplot::geom_point
#' 
#' @return A ggplot of the reduced cellular embbedings and trajectories. 
#'
#' @export plot.slingshot

plot.slingshot <- function(result, 
                           clusters = NULL, 
                           lineages = FALSE,
                           pt_size = 0.1, 
                           line_size = 0.1,
                           relevant = TRUE,
                           Pseudotime = FALSE,
                           Expression = FALSE,
                           object = NULL,
                           assay = NULL,
                           slot = 'norm.scaled',
                           feature = NULL,
                           ...) {
  
  if(isFALSE(is(object = result$assignments, 'SlingshotDataSet'))) {
    
    stop('result must be a SlingshotDataSet \n')

  }
  
  if(!is.null(clusters)) {
    
    clusters <- as.character(clusters)
    
    if(length(clusters) != nrow(test@clusterLabels)) {
      
      stop('supplied clusters must be the same length as number of cells contained in results \n')
      
    }
    
    if(!is.character(clusters)) {
      
      stop('unable to convert supplied clusters to character string \n')
      
    }
    
  }
  
  if(!is.logical(lineages)) {
    
    stop('lineages should be boolean, TRUE/FALSE \n')
    
  }
  
  if(!is.numeric(pt_size)) {
    
    stop('pt_size msut be numerical \n')
    
  }
  
  if(!is.numeric(line_size)) {
    
    stop('line_size msut be numerical \n')
    
  }
  
  if(is.null(clusters)) {
    
    clusters <- slingshot::slingClusterLabels(result$assignments)
    clusters <- apply(clusters, 1, which.max)
    
  }
  
  if(!is.null(Pseudotime)) {
    
    pseudotime <- result$pseudotimes
    
  }
  
  if(!feature %in% rownames(object@methods[[assay]][[slot]])) {
    
    cat(crayon::cyan(paste0(Sys.time(), ': ',feature, ' not present in ', slot, ', switching to count matrix \n')))
    slot <- 'counts'
    
  }
  
  plot.list <- list()
  
  red_whole <- slingshot::reducedDim(result$assignments)[,1:2]

  red_whole <- cbind(red_whole,clusters)

  for(x in 1:length(result$assignments@curves)) {
    
    red <- slingshot::reducedDim(result$assignments)[,1:2]
    
    if(isTRUE(relevant)) {
      
      rel.cells <- rownames(result$pseudotimes[!is.na(result$pseudotimes[,x]),])
      
      red <- cbind(red,clusters)
      
      red <- red[rel.cells,]
      
      pseudotime_x <- pseudotime[rel.cells,x]
      
      expr <- object@methods[[assay]][[slot]][feature,rel.cells]
      
    } else {
      
      pseudotime_x <- pseudotime[,x]
      pseudotime_x[is.na(pseudotime_x)] <- 0
      expr <- object@methods[[assay]][[slot]][feature,]
      red <- cbind(red,clusters)
      
    }
    
    clust_centres <- data.frame(clusters = unique(red[,3]))
    
    centre_1 <- list()
    centre_2 <- list()
    
    count <- 1
    
    for(d in clust_centres$clusters) {
      
      centre_1[[count]] <- mean(red_whole[red_whole[,'clusters'] == d,][,1])
      centre_2[[count]] <- mean(red_whole[red_whole[,'clusters'] == d,][,2])
      
      count <- count + 1

    }
    
    clust_centres[,'dim1'] <- unlist(centre_1)
    clust_centres[,'dim2'] <- unlist(centre_2)
    
    if(isFALSE(Pseudotime) && isFALSE(Expression)) {

      df <- data.frame(dim1 = red[,1], dim2 = red[,2], cluster = as.character(red[,3]))
      
      p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = dim1, y = dim2, col = cluster)) + 
        ggplot2::geom_point(size = pt_size, ...) + 
        ggplot2::theme_classic() +
        ggplot2::theme(legend.position = 'none') +
        ggplot2::scale_color_manual(values = colorspace::qualitative_hcl(n = length(unique(df[,3])), palette = 'Dark 3'))
      
      if (isFALSE(lineages)) {
        
        curve_i <- slingshot::slingCurves(result$assignments)[[x]]
        curve_i <- curve_i$s[curve_i$ord, ][,1:2]
        colnames(curve_i) <- c("dim1", "dim2")
        p <- p + ggplot2::geom_path(data = as.data.frame(curve_i), col = "black", size = line_size)
        
      } else if (isTRUE(lineages)) {
        
        print(x)
        
        lineage_i <- slingshot::slingLineages(result$assignments)[[x]]
        
        print(lineage_i)
        
        print(clust_centres)
        
        coord <- clust_centres[lineage_i %in% clust_centres[,'clusters'],]

        print(coord)
        
        p <- p + ggplot2::geom_line(data = coord, mapping = ggplot2::aes(group = 1), color = 'black')
        
      }
      
      plot.list[[x]] <- p
      
    } else if(isTRUE(Pseudotime) && isFALSE(Expression)) {
      
      df <- data.frame(dim1 = red[,1], dim2 = red[,2], pseudotime = pseudotime_x, cluster = as.character(red[,3]))
      p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = dim1, y = dim2, col = pseudotime)) + 
        ggplot2::geom_point(size = pt_size, ...) + 
        ggplot2::theme(legend.title.align=0.5) +
        ggplot2::theme_classic() +
        ggplot2::scale_color_continuous()
      
    if (isFALSE(lineages)) {
      
      curve_i <- slingshot::slingCurves(result$assignments)[[x]]
      curve_i <- curve_i$s[curve_i$ord, ][,1:2]
      colnames(curve_i) <- c("dim1", "dim2")
      p <- p + ggplot2::geom_path(data = as.data.frame(curve_i), col = "black", size = line_size)
      
    } else if (isTRUE(lineages)) {
      
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
      
      lineage_i <- slingshot::slingLineages(result$assignments)[[x]]
      coord <- clust_centres[clust_centres[lineage_i,'clusters'],]
      p <- p + ggplot2::geom_line(data = clust_centres, mapping = ggplot2::aes(group = 1), color = 'black')
      
    }

    plot.list[[x]] <- p
    
  } else if (isFALSE(Pseudotime) && isTRUE(Expression)) {
    
    df <- data.frame(dim1 = red[,1], dim2 = red[,2], expression = expr, cluster = as.character(red[,3]))
    p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = dim1, y = dim2, col = expression)) + 
      ggplot2::geom_point(size = pt_size, ...) + 
      ggplot2::theme(legend.title.align=0.5) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_continuous()
    
    if (isFALSE(lineages)) {
      
      curve_i <- slingshot::slingCurves(result$assignments)[[x]]
      curve_i <- curve_i$s[curve_i$ord, ][,1:2]
      colnames(curve_i) <- c("dim1", "dim2")
      p <- p + ggplot2::geom_path(data = as.data.frame(curve_i), col = "black", size = line_size)
      
    } else if (isTRUE(lineages)) {
      
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
      
      lineage_i <- slingshot::slingLineages(result$assignments)[[x]]
      coord <- clust_centres[clust_centres[lineage_i,'clusters'],]
      p <- p + ggplot2::geom_line(data = clust_centres, mapping = ggplot2::aes(group = 1), color = 'black')
      
      }
    
    plot.list[[x]] <- p
    
    }
    
  }

  return(plot.list)
  
}
