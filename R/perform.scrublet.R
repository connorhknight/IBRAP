#' @name perform.scrublet
#' @aliases perform.scrublet
#' 
#' @title Python module: scrublet
#'
#' @description Removes doublets from dataset.
#' 
#' @param counts Counts matrix
#' @param total_counts Total number of cells. NULL = automatically counts.
#' @param sim_doublet_ratio Number of doublets to simulate relative to observed
#' @param n_neighbors Expected number of neighbours per cell
#' @param expected_doublet_rate Expected percentage of doublets to be present in the dataset
#' @param stdev_doublet_rate Uncertainty in expected doublet rate
#' @param random_state Random state for doublet simulation, approximate nearest neighbour search, nd PCA/Truncated PCA
#' @param synthetic_doublet_umi_subsampling Sampling rate for UMIs in a cell when synthesising doublets
#' @param use_approx_neighbors Use approximate nearest neighbor method `(annoy)` for the KNN classifier
#' @param distance_metric Define distance metric for nearest neighbour calculation: 'angular', 'euclidean', 'manhattan', 'hamming', 'dot'.
#' @param get_doublet_neighbor_parents return the transcriptomes of the parent cells for simulated doublets
#' @param min_counts Minimum counts per cell
#' @param min_cells Minimum number of cells per gene
#' @param min_gene_variability_pctl Variability cutoff when deducing highly variable genes prior to PCA reduction
#' @param log_transform Log transforms the data 
#' @param mean_center Should the dataset be centred around the mean
#' @param normalize_variance Should the genes have a total variance of 1
#' @param n_prin_comps Number of principal components to retain
#' @param svd_solver Which SVD solver to use: 'auto', 'full', 'arpack', 'randomized'.
#' 
#' @usage perform.scrublet(counts = counts, expected_doublet_rate = 0.025)
#' 
#' @return Doublet-omitted sparse matrix
#'
#' @export

perform.scrublet <- function(counts,
                             total_counts = NULL, 
                             sim_doublet_ratio = 2.0, 
                             n_neighbors = NULL, 
                             expected_doublet_rate = 0.075, 
                             stdev_doublet_rate = 0.02, 
                             random_state = 0L,
                             synthetic_doublet_umi_subsampling = 1.0,
                             use_approx_neighbors = TRUE,
                             distance_metric = 'euclidean',
                             get_doublet_neighbor_parents = FALSE,
                             min_counts = 3L,
                             min_cells = 3L,
                             min_gene_variability_pctl = 85L,
                             log_transform = FALSE,
                             mean_center = TRUE, 
                             normalize_variance = TRUE,
                             n_prin_comps = 30L,
                             svd_solver = 'arpack') {
  
  if(!is(object = counts, class2 = 'matrix')) {
    
    if (!is(object = counts, class2 = 'dgCMatrix')) {
      
      cat(crayon::cyan('counts must be in matrix or dgCMatrix format\n'))
      return(counts)
      
    }
    
  } 
  
  if(!is.null(total_counts)) {
    
    if(!is.numeric(total_counts)) {
      
      cat(crayon::cyan('total_counts must be numerical\n'))
      return(counts)
      
    }
    
  }
  
  if(!is.numeric(sim_doublet_ratio)) {
    
    cat(crayon::cyan('sim_doublet_ratio must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.null(n_neighbors)) {
    
    if(!is.numeric(n_neighbors)) {
      
      cat(crayon::cyan('n_neighbors must be numerical\n'))
      return(counts)
      
    }
    
  }
  
  if(!is.numeric(expected_doublet_rate)) {
    
    cat(crayon::cyan('expected_doublet_rate must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(stdev_doublet_rate)) {
    
    cat(crayon::cyan('stdev_doublet_rate must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(random_state)) {
    
    cat(crayon::cyan('random_state must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(synthetic_doublet_umi_subsampling)) {
    
    cat(crayon::cyan('synthetic_doublet_umi_subsampling must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.logical(use_approx_neighbors)) {
    
    cat(crayon::cyan('use_approx_neighbors must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.character(distance_metric)) {
    
    cat(crayon::cyan('distance_metric must be character string\n'))
    return(counts)
    
  }
  
  if(!is.logical(get_doublet_neighbor_parents)) {
    
    cat(crayon::cyan('get_doublet_neighbor_parents must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.numeric(min_counts)) {
    
    cat(crayon::cyan('min_counts must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(min_cells)) {
    
    cat(crayon::cyan('min_cells must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.numeric(min_gene_variability_pctl)) {
    
    cat(crayon::cyan('min_gene_variability_pctl must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.logical(log_transform)) {
    
    cat(crayon::cyan('log_transform must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.logical(mean_center)) {
    
    cat(crayon::cyan('mean_center must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.logical(normalize_variance)) {
    
    cat(crayon::cyan('normalize_variance must be logical: TRUE/FALSE\n'))
    return(counts)
    
  }
  
  if(!is.numeric(n_prin_comps)) {
    
    cat(crayon::cyan('n_prin_comps must be numerical\n'))
    return(counts)
    
  }
  
  if(!is.character(svd_solver)) {
    
    cat(crayon::cyan('svd_solver must be numerical\n'))
    return(counts)
    
  }
  
  cat(crayon::cyan('Initialising scrublet\n'))
  scrublet <- reticulate::import('scrublet', convert = FALSE)
  cat(crayon::cyan('Python modules loaded\n'))
  
  scrub1 <- scrublet$Scrublet(counts_matrix = as.data.frame(as.matrix(t(counts))), 
                              total_counts = total_counts, 
                              sim_doublet_ratio = sim_doublet_ratio, 
                              n_neighbors = n_neighbors, 
                              expected_doublet_rate = expected_doublet_rate, 
                              stdev_doublet_rate = stdev_doublet_rate, 
                              random_state = random_state)
  
  cat(crayon::cyan('scrublet object created\n'))
  
  res1 <- reticulate::py_to_r(scrub1$scrub_doublets(synthetic_doublet_umi_subsampling = synthetic_doublet_umi_subsampling,
                                                    use_approx_neighbors = use_approx_neighbors, 
                                                    distance_metric = distance_metric, 
                                                    get_doublet_neighbor_parents = get_doublet_neighbor_parents, 
                                                    min_counts = min_counts,
                                                    min_cells = min_cells, 
                                                    min_gene_variability_pctl = min_gene_variability_pctl,
                                                    log_transform = log_transform,
                                                    mean_center = mean_center,
                                                    normalize_variance = normalize_variance,
                                                    n_prin_comps = n_prin_comps,
                                                    svd_solver = svd_solver,
                                                    verbose = TRUE))
  
  sim.plot <- ggplot2::qplot(as.vector(reticulate::py_to_r(scrub1$doublet_scores_sim_)), 
                             geom = 'histogram') + 
    ggplot2::stat_bin(bins = 100) + 
    ggplot2::xlab('doublet scores') + 
    ggplot2::ylab('frequency') + 
    ggplot2::ggtitle(paste0('simulated_doublets')) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
  
  obs.plot <- ggplot2::qplot(as.vector(res1)[[1]], 
                             geom = 'histogram') + 
    ggplot2::stat_bin(bins = 80) + 
    ggplot2::xlab('doublet scores') + 
    ggplot2::ylab('frequency') + 
    ggplot2::ggtitle(paste0('observed doublets')) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
  
  comb.plot <- cowplot::plot_grid(sim.plot, obs.plot, ncol = 2, nrow = 1)
  print(comb.plot)
  
  cat(crayon::cyan('doublets detected\n'))
  counts <- as.matrix(counts)
  counts <- counts[,!res1[[2]]]
  counts <- Matrix::Matrix(data = counts, sparse = T)
  cat(crayon::cyan('matrix scrubbed\n'))
  
  return(counts)
}