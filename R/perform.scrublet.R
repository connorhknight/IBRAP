#' Perform Scrublet - a python-based package that omits doublets
#'
#' @param object A SingleCellExperiment
#' @param split.by Which metadata should be used to split the data, i.e. individual samples
#' @param min_counts Minimum number of counts to be present in a cell
#' @param min_cells Minimum number of cells to be present in a feature
#' @param ... 
#' @keywords
#' @export
#' @examples Read10X_output(directory = 'Users/user/10x_output')
#'
#'

perform.scrublet <- function(object, 
                             split.by, 
                             min_counts = 1,
                             min_cells = 1,
                             min_gene_variability_pctl = 85,
                             ...) {
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
      res1 <- scrub1$scrub_doublets(min_counts = min_counts, 
                                    min_cells = min_cells, 
                                    min_gene_variability_pctl = min_gene_variability_pctl, 
                                    verbose = TRUE,
                                    ...)
      cat(crayon::cyan('doublets detected\n'))
      raw_counts <- t(as.data.frame(raw_counts))
      raw_counts <- raw_counts[,!reticulate::py_to_r(res1)[[2]] == TRUE]
      cat(crayon::cyan('matrix scrubbed\n'))
      raw_counts_list[[l]] <- raw_counts
    }
    raw_counts <- do.call('cbind', raw_counts_list)
    object <- object[,colnames(raw_counts)]
    assay(object, 'counts') <- raw_counts
    return(object)
    rm(obj, raw_counts, scrub1, res1, scrubbed)
  }
}