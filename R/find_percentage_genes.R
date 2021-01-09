#' Calculates the total percentage of genes matching a distinct character pattern.
#'
#' This function calculates the percentage contribution of all transcripts associated with a specified regex pattern
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param pattern A regex pattern to identify genes
#' @param column.name Please name the column containing this feature information.
#' @export

find_percentage_genes <- function(object, pattern='^MT', column.name) {
  temp <- data.frame(tmp=colSums(assay(object, 'counts')[grep(as.character(pattern),rownames(assay(object, 'counts'))),]) / colData(object)$total.counts * 100)
  colnames(temp) <- c(column.name)
  colData(object) <- cbind(colData(object), temp)
  return(object)
}
