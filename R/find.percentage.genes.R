#' Calculates the total percentage of genes matching a distinct character pattern.
#'
#' This function calculates the percentage contribution of all transcripts associated with a specified regex pattern
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param pattern A regex pattern to identify genes
#' @param which.assay Which assay to find percentage genes from
#' @param prefix which prefix to add to metadata column name
#' @export

find.percentage.genes <- function(object, 
                                  pattern='^MT-', 
                                  which.assay='counts', 
                                  prefix) {
  temp <- data.frame(tmp=colSums(assay(object, which.assay)[grep(pattern = pattern, x = rownames(assay(object, which.assay))),]) / colData(object)$total.counts * 100)
  colnames(temp) <- c(column.name)
  if(column.name %in% colnames(colData(object))) {
    colData(object) <- colData(object)[,colnames(colData(object)) != column.name]
  }
  colData(object) <- cbind(colData(object), temp)
  return(object)
}
