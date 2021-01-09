#' Produces SingleCellExperiment object
#'
#' This function inputs an expression matrix into a singlecellexperiment object whilst calculating total RNA per cell and total feature appearance metadata.
#'
#' @import SingleCellExperiment
#' @param counts An expression MATRIX containing the raw counts of individual samples.
#' @param project.name Here name the individual sample/project.
#' @export
#' @examples create_sce_object(counts = counts.matrix, project.name = "sample_1")
#'

create_sce_object <- function(counts, project.name, min.cells) {
  matrix <- matrix[names(which(rowSums(matrix) > as.numeric(min.cells))),]
  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$original.project <- project.name
  colData(sce)$total.counts <- colSums(assay(sce, 'counts'))
  colData(sce)$total.features <- colSums(assay(sce, 'counts') != 0)
  sce <- sce[!duplicated(rownames(assay(sce, 'counts'))), ]
  return(sce)
}
