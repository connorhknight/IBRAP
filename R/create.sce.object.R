#' Produces SingleCellExperiment object
#'
#' This function inputs an expression matrix into a SingleCellExperiment S4 class object whilst annotating fundemental metadata on raw data.
#' 
#' @param counts An expression MATRIX containing the raw counts of individual samples.
#' @param project.name Here name the individual sample/project.
#' @param min.cells A feature must be present in an x number of cells to be retained.
#' @param metadata A dataframe of additional metadata to add to the SingleCellExperiment object.
#' @export
#' @examples create_sce_object(counts = matrix, project.name = "pancreas_1", min.cells = 3, metadata = new.metadata)
#'

create.sce.object <- function(counts, 
                              project.name, 
                              min.cells, 
                              metadata) {
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  cat(crayon::cyan('SCE object created\n'))
  colData(sce)$original.project <- project.name
  sce <- metadata.generator(object = sce, assay = 'counts', prefix = 'raw_')
  if(!is.null(metadata))  {
    colData(sce) <- cbind(colData(sce), metadata)
    cat(crayon::cyan('Custom metadata appended\n'))
  }
  sce <- sce[!duplicated(rownames(assay(sce, 'counts'))), ]
  cat(crayon::cyan('Completed\n'))
  return(sce)
}
