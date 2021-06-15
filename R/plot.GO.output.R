#' @name plot.GO.output
#' @aliases plot.GO.output
#' 
#' @title Create an IBRAP class object 
#'
#' @description Creates and produces project metadata into an IBRAP S4 class object
#'
#' @usage createIBRAPobject(counts = counts, original.project = 'project_1', method.name = 'RAW', meta.data = df, min.cells = 3, min.features = 200)
#' 
#' @return IBRAP S4 class object containing raw counts and metadata
#'
#' @export

plot.GO.output <- function(result) {
  
  result$rank <- -log10(result$rank)
  result$cluster <- as.character(result$cluster)
  
  ggplot2::ggplot(data = result, ggplot2::aes(y = Term, x = cluster, color = rank)) + 
    ggplot2::geom_point(ggplot2::aes(size = rank)) + ggplot2::scale_size_continuous() + 
    ggplot2::theme_bw() + ggplot2::guides(fill=ggplot2::guide_legend(title="-log10(rank)"))
  
}
