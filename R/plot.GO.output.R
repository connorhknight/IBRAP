#' @name plot.GO.output
#' @aliases plot.GO.output
#' 
#' @title Plot GO enrichment output
#'
#' @description Creates a dot plot of GO terms
#'
#' @param result The result file from the perform.GO.enrichment
#' 
#' @return IBRAP S4 class object containing raw counts and metadata
#'
#' @export plot.GO.output

plot.GO.output <- function(result) {
  
  result$rank <- -log10(result$rank)
  result$cluster <- as.character(result$cluster)
  
  ggplot2::ggplot(data = result, ggplot2::aes(y = Term, x = cluster, color = rank)) + 
    ggplot2::geom_point(ggplot2::aes(size = rank)) + ggplot2::scale_size_continuous() + 
    ggplot2::theme_bw() + ggplot2::guides(color=ggplot2::guide_legend(title="-log10(rank)"), size=ggplot2::guide_legend(title="-log10(rank)"))
  
}
