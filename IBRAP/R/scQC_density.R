#' Produces a density plot displaying specific metadata information.
#'
#' This function produces a plot with ggplots2, where y correlates with n of cells and x with a metadata variable.
#'
#' Remember, you can view your available  metadata by specifying colData(sce_object)
#'
#' @import ggplot2
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param column Which metadata column would you like to
#' @param title What should the title be?
#' @param cutoff Where should a cutoff point be visualised. This is important for comprehending the effect of the filtration stage.
#' @export

scQC_density <- function(object, column, title, cutoff){
  metadata <- as.data.frame(colData(object))
  tmp <- as.data.frame(metadata[ , grepl(column, names(metadata)) ])
  rownames(tmp) <- colnames(object)
  colnames(tmp) <- 'selected_feature'
  tmp$original.project <- metadata$original.project
  plot <- ggplot(tmp, aes(color = original.project, x = selected_feature, fill = original.project)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    theme_classic() +
    ylab("Cell density") +
    xlab(column) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", face = "bold", size = (15))) +
    geom_vline(xintercept = cutoff)
  print(plot)
}
