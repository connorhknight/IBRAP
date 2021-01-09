#' Produces a 3D scatterplot displaying specific metadata information.
#'
#' This function produces a plot with ggplots2, where three metadata variables can be visualised simultaneously.
#' The y-axis and x-axis will display two metadata variables, and the thrid dimension is introduced with a colour gradient.
#'
#' Remember, you can view your available  metadata by specifying colData(sce_object)
#'
#' @import ggplot2
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param x x-axis
#' @param y y-axis
#' @param z colour gradient
#' @export

scQC_correlate <- function(object, x, y, z){
  metadata <- as.data.frame(colData(object))
  tmp <- data.frame(x=as.numeric(metadata[ , grepl(x, names(metadata)) ]))
  tmp$y <- as.numeric(metadata[ , grepl(y, names(metadata)) ])
  tmp$z <- as.numeric(metadata[ , grepl(z, names(metadata)) ])
  tmp$original.project <- object$original.project
  rownames(tmp) <- colnames(object)
  plot <- ggplot(tmp, aes(x=x, y=y, color=z)) +
    geom_point() +
    scale_colour_gradient(low = "blue", high = "red") +
    stat_smooth(method=lm, formula = y ~ x + z) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 3000) +
    geom_hline(yintercept = 1000) +
    facet_wrap(~original.project) +
    labs(x=x, y=x, color=z)
  print(plot)
}
