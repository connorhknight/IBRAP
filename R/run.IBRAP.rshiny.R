#' @name run.IBRAP.rshiny
#' @aliases run.IBRAP.rshiny
#' 
#' @title Rshiny application initiator
#'
#' @description Initiates the Rshiny application for data interpretation
#'
#' @param ... arguments to pass to shiny::runApp
#'
#' @export run,IBRAP.rshiny

run.IBRAP.rshiny <- function(...) {
  
  shiny::runApp(appDir = system.file('rshiny', 'app.R', package = 'IBRAP'), ...)
  
}
