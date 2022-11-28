# Set up environment

.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("Welcome to IBRAP")
  
}

.onLoad <- function(libname, pkgname) {
  
  if(!suppressWarnings(suppressPackageStartupMessages(require('celltalker', quietly=TRUE,character.only=TRUE)))){
    remotes::install_github("arc85/celltalker")
    suppressPackageStartupMessages(library('celltalker',character.only=TRUE))
  }
  
}
