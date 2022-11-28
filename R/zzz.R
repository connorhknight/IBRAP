# Set up environment

.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("Welcome to IBRAP")
  
}

.onLoad <- function(libname, pkgname) {
  
  for(pkg in pkglist)
    if(!suppressWarnings(suppressPackageStartupMessages(require('celltalker', quietly=TRUE,character.only=TRUE)))){
      devtools::install_github("arc85/celltalker")
      suppressPackageStartupMessages(library('celltalker',character.only=TRUE))
    }
  
}
