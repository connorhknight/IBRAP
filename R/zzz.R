# Set up environment

.onAttach <- function(libname, pkgname) {
  
  packageStartupMessage("Welcome to IBRAP")
  
}

.onLoad <- function(libname, pkgname) {
  
  library(Matrix)
  
}