# Set up environment

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to IBRAP")
}

.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
}