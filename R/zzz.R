# Set up environment

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to IBRAP")
}

.onLoad <- function(libname, pkgname) {
  reticulate::configure_environment(pkgname)
  Sys.setenv(R_MAX_VSIZE=250Gb)
}