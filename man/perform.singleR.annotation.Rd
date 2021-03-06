% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform.singleR.annotation.R
\name{perform.singleR.annotation}
\alias{perform.singleR.annotation}
\title{Performs automated cell annotation on query datasets using reference data}
\usage{
perform.singleR.annotation(
  object,
  assay = "RAW",
  slot = "counts",
  ref,
  log.transform = TRUE,
  tpm.transform = FALSE,
  ref.labels,
  column.suffix = "1",
  ...
)
}
\arguments{
\item{object}{IBRAP S4 class object}

\item{assay}{Character. String containing indicating which assay to use}

\item{slot}{Character. String indicating which slot within the assay should be sourced}

\item{ref}{AnyMatrix. A matrix of the reference datasets, if data is end-bias then it should be log normalised, if it is full-length then it requires tpm normalisation. Both can be completed within this function.}

\item{log.transform}{Boolean. Should the reference data be log transformed. Default = TRUE}

\item{tpm.transform}{Boolean. Should the reference data be tpm normalised. Default = FALSE}

\item{ref.labels}{Vector. The cluster assignments for the reference data. Default = NULL}

\item{column.suffix}{Character. A suffix to append the end of the new metadata columns if this functiuons is to be used multiple times. Default = '1'}

\item{...}{arguments to be passed to singleR::SingleR}
}
\value{
Produces a new 'methods' assay containing normalised, scaled and HVGs.
}
\description{
SingleR iterates through singular cells and iterates through probabilitiy comparisons to identify which cell type the query cell is likely to be. If a probably cell type cannot be discovered then
}
