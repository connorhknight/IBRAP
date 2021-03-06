% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.feature.score.R
\name{add.feature.score}
\alias{add.feature.score}
\title{Provides scores for a given vector of features}
\usage{
add.cell.cycle(object = obj, assay = 'RAW', slot = 'counts')
}
\arguments{
\item{object}{IBRAP S4 class object}

\item{assay}{A character string containing indicating which assay to use}

\item{slot}{String indicating which slot within the assay should be sourced}

\item{transform}{Boolean. If raw counts are supplied, this must be TRUE to normalise data}

\item{features}{A character vector of genes to be scored}

\item{column.name}{Character naming the column containing the scores in the metadata dataframe}
}
\value{
IBRAP S4 class object containing module scores for each cell in the metadata
}
\description{
Produces a module score per cell for the supplied genes
}
