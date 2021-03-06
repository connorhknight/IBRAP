% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform.seurat.subclustering.R
\name{perform.seurat.subclustering}
\alias{perform.seurat.subclustering}
\title{Seurat subclustering}
\usage{
perform.seurat.subclustering(
  object,
  assay,
  clust.method,
  column,
  clusters,
  neighbours,
  algorithm,
  res = 0.6,
  ...
)
}
\arguments{
\item{object}{An IBRAP S4 class object}

\item{assay}{Character. Which assay within the object to access}

\item{clust.method}{Character. Which cluster_assignments dataframe to access}

\item{column}{Character. Which column to access within the cluster_assignment dataframe}

\item{clusters}{Which cluster(s) would you like to subcluster}

\item{neighbours}{Character. String indicating which neighbourhood graphs should be used.}

\item{algorithm}{Numerical. Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python. Default = 1 Default = NULL}

\item{res}{Numerical vector. Which resolution to run the clusterign algorithm at, a smaller and larger value identified less and more clusters, respectively. Default = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)}

\item{...}{arguments to be passed to Seurat::FindClusters}

\item{cluster.df.name}{Character. What to call the df contained in clusters. Default = 'seurat}
}
\value{
A new column within the defined cluster_assignment dataframe containing original and new subclusters
}
\description{
Seurat subclustering
}
