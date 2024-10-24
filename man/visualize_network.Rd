% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize_network.R
\name{visualize_network}
\alias{visualize_network}
\title{Visualize Network}
\usage{
visualize_network(final_weights, protein_file, gene_list, use_full_data = TRUE)
}
\arguments{
\item{final_weights}{The result of function calculate-final-weights}

\item{protein_file}{A total file of protein relationship.}

\item{gene_list}{The gene list for analysis.}

\item{use_full_data}{The selection to determine whether to use all proteins.}
}
\value{
A graph of gene network.
}
\description{
This function is used to create a network of selected genes.
}
\details{
Title
}
\examples{
library(data.table)
library(igraph)
library(biomaRt)
library(ggraph)
library(R.utils)
library(ggplot2)
#
protein_file <- system.file("extdata", "9606.protein.links.v12.0.txt.gz", package = "XYlabmiR")
mirna_file <- system.file("extdata", "miRTarget.zip", package = "XYlabmiR")
final_weights <- system.file("extdata", "final_weights.zip", package = "XYlabmiR")
#
gene_list <- c("EGFR", "EPCAM", "VIM")
visualize_network(final_weights, protein_file, gene_list, use_full_data = FALSE)
}
