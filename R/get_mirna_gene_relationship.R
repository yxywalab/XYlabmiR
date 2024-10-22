#' Title
#'
#' @param mirna_file The miRNA-mRNA relationship file to obtain vital miRNAs.
#'
#' @import R.utils
#' @import biomaRt
#' @import ggraph
#' @import igraph
#' @import data.table
#' @importFrom data.table :=
#' @importFrom data.table as.data.table
#' @importFrom biomaRt useEnsembl
#' @importFrom biomaRt getBM
#' @importFrom igraph graph_from_data_frame
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @return The result of important miRNAs counts.
#' @export
#'
#' @examples
#' library(data.table)
#' library(igraph)
#' library(biomaRt)
#' library(ggraph)
#' library(R.utils)
#' #
#' mirna_file <- system.file("extdata", "miRTarget.zip", package = "XYlabmiR")
#'
#' #
#' gene_list <- c("EGFR", "EPCAM", "VIM")
#' get_mirna_gene_relationship(mirna_file)
get_mirna_gene_relationship <- function(mirna_file) {
  .datatable.aware <- TRUE
  mirna_data <- fread(mirna_file, header = TRUE)
  mirna_count <- mirna_data[, .N, by = mRNA]
  setnames(mirna_count, "N", "mirna_count")
  return(mirna_count)
}
