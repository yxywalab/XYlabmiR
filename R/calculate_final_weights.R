#' Title
#'@description
#'To calculate genes weights based on protein and miRNA relation.
#'
#' @param gene_list The gene list for analysis.
#' @param protein_file A total file of protein relationship.
#' @param mirna_file The miRNA-mRNA relationship file to obtain vital miRNA.
#' @param use_full_data The selection to determine whether to use all proteins.
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
#' @return Final_weights
#' @export
#'
#' @examples
#' library(data.table)
#' library(igraph)
#' library(biomaRt)
#' library(ggraph)
#' library(R.utils)# #
#' protein_file <- system.file("extdata", "9606.protein.links.v12.0.txt.gz", package = "XYlabmiR")
#' mirna_file <- system.file("extdata", "miRTarget.zip", package = "XYlabmiR")#
#' #
#' gene_list <- c("EGFR", "EPCAM", "VIM")#
#' calculate_final_weights(gene_list, protein_file, mirna_file, use_full_data = TRUE)


calculate_final_weights <- function(gene_list, protein_file, mirna_file, use_full_data = TRUE) {
  .datatable.aware <- TRUE
  initial_weights <- calculate_node_importance(gene_list, protein_file, use_full_data)
  mirna_counts <- get_mirna_gene_relationship(mirna_file)

  final_weights <- merge(initial_weights, mirna_counts, by.x = "gene", by.y = "mRNA", all.x = TRUE)
  final_weights[, final_weight := centrality * mirna_count]
  final_weights <- final_weights[order(-final_weight)]

  return(final_weights)
}
