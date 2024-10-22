#' Title
#'@description
#'To obtain aim genes.
#'
#' @param gene_list The gene list for analysis.
#' @param protein_file A total file of protein relationship.
#' @param use_full_data The miRNA-mRNA relationship file to obtain vital miRNA.
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
#' @importFrom data.table .N
#' @return A result containing genes importance.
#' @export
#'
#' @examples
#' library(data.table)
#' library(igraph)
#' library(biomaRt)
#' library(ggraph)
#' library(R.utils)
#' #
#' protein_file <- system.file("extdata", "9606.protein.links.v12.0.txt.gz", package = "XYlabmiR")
#'
#' #
#' gene_list <- c("EGFR", "EPCAM", "VIM")
#'
#' calculate_node_importance(gene_list, protein_file, use_full_data = TRUE)
calculate_node_importance <- function(gene_list, protein_file, use_full_data = TRUE) {
  #protein_data = data.table::fread(protein_file, header = TRUE)
  .datatable.aware <- TRUE
  protein_data <- as.data.table(data.table::fread(protein_file, header = TRUE))

  protein_data[, protein1 := sub("9606.", "", protein1)]
  protein_data[, protein2 := sub("9606.", "", protein2)]

  mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  gene_mapping <- getBM(attributes = c('ensembl_peptide_id', 'external_gene_name'),
                        filters = 'ensembl_peptide_id',
                        values = unique(c(protein_data$protein1, protein_data$protein2)),
                        mart = mart)

  protein_data <- merge(protein_data, gene_mapping, by.x = "protein1", by.y = "ensembl_peptide_id", all.x = TRUE)
  setnames(protein_data, "external_gene_name", "gene1")
  protein_data <- merge(protein_data, gene_mapping, by.x = "protein2", by.y = "ensembl_peptide_id", all.x = TRUE)
  setnames(protein_data, "external_gene_name", "gene2")

  if (!use_full_data) {
    protein_data <- protein_data[gene1 %in% gene_list & gene2 %in% gene_list]
  }

  graph <- graph_from_data_frame(d = protein_data[, .(gene1, gene2, combined_score)], directed = FALSE)
  weighted_centrality <- graph.strength(graph, mode = "all", loops = FALSE)

  centrality_df <- data.table(gene = names(weighted_centrality), centrality = weighted_centrality)
  result <- centrality_df[gene %in% gene_list, ]

  return(result)
}
