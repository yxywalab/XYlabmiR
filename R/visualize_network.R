#' Title
#' @title Visualize Network
#'
#' @description This function is used to create a network of selected genes.
#'
#' @param final_weights The result of function calculate-final-weights
#' @param protein_file A total file of protein relationship.
#' @param gene_list The gene list for analysis.
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
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_void
#' @importFrom data.table .N
#'
#' @return A graph of gene network.
#' @export
#'
#' @examples
#' library(data.table)
#' library(igraph)
#' library(biomaRt)
#' library(ggraph)
#' library(R.utils)
#' library(ggplot2)
#' #
#' protein_file <- system.file("extdata", "9606.protein.links.v12.0.txt.gz", package = "XYlabmiR")
#' mirna_file <- system.file("extdata", "miRTarget.zip", package = "XYlabmiR")
#' final_weights <- system.file("extdata", "final_weights.zip", package = "XYlabmiR")
#' #
#' gene_list <- c("EGFR", "EPCAM", "VIM")
#' visualize_network(final_weights, protein_file, gene_list, use_full_data = FALSE)
visualize_network <- function(final_weights, protein_file, gene_list, use_full_data = TRUE) {
  .datatable.aware <- TRUE
  protein_data <- fread(protein_file, header = TRUE)


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

  final_weights <- fread(final_weights, header = TRUE)

  V(graph)$size <- final_weights$final_weight
  V(graph)$label <- final_weights$gene


  mirna_counts <- final_weights$mirna_count
  names(mirna_counts) <- final_weights$gene
  V(graph)$mirna_count <- ifelse(V(graph)$label %in% names(mirna_counts),
                                 mirna_counts[match(V(graph)$label, names(mirna_counts))],
                                 0)


  V(graph)$size <- (V(graph)$size^(1/3) - min(V(graph)$size^(1/3))) /
    (max(V(graph)$size^(1/3)) - min(V(graph)$size^(1/3))) * 5


  E(graph)$width <- E(graph)$combined_score / max(E(graph)$combined_score, na.rm = TRUE) * 5


  ggraph(graph, layout = "fr") +
    geom_edge_link(aes(width = width), alpha = 0.5, color = "gray") +
    geom_node_point(aes(size = size), color = "steelblue") +
    geom_node_text(aes(label = label, filter = size >= 0), vjust = 1.5) +
    geom_node_text(aes(label = mirna_count), vjust = -1, color = "red") +
    scale_size(range = c(1, 10)) +
    labs(title = "Gene Interaction Network",
         subtitle = "Node size indicates gene weight (normalized), edge width indicates combined score") +
    theme_void()
}
