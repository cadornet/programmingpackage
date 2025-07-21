#' Filter ribosomal, mitochondrial and pseudogenes from a Seurat object
#'
#' This function removes ribosomal genes, ribosomal pseudogenes, and mitochondrial genes from a Seurat object
#' based on gene name patterns (e.g. "RPS", "RPL", "MT-").
#'
#' @param seurat_obj A Seurat object
#' @param out_list_file Path to save the list of removed gene names (default: "removed_genes_list.txt")
#' @param out_summary_file Path to save a summary table of removed genes by category (default: "removed_genes_summary.txt")
#'
#' @return A filtered Seurat object
#' @export
#'
#' @import Seurat
filter_unwanted_genes <- function(seurat_obj,
                                  out_list_file = "removed_genes_list.txt",
                                  out_summary_file = "removed_genes_summary.txt") {
  gene_names <- rownames(seurat_obj)

  ribosomal_genes <- grep("^RPS|^RPL", gene_names, value = TRUE)
  mitochondrial_genes <- grep("^MT-", gene_names, value = TRUE)
  pseudogenes <- grep("RPS.*P|RPL.*P", gene_names, value = TRUE)

  removed_genes <- data.frame(
    Categories = c("Ribosomal", "Mitochondrial", "Pseudogenes"),
    N_removed_genes = c(
      length(ribosomal_genes),
      length(mitochondrial_genes),
      length(pseudogenes)
    )
  )

  genes_to_remove <- unique(c(ribosomal_genes, mitochondrial_genes, pseudogenes))
  filtered_seurat <- subset(seurat_obj, features = setdiff(gene_names, genes_to_remove))

  
  write.table(genes_to_remove, out_list_file, row.names = FALSE, quote = FALSE, col.names = FALSE)
  write.table(removed_genes, out_summary_file, row.names = FALSE, quote = FALSE, sep = "\t")

  print(removed_genes)

  return(filtered_seurat)
}
