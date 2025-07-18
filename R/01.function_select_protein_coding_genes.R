#' Select protein-coding genes from a Seurat object using a GTF annotation
#'
#' This function filters a Seurat object to retain only protein-coding genes based on a GTF annotation file.
#'
#' @param seurat_obj A Seurat object created from 10X data
#' @param gtf_path Path to the GTF annotation file (e.g., from Ensembl)
#' @return A filtered Seurat object containing only protein-coding genes
#' @export
#'
#' @importFrom rtracklayer import
#' @import Seurat
select_protein_coding_genes <- function(seurat_obj, gtf_path) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required but not installed.")
  }

  genes_gtf <- rtracklayer::import(gtf_path, format = "gtf", feature.type = "gene")
  protein_coding_genes <- genes_gtf[genes_gtf$gene_biotype == "protein_coding"]
  keep_gene_names <- protein_coding_genes$gene_name

  filtered_seurat <- subset(seurat_obj, features = keep_gene_names)
  return(filtered_seurat)
}
