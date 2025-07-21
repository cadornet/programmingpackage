#' Compute genes with UMI ≥3 using base matrix and plot violin
#'
#' This function calculates the number of genes with UMI ≥3 per cell and plots a violin plot.
#' It uses colSums(counts >= 3) exactly as in the original script.
#'
#' @param seurat_obj A Seurat object
#' @param feature_name Name of metadata column to store (default: "n_genes_≥3UMI")
#' @param save_path Path to save the violin plot (default: "violinplot_umi3.png")
#'
#' @return A Seurat object with the new metadata column
#' @export
#'
#' @importFrom Seurat GetAssayData VlnPlot
#' @import ggplot2
summarize_expression_umi3_exact <- function(seurat_obj,
                                            feature_name = "n_genes_≥3UMI",
                                            save_path = "violinplot_umi3.png") {
  counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")


  counts_mat <- as.matrix(counts)
  genes_above_3 <- colSums(counts_mat >= 3)

  seurat_obj[[feature_name]] <- genes_above_3

  p <- Seurat::VlnPlot(seurat_obj, features = feature_name) +
    ggplot2::ggtitle("Number of genes with expression ≥3 UMIs")

  ggplot2::ggsave(save_path, plot = p, width = 6, height = 4)

  return(seurat_obj)
}

