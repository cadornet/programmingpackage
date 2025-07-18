#' Run clustering and plot clusters on UMAP
#'
#' This function runs clustering on a Seurat object using the specified PCA dimensions and resolution,
#' plots the resulting clusters on a UMAP embedding, and saves the plot.
#'
#' @param seurat_obj A Seurat object with PCA and UMAP already computed
#' @param dims PCA dimensions to use for clustering (default: 1:10)
#' @param resolution Clustering resolution parameter (default: 0.5)
#' @param save_path Path to save the UMAP cluster plot (default: "umap_clusters_res05.png")
#'
#' @return A list containing:
#' \describe{
#'   \item{seurat_obj}{The Seurat object with clustering information}
#'   \item{plot}{The UMAP plot with clusters}
#'   \item{cluster_table}{A table with number of cells per cluster}
#' }
#' @export
#'
#' @importFrom Seurat FindNeighbors FindClusters DimPlot Idents
#' @import ggplot2
run_clustering <- function(seurat_obj,
                           dims = 1:10,
                           resolution = 0.5,
                           save_path = "umap_clusters_res05.png") {
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution)
  
  cluster_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggplot2::ggtitle(sprintf("Cellular clusters (resolution %.2f)", resolution))
  
  ggplot2::ggsave(save_path, plot = cluster_plot, width = 7, height = 5)
  
  cluster_table <- table(Seurat::Idents(seurat_obj))
  
  return(list(
    seurat_obj = seurat_obj,
    plot = cluster_plot,
    cluster_table = cluster_table
  ))
}
