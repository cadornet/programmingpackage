#' Run UMAP and generate a DimPlot
#'
#' This function runs UMAP on the provided Seurat object using the specified PCA dimensions,
#' generates a UMAP plot and saves it as a PNG file.
#'
#' @param seurat_obj A Seurat object that has already undergone PCA
#' @param dims Numeric vector of PCA dimensions to use (default: 1:10)
#' @param save_path File path to save the UMAP plot (default: "umap_10PC.png")
#'
#' @return A list containing:
#' \describe{
#'   \item{seurat_obj}{The Seurat object with UMAP embedding}
#'   \item{plot}{The ggplot2 object of the UMAP plot}
#' }
#' @export
#'
#' @importFrom Seurat RunUMAP DimPlot
#' @import ggplot2
run_umap_plot <- function(seurat_obj,
                          dims = 1:10,
                          save_path = "umap_10PC.png") {
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = dims)
  
  umap_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap") +
    ggplot2::ggtitle(sprintf("UMAP based on first %d PCA components", max(dims)))
  
  ggplot2::ggsave(save_path, plot = umap_plot, width = 7, height = 5)
  
  return(list(
    seurat_obj = seurat_obj,
    plot = umap_plot
  ))
}
