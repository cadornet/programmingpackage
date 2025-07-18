#' Run PCA on a Seurat object and plot explained variance
#'
#' This function normalizes, scales, and runs PCA on a Seurat object. It then computes
#' the percentage of variance explained by the first 20 principal components and
#' generates a barplot, which is saved as a PNG file.
#'
#' @param seurat_obj A Seurat object
#' @param filename Name of the PNG file to save the plot (default: "pca_20PC.png")
#' @return A list containing:
#' \describe{
#'   \item{seurat_obj}{The Seurat object with PCA results}
#'   \item{pca_var_perc}{Vector with percentage variance explained per PC}
#'   \item{pca_var_cumsum}{Vector with cumulative variance explained}
#' }
#' @export
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA VariableFeatures
run_pca_variance_plot <- function(seurat_obj, filename = "pca_20PC.png") {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  
  # Varianza spiegata
  pca_var <- seurat_obj[["pca"]]@stdev^2
  pca_var_perc <- pca_var / sum(pca_var) * 100
  pca_var_cumsum <- cumsum(pca_var_perc)
  
  # Plot
  png(filename = filename, width = 800, height = 600)
  barplot(pca_var_perc[1:20],
          names.arg = 1:20,
          xlab = "Principal Component",
          ylab = "Percent Variance Explained",
          main = "PCA - Explained Variance (first 20 PCs)")
  dev.off()
  
  return(list(
    seurat_obj = seurat_obj,
    pca_var_perc = pca_var_perc,
    pca_var_cumsum = pca_var_cumsum
  ))
}
