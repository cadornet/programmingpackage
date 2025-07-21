#' Annotate cell types using SingleR and HumanPrimaryCellAtlas
#'
#' This function annotates cells in a Seurat object using SingleR and the Human Primary Cell Atlas reference.
#' It stores the predicted labels in the metadata, plots a UMAP with labels, and saves it.
#'
#' @param seurat_obj A Seurat object
#' @param save_path File path to save the UMAP plot (default: "umap_celltypes_SingleR.png")
#'
#' @return A list containing:
#' \describe{
#'   \item{seurat_obj}{The Seurat object with SingleR labels added}
#'   \item{plot}{The UMAP plot grouped by predicted cell type}
#'   \item{confusion_table}{Contingency table of clusters vs. predicted labels}
#' }
#' @export
#'
#' @importFrom Seurat GetAssayData DimPlot Idents
#' @importFrom celldex HumanPrimaryCellAtlasData
#' @importFrom SingleR SingleR
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @import ggplot2
annotate_cells_singleR <- function(seurat_obj,
                                   save_path = "umap_celltypes_SingleR.png") {
  expr_matrix <- Seurat::GetAssayData(seurat_obj, slot = "data")
  
e
  hpca <- celldex::HumanPrimaryCellAtlasData()
  
 
  singleR_results <- SingleR::SingleR(test = expr_matrix,
                                      ref = hpca,
                                      labels = hpca$label.main)
  
  seurat_obj$SingleR_label <- singleR_results$labels
  
  # Plot
  umap_plot <- Seurat::DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE) +
    ggplot2::ggtitle("Cell type annotations (SingleR)")
  
  ggplot2::ggsave(save_path, plot = umap_plot, width = 8, height = 6)
  

  confusion <- table(Cluster = Seurat::Idents(seurat_obj),
                     CellType = seurat_obj$SingleR_label)
  
  return(list(
    seurat_obj = seurat_obj,
    plot = umap_plot,
    confusion_table = confusion
  ))
}
