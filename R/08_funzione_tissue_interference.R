#' Infer the tissue of origin based on SingleR annotations
#'
#' This function analyzes the most abundant cell types predicted by SingleR
#' and proposes a tissue of origin based on common marker patterns.
#' The result is printed and saved to a text file.
#'
#' @param seurat_obj A Seurat object with SingleR annotations (SingleR_label)
#' @param output_file Path to save the output text file (default: "tissue_origin_inference_summary.txt")
#'
#' @return A list containing:
#' \describe{
#'   \item{tissue_guess}{The inferred tissue of origin}
#'   \item{celltype_counts}{The table of predicted cell types}
#'   \item{text}{The inference summary as character}
#' }
#' @export
infer_tissue_origin <- function(seurat_obj, output_file = "tissue_origin_inference_summary.txt") {
  # Extract cell types and cluster IDs
  cell_types <- seurat_obj$SingleR_label
  cluster_ids <- Idents(seurat_obj)
  

  celltype_counts <- sort(table(cell_types), decreasing = TRUE)
  top_labels <- names(celltype_counts)[1:3]
  

  tissue_guess <- "Unknown"
  if (any(grepl("T cell|B cell|NK cell", top_labels))) {
    tissue_guess <- "Blood or lymphoid tissue"
  } else if (any(grepl("Epithelial", top_labels))) {
    tissue_guess <- "Epithelial tissue (e.g. gut, skin, lung)"
  } else if (any(grepl("Neural", top_labels))) {
    tissue_guess <- "Neural tissue"
  } else if (any(grepl("Fibroblast|Endothelial", top_labels))) {
    tissue_guess <- "Connective or stromal tissue"
  }
  

  output_text <- paste0(
    "### Tissue Origin Inference\n\n",
    "Top predicted cell types:\n",
    paste0(" - ", names(head(celltype_counts, 5)), ": ", head(celltype_counts, 5), collapse = "\n"),
    "\n\nInferred tissue of origin:\n",
    "**", tissue_guess, "**\n\n",
    "This hypothesis is based on the frequency of annotated cell types and the cluster structure observed in the UMAP."
  )
  

  writeLines(output_text, con = output_file)
  cat("✅ Tissue inference saved to", output_file, "\n")
  

  return(list(
    tissue_guess = tissue_guess,
    celltype_counts = celltype_counts,
    text = output_text
  ))
}

