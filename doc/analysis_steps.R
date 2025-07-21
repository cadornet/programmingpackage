## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

library(Seurat)
library(programmingpackage)
library(SingleR)
library(celldex)
library(SummarizedExperiment)

## ----step0--------------------------------------------------------------------
data_dir <- system.file("extdata/filtered_feature_bc_matrix", package = "programmingpackage")
gtf_path <- system.file("extdata/Homo_sapiens.GRCh38.111.gtf", package = "programmingpackage")

sc_data <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = sc_data)

## ----step1--------------------------------------------------------------------
seurat_obj <- select_protein_coding_genes(seurat_obj, gtf_path)

## ----step2--------------------------------------------------------------------
seurat_obj <- summarize_expression_umi3_exact(seurat_obj)

## ----show-violinplot, echo = FALSE--------------------------------------------
knitr::include_graphics("images/violinplot_umi3.png")

## ----step3--------------------------------------------------------------------
seurat_obj <- filter_unwanted_genes(seurat_obj)
removed_table <- read.table(
  system.file("extdata/my_removed_genes.txt", package = "programmingpackage"),
  header = TRUE, sep = "\t"
)
knitr::kable(removed_table, caption = "Number of removed genes by category")

## ----step4--------------------------------------------------------------------
result_pca <- run_pca_variance_plot(seurat_obj)
seurat_obj <- result_pca$seurat_obj

## ----show-pca, echo = FALSE---------------------------------------------------
knitr::include_graphics("images/pca_20PC.png")

## ----step5--------------------------------------------------------------------
result_umap <- run_umap_plot(seurat_obj, dims = 1:10)
seurat_obj <- result_umap$seurat_obj

## ----show-umap, echo = FALSE--------------------------------------------------
knitr::include_graphics("images/umap_10PC.png")

## ----step6--------------------------------------------------------------------
result_cluster <- run_clustering(seurat_obj, dims = 1:10, resolution = 0.5)
seurat_obj <- result_cluster$seurat_obj

## ----show-clustering, echo = FALSE--------------------------------------------
knitr::include_graphics("images/umap_clusters_res05.png")

## ----step7--------------------------------------------------------------------
result_annot <- annotate_cells_singleR(seurat_obj)
seurat_obj <- result_annot$seurat_obj

## ----show-annotation, echo = FALSE--------------------------------------------
knitr::include_graphics("images/umap_celltypes_SingleR.png")

