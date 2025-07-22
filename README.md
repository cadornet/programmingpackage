# Programming approaches for Bioinformatics's Exam – July 2025

This repository contains the instructions to complete a single-cell RNA-seq analysis pipeline using **Seurat** and **SingleR**.  
It includes functions for filtering protein-coding genes, excluding unwanted features, computing PCA/UMAP, clustering, and annotating cell types with external references.  
The project provides a complete R package with modular functions needed to perform the full analysis described in this [Rmarkdown](https://github.com/cadornet/programmingpackage/blob/main/correct_rmd_scRNAseq_analysys_for_programming.Rmd)


## Installation

Clone the repository and load the package locally using `devtools`:
```r
# install.packages("devtools")
devtools::install_github("cadornet/programmingpackage")
```
## **Features**

- Protein-coding gene filtering via GTF file
- [UMI ≥3 expression summarization](https://github.com/cadornet/programmingpackage/blob/main/results/violinplot_umi3.png?raw=true)
- [Filtering of ribosomal, mitochondrial and pseudogenes](https://github.com/cadornet/programmingpackage/blob/main/results/my_summary.txt?raw=true)
- [PCA and variance plot](https://github.com/cadornet/programmingpackage/blob/main/results/pca_20PC.png?raw=true)
- [UMAP embedding and visualization](https://github.com/cadornet/programmingpackage/blob/main/results/umap_10PC.png?raw=true)
- [Clustering](https://github.com/cadornet/programmingpackage/blob/main/results/umap_clusters_res05.png?raw=true)
- [Cell type annotation via SingleR and Human Primary Cell Atlas](https://github.com/cadornet/programmingpackage/blob/main/results/umap_celltypes_SingleR.png?raw=true)


## Example Workflow
suggestion --> A deeper and more detailed explanation of each function including can be found [here](https://github.com/cadornet/programmingpackage/blob/main/docs/detailed_explanation.Rmd)

```r
library(programmingpackage)
library(Seurat) 

# Load data --> load the raw single-cell dataset
sc_data <- Read10X(data.dir = "your_path/filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = sc_data)

# Filter protein-coding genes --> This function filters the Seurat object to retain only protein-coding genes, based on 
annotations from a provided GTF file.
  select_protein_coding_genes <- function(seurat_obj, gtf_path) {
  genes_gtf <- rtracklayer::import(gtf_path, format = "gtf", feature.type = "gene")
  protein_coding_genes <- genes_gtf[genes_gtf$gene_biotype == "protein_coding"]
  keep_gene_names <- protein_coding_genes$gene_name
  filtered_seurat <- subset(seurat_obj, features = keep_gene_names)
  return(filtered_seurat)
}

# Summarize UMI expression --> This function calculates, for each cell, the number of genes with expression levels of 
at least 3 UMIs, stores the result as a new metadata column, and visualizes the distribution using a violin plot.
  summarize_expression_umi3_exact <- function(seurat_obj,
  feature_name = "n_genes_≥3UMI") {
  counts <- GetAssayData(seurat_obj, layer = "counts")
  counts_mat <- as.matrix(counts)
  seurat_obj[[feature_name]] <- colSums(counts_mat >= 3) 
  VlnPlot(seurat_obj, features = feature_name) + 
    ggtitle("Genes with ≥3 UMI per cell")
}

# Filter unwanted features --> This function removes unwanted genes from the Seurat object, specifically mitochondrial,
 ribosomal, and pseudogenes, based on gene name patterns, to improve downstream analysis quality.
  filter_unwanted_genes <- function(seurat_obj) {
  genes <- rownames(seurat_obj)
  genes_to_remove <- unique(c( 
    grep("^RPS|^RPL", genes, value = TRUE),
    grep("^MT-", genes, value = TRUE),
    grep("RPS.*P|RPL.*P", genes, value = TRUE)
  ))
  subset(seurat_obj, features = setdiff(genes, to_remove))
}

# Run PCA -> This function performs normalization, variable feature selection, scaling, and PCA on the Seurat object. 
It then computes the variance explained by each principal component and visualizes the first 20 components in a barplot.
run_pca_variance_plot <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  pca_var <- seurat_obj[["pca"]]@stdev^2
  pca_perc <- pca_var / sum(pca_var) * 100
  barplot(pca_perc[1:20],
          names.arg = 1:20,
          xlab = "PC",
          ylab = "% Variance explained",
          main = "PCA - explained variance of the first 20 PCs")
}

# UMAP --> This function computes a UMAP embedding using the specified PCA dimensions and generates a UMAP plot with a
dynamic title indicating the number of components used.
  run_umap_plot <- function(seurat_obj, dims = 1:10) {
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  DimPlot(seurat_obj, reduction = "umap") +
    ggtitle(sprintf("UMAP based on first %d PCA components", max(dims)))
}

# Clustering --> The function constructs a nearest-neighbor graph based on the selected PCA dimensions and applies a 
community detection algorithm to identify cell clusters, then it returns a UMAP plot for visualization.
  run_clustering <- function(seurat_obj, dims = 1:10, resolution = 0.5) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  DimPlot(seurat_obj, reduction = "umap", label = TRUE) +
    ggtitle(sprintf("Clusters (resolution %.2f)", resolution))
}

#Annotation --> This function automatically assigns a cell type identity to each cell by comparing its gene expression 
profile with reference profiles from the HPCA. The predicted labels are added to the Seurat object metadata and visualized on a UMAP plot.
annotate_cells_singleR <- function(seurat_obj) {
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  ref <- celldex::HumanPrimaryCellAtlasData()
  singleR_labels <- SingleR(test = expr_matrix,
                            ref = ref,
                            labels = ref$label.main)$labels
  seurat_obj$SingleR_label <- singleR_labels
  DimPlot(seurat_obj, group.by = "SingleR_label", label = TRUE) +
    ggtitle("Cell type annotations (SingleR)")
}

# Infer Tissue of Origin --> This function estimates the tissue of origin of the sample by analyzing the most represented annotated
cell types. It uses regular expressions to detect characteristic patterns within the top three cell types and returns a biologically
 plausible tissue classification.
infer_tissue_origin <- function(seurat_obj) {
  cell_types <- seurat_obj$SingleR_label
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
message("Inferred tissue of origin: ", tissue_guess)
}

```
# The final interpretation is saved in:
[tissue_origin_inference_summary.txt](https://github.com/cadornet/programmingpackage/blob/main/vignettes/images/tissue_origin_inference_summary.txt)


## Vignette
A complete step-by-step explanation of the analysis is available [here](https://cadornet.github.io/programmingpackage/correct_vignette_analysis_steps.html).

## Docker Support
This project also includes a ready-to-use Docker Container to ensure reproducility and compatibility
The **Dockerfile** used to build the container is included in this repository (or you can find it [here](https://github.com/cadornet/programmingpackage/blob/main/Dockerfile?raw=true)) and installs all required system and R dependencies, including:
Seurat
rtracklayer
SingleR
celldex
SummarizedExperiment
ggplot2, and others
### Build the Docker Image
To build the Docker image:
```bash
docker build -t programmingpackage .
```
To start an interactive R session inside the container:
```bash
docker run -it programmingpackage
```
Once inside, load the package and start using its functions:
```r
library(programmingpackage)
?select_protein_coding_genes
```

## Docker Image
To ensure reproducibility, the entire analysis pipeline is also available as a Docker image.
You can pull and run the container from Docker Hub:
```bash
docker pull chiara9/programmingpackage
docker run -it chiara9/programmingpackage
```
This image includes all necessary dependencies, R packages (such as Seurat, SingleR, rtracklayer, etc.), and the custom R package programmingpackage used in this project.
🔗 Docker Hub Link: https://hub.docker.com/r/chiara9/programmingpackage





## License
MIT © 2025 Chiara Adornetto

