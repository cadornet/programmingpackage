# Programming approaches for Bioinformatics's Exam – July 2025

This repository contains the instructions to complete a single-cell RNA-seq analysis pipeline using **Seurat** and **SingleR**.  
It includes functions for filtering protein-coding genes, excluding unwanted features, computing PCA/UMAP, clustering, and annotating cell types with external references.  
The project provides a complete R package with modular functions needed to perform the full analysis described in [this vignette](https://github.com/cadornet/programmingpackage/blob/main/docs/analysis_steps.html).

## Installation

Clone the repository and load the package locally using `devtools`:
```r
# install.packages("devtools")
devtools::install_github("cadornet/programmingpackage")
```
## **Features**

- Protein-coding gene filtering via GTF file
- UMI ≥3 expression summarization
- Filtering of ribosomal, mitochondrial and pseudogenes
- PCA and variance plot
- UMAP embedding and visualization
- Clustering with SNN and resolution control
- Cell type annotation via SingleR and Human Primary Cell Atlas


## Example Workflow


```r
library(programmingpackage)
library(Seurat) 

# Load data --> load the raw single-cell dataset
sc_data <- Read10X(data.dir = "your_path/filtered_feature_bc_matrix/")
seurat_obj <- CreateSeuratObject(counts = sc_data)

# Filter protein-coding genes --> Using the Ensembl GTF annotation file, filter the dataset to retain only protein-coding genes. The filtering ensures that downstream analysis focuses only on protein-coding gene expression profiles.
gtf_path <- "your_path/Homo_sapiens.GRCh38.111.gtf"
seurat_obj <- select_protein_coding_genes(seurat_obj, gtf_path)

# Summarize UMI expression --> To assess the richness of gene expression per cell, calculate the number of genes with expression ≥3 UMIs for each cell.
seurat_obj <- summarize_expression_umi3_exact(seurat_obj)

# Filter unwanted features --> Remove unwanted gene categories from the dataset. The filtering ensures that low-complexity or confounding gene types are excluded from downstream analysis.
seurat_obj <- filter_unwanted_genes(seurat_obj)

# Run PCA and UMAP --> Performe dimensionality reduction using PCA to capture the major sources of variation in the dataset. Then applied UMAP for visualization of the dataset in a low-dimensional space to provide an intuitive representation of cell relationships.
pca_result <- run_pca_variance_plot(seurat_obj)
seurat_obj <- pca_result$seurat_obj

umap_result <- run_umap_plot(seurat_obj)
seurat_obj <- umap_result$seurat_obj

# Clustering and annotation --> Performe clustering on the dataset to identify groups of transcriptionally similar cells.
cluster_result <- run_clustering(seurat_obj)
annot_result <- annotate_cells_singleR(cluster_result$seurat_obj)
```
## Vignette
A complete step-by-step explanation of the analysis is available [here](https://github.com/cadornet/programmingpackage/blob/main/docs/analysis_steps.html).

## License
MIT © 2025 Chiara Adornetto

