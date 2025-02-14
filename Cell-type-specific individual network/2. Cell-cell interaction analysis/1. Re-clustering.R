library(Seurat)
library(celldex)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)

####Re-clustering of T cells/epithelial cells
##load data
pbmc <- readRDS("cancer_batch_cluster.rds") 

##Select the T cells or epithelial cells of all samples
#T cell
seurat_obj <- subset(x = pbmc, idents = 0 )
#epithelial cell
#seurat_obj <- subset(x = pbmc,idents = c(2,4,6,10))
# ident == "fill in the cell cluster number annotated as T cells or epithelial cells"

pbmc <- seurat_obj

## re-clustering
library(tidyverse)
combined.integrated <- FindVariableFeatures(pbmc, verbose = FALSE, assay = "RNA") %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = 'pca', dims = 1:20)
combined.integrated <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
combined.integrated <- FindClusters(combined.integrated, resolution = 0.5)

pbmc <- combined.integrated 

# Select an appropriate resolution
for (i in c(0.1,0.2,0.3)){
  seurat_obj = FindClusters(pbmc,resolution = i,algorithm =1)
}

library(clustree)
clustree(seurat_obj@meta.data,prefix= "RNA_snn_res.")

pbmc <- seurat_obj
saveRDS(pbmc,"Tcell.rds")
#saveRDS(pbmc,"epithelial.rds")











