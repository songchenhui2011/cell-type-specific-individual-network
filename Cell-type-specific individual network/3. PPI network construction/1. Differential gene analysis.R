devtools::install_github('satijalab/seurat-data')
library(Seurat)
library(SeuratData)
library(ggplot2)

#load data
pbmc <- readRDS("ancer_batch_cluster.rds") 

## Select the T cells or epithelial cells of all samples
# T cells
seurat_obj <- subset(x = pbmc, idents = 0 )
# epithelial cells of
#seurat_obj <- subset(x = pbmc, idents = c(2,4,6,10))
#ident == "fill in the cell cluster number annotated as T cells or epithelial cells"

ifnb <- seurat_obj

# Normalize the data
ifnb <- NormalizeData(ifnb)

# Find differentially expressed genes between two samples from the same patient
Idents(ifnb) <- "orig.ident"
monocyte.de.markers <- FindMarkers(ifnb, ident.1 ="C1", ident.2 ="A1")
#ident.1 ="one sample name",ident.2 ="Another sample name from the same patient"

## Save differentially expressed genes data for constructing PPI network in Cytoscape
write.table(monocyte.de.markers,"C1_A1_Tcell.xls",sep='\t',quote = F)
#write.table(monocyte.de.markers,"C1_A1_epithelial.xls",sep='\t',quote = F)
