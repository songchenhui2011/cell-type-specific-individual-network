library(SingleR)
library(Seurat)
library(celldex)

# Download scRNA-seg dataset from GEO
# Load seurat object which was downloaded from GEO
seuratObj <- readRDS("data location/seuratobj file.rds") 
# Load cell cluster data which was downloaded from GEO
group <- read.csv("data location/cluster file.csv")

# Select a sample data
seuratObj <- subset(seuratObj,subset =orig.ident =="sample name")

pbmc <- seuratObj

pbmc@meta.data$sample_cell_barcode <- rownames(pbmc@meta.data)

pbmc <- subset(pbmc, subset= nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <15)

ncol(as.data.frame(pbmc[["RNA"]]@counts))
pbmc <- NormalizeData(pbmc, normalization.method="LogNormalize", scale.factor=10000)

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nFeatures = 2000)
pbmc = FindVariableFeatures(pbmc)

pbmc <-ScaleData(pbmc, features=rownames(pbmc))

pbmc <- RunPCA(pbmc, features=VariableFeatures(object=pbmc))
print(pbmc[["pca"]], dims=1:5, nfeatures=5)
VizDimLoadings(pbmc, dims=1:2, reduction="pca")
DimPlot(pbmc, reduction="pca")
DimHeatmap(pbmc, dims=1, cells=500, balanced=TRUE)

pbmc <-JackStraw(pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(pbmc, dims=1:20)

pbmc <-FindNeighbors(pbmc, dims=1:10)
pbmc <-FindClusters(pbmc, resolution=0.5)

pbmc <- RunUMAP(pbmc, dims=1:10)

pbmc@meta.data <- merge(pbmc@meta.data,group,by="sample_cell_barcode",all.x=TRUE)

rownames(pbmc@meta.data) <- pbmc@meta.data$sample_cell_barcode

# Extract the cell count for each cell type in each sample.
library(reshape)
library(tidyverse)
plotC <- table(pbmc@meta.data$renamed_cellactivity_clusters) %>% melt()
colnames(plotC) <- c("CellType","Number")
write.csv(plotC,"cell_count.csv")

# Save the seurat object for each sample
seuratObj <- readRDS("sample name.RDS") 
