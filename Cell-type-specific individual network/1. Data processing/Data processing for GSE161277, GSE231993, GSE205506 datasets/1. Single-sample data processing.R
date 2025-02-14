if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
if(!require(ggplot2))install.packages("ggplot2")

library(multtest)
library(Seurat)
library(dplyr)
library(patchwork)
library(R.utils)
library(ggplot2)

#Download GSE161277 dataset from GEO 
#Each sample's three data files are placed in one folder and are respectively named as "barcode.tsv.gz","features.tsv.gz","matrix.mtx.gz"
#Load one sample data
pbmc.data <-Read10X(data.dir ="data location/patient3_carcinoma")

pbmc <- CreateSeuratObject(counts= pbmc.data, project="pbmc3k", min.cells = 3, min.features=200)

## Initialize the seurat object with the raw (non-normalized data). 
lalala <- as.data.frame(pbmc[["RNA"]]@counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1="nCount_RNA", feature2 ="percent.mt")
plot2 <- FeatureScatter(pbmc, feature1="nCount_RNA", feature2 ="nFeature_RNA")
CombinePlots(plots =list(plot1, plot2))

pbmc <- subset(pbmc, subset= nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <15)

ncol(as.data.frame(pbmc[["RNA"]]@counts))
pbmc <- NormalizeData(pbmc, normalization.method="LogNormalize", scale.factor=10000)

pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nFeatures = 2000)
pbmc = FindVariableFeatures(pbmc)

top10 <- head(VariableFeatures(pbmc), 10) 
plot1 <-VariableFeaturePlot(pbmc)
plot2 <-LabelPoints(plot=plot1, points=top10, repel= TRUE)

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

# Save data
saveRDS(pbmc,"pbmc_doublet.rds")



