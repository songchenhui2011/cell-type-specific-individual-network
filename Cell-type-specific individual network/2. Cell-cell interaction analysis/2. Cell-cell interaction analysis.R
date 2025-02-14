# Load re-clustered T cells or epithelial cells data
Tcell <- readRDS("Tcell.rds")
epithelial <- readRDS("epithelial.rds")

Tcell$clusters <- paste0("T",Tcell$seurat_clusters)
epithelial$clusters <- paste0("E",epithelial$seurat_clusters)

# Merge data
sce2 <-merge(Tcell,epithelial)

# Load the required libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
devtools::install_github("sqjin/CellChat")

library(CellChat)
library(tidyverse)
library(Seurat)
options(stringsAsFactors = FALSE)

# Create a CellChat object
cellchat <- createCellChat(object = sce2,
                           meta = sce2@meta.data,
                           group.by = "clusters")

# Set the ligand-receptor interaction database  
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)

# Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
# use all CellChatDB for cell-cell communication analysis:
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

## Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)  # This step is necessary even if using the whole database
future::plan("multisession", workers = 1) # do parallel 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#project gene expression data onto PPI 
cellchat <- projectData(cellchat, PPI.human)

### Inference of cell-cell communication network
#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network 
cellchat <- aggregateNet(cellchat)

### Visualization of cell-cell communication network
#Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat,measure="count",color.heatmap = "Reds")
netVisual_heatmap(cellchat,measure="weight",color.heatmap = "Reds")



