library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(clustree)

###Merge data from multiple samples
# Load sample data
normal1 <- readRDS('patient1_normal.rds')
normal2 <- readRDS('patient2_normal.rds')
normal3 <- readRDS('patient3_normal.rds')
adenoma1 <- readRDS('patient1_adenoma.rds')
adenoma2 <- readRDS('patient2_adenoma.rds')
adenoma3_1 <- readRDS('patient3_adenoma1.rds')
adenoma3_2 <- readRDS('patient3_adenoma2.rds')
carcinoma1 <- readRDS('patient1_carcinoma.rds')
carcinoma2 <- readRDS('patient2_carcinoma.rds')
carcinoma3 <- readRDS('patient3_carcinoma.rds')
paracancer2 <- readRDS('patient2_paracancer.rds')

#Modify orig.ident to the sample name
normal1@meta.data$orig.ident <- "N1"
normal2@meta.data$orig.ident <- "N2"
normal3@meta.data$orig.ident <- "N3"
adenoma1@meta.data$orig.ident <- "A1"
adenoma2@meta.data$orig.ident <- "A2"
adenoma3_1@meta.data$orig.ident <- "A3"
adenoma3_2@meta.data$orig.ident <- "A4"
carcinoma1@meta.data$orig.ident <- "C1"
carcinoma2@meta.data$orig.ident <- "C2"
carcinoma3@meta.data$orig.ident <- "C3"
paracancer2@meta.data$orig.ident <- "P2"

# Merge data
CancerMerge <- merge(normal1, 
              y=c(normal2,normal3,adenoma1,adenoma2,adenoma3_1,adenoma3_2,carcinoma1,carcinoma2,carcinoma3,paracancer2), 
              add.cell.ids = c("normal1","normal2","normal3","adenoma1","adenoma2","adenoma3_1","adenoma3_2","carcinoma1","carcinoma2","carcinoma3","paracancer2"),
              project="ALL") 

table (CancerMerge$orig.ident)


#Batch effects between samples were removed using the Harmony 
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
library(harmony)
test.seu <- CancerMerge
test.seu <-  test.seu%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()

test.seu <- RunPCA(test.seu, features = VariableFeatures(object = test.seu),npcs = 50, verbose = FALSE)

seuratObj <- RunHarmony(test.seu,"orig.ident")

names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,dims=1:15,
                     reduction = "harmony")
saveRDS(seuratObj,"cancer_merged_batch.rds")

## clustering analysis
cancer_cluster =seuratObj
cancer_cluster <- FindNeighbors(cancer_cluster,reduction = "harmony",
                             dims=1:15)

##Select an appropriate resolution based on the cell clustering results
for (i in c(0.1,0.2,0.3)){
  cancer_cluster = FindClusters(cancer_cluster,resolution = i,algorithm =1)
}

clustree(cancer_cluster@meta.data,prefix= "RNA_snn_res.")
markers_genes <- FindAllMarkers(cancer_cluster,
                                logfc.threshold = 0.25,
                                min.pct = 0.25,
                                essay = "RNA")
write.csv(markers_genes,"markers_genes.csv")
saveRDS(cancer_cluster,"cancer_batch_cluster.rds")

top10 <- markers_genes %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
top5 <- markers_genes %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC)

write.csv(top10,"markers_genes_top10.csv")
write.csv(top5,"markers_genes_top5.csv")

## Visualize the marker genes of cell clusters
DotPlot(cancer_cluster,features = unique(top5$gene),assay="RNA")+RotatedAxis(
)+ ggplot2::coord_flip()


