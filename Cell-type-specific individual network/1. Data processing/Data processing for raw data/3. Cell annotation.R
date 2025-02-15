library(SingleR)
library(Seurat)
library(celldex)

#Load the clustered data
cancer_cluster <- readRDS("cancer_batch_cluster.rds") 

# Cell type annotation
ref_use <- HumanPrimaryCellAtlasData() 
pbmc <- cancer_cluster

pred <- SingleR(test=as.matrix(pbmc@assays$RNA@data),      
                ref=ref_use,                                
                labels=ref_use$label.fine,                  
                clusters = pbmc$seurat_clusters,            
                method = "cluster")                       

# Save the SingleR annotation results
write.table(pred,"SinglerR_annotation.xls",sep='\t',quote = F)
