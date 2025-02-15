install.packages("Seurat")
install.packages("tidyverse")
install.packages("devtools")
install.packages("WGCNA")

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)

if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
   BiocManager::install(c("GO.db", "impute", "preprocessCore"))
  install.packages("tester")

devtools:: install_github ('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)
#If the installation of the hdWGCNA package fails, you can try switching to a different version of R.

# load the clustered data
pbmc <- readRDS("cancer_batch_cluster.rds") 

## Select the T cell or epithelial cell data from one sample
#T cells
seurat_obj <- subset(x = pbmc, subset =orig.ident =="C1",idents = 0 )
#epithlial cells
#seurat_obj <- subset(x = pbmc, subset =orig.ident =="C1",idents = c(2,4,6,10))
#orig.ident == "fill in the sample name"
#ident == "fill in the cell cluster number annotated as T cells or epithelial cells"

# Setup data for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# Construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("seurat_clusters", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # 'umap' select the dimensionality reduction to perform KNN on
  k = 30, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'orig.ident' # set the Idents of the metacell seurat object
)

# Normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

## Set up the expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "C1", # the name of the group of interest in the group.by column
  group.by='orig.ident', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)


## Select soft-power threshold

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)


# Construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=NULL,  
  setDatExpr=FALSE,
  tom_name = '3seddvs' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='INH hdWGCNA Dendrogram')
TOM <- GetTOM(seurat_obj)

# Need to run ScaleData first or else harmony throws an error
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# Compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=NULL 
)

# Harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# Module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

## Compute module connectivity
# Compute eigengene-based connectivity (kME):
install.packages("qlcMatrix")
library(qlcMatrix)
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = NULL, group_name = 'INH'
)

## Individual module network plots
ModuleNetworkPlot(seurat_obj)


## Combined hub gene network plots
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 6, n_other=3,
  edge_prop = 0.5,
  mods = 'all'
)

## Applying UMAP to co-expression networks
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# Get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, 
  label_hubs=2 ,
  keep_grey_edges=FALSE
  
)

network_igraph <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
saveRDS(network_igraph,"C1_Tcell.RDS")

## Convert data types for network topology analysis
# Convert the igraph object into a graphNEL object for QUACN R package analysis
install.packages("igraph")
library(igraph)
network_igraph <- simplify(
  network_igraph,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = igraph_opt("edge.attr.comb")
)

network_graphnel <- as_graphnel(network_igraph)

if (!require ("QuACN", quietly = TRUE)) install.packages('QuACN') ##本地安装的
if (!require ("BiocManager", quietly = TRUE)) install.packages ("BiocManager") BiocManager::install ("RBGL")
if (!require ("combinat", quietly = TRUE)) install.packages('combinat')
BiocManager::install("graph")
library(RBGL)
library(QuACN) 
library(graph)
library(combinat) 

graphIndex <- graphIndexComplexity(network_graphNEL)

# Convert the igraph object into a data frame and imported into Cytoscape for network analysis
network_frame <- as_data_frame(network_igraph)
write.csv(network_frame,"N1_Tcell_allgenes.csv")
