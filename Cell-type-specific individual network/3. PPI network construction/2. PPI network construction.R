library(dynamicTreeCut)
library(openxlsx)
library(stringr)
library(Matrix)
library(WGCNA)
library(data.table)

Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

# Download protein network data and list of STRING proteins from https://cn.string-db.org/cgi/download (Restrict data to a specific organism:Homo sapiens)
idmapping <- fread(file="9606.protein.info.v12.0.txt.gz")
ppi <- fread(file="9606.protein.links.v12.0.txt.gz")

# Load the differential expression gene data of all samples
gene <- read.csv("genelist.csv")
# Select the differential expression genes between the epithelial cells/T cells of two samples from one patient
gene <- gene["C3_A4_Tcell.upno"]

# Match SYMBOL，ENSP ID
colnames(gene)[1]<-"preferred_name"
n<-merge(gene,idmapping[,c(1,2)],by="preferred_name",all.x=T)

prots_mat <- n

# Construct PPI networks of target genes
colnames(prots_mat)[2]<-"string_protein_id"
ppi3 <-  ppi[ppi$protein1 %in% prots_mat$string_protein_id, ]
ppi3 <-  ppi3[ppi3$protein2 %in% prots_mat$string_protein_id, ]
ppi <- ppi3

# Select combined score
ppi <- subset(ppi,combined_score>600)
 
# Output file for Cytoscape analysis: save as CSV file with protein column names of node1 and node2
colnames(ppi)[1] <- "node1"
colnames(ppi)[2] <- "node2"
write.csv(ppi,"C3_A4_Tcell.upno.csv")
#Gene list can also be input into the string database（https://cn.string-db.org）to construct PPI networks

# Convert PPI data type for graphIndex analysis
edges <- data.frame(from=ppi[,1],to=ppi[,2])
library(igraph)
g <- graph.data.frame(edges, directed = FALSE)  

network_igraph <- simplify(
  g,
  remove.multiple = TRUE,
  remove.loops = TRUE,
  edge.attr.comb = igraph_opt("edge.attr.comb")
)

network_graphNEL <- as_graphnel(network_igraph)

# GraphIndex analysis
library(QuACN)
graphIndex <- graphIndexComplexity(network_graphNEL)






