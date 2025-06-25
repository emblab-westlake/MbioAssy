###@email zhaoze@westlake.edu.cn
###update: add the function of using RMT theory to find the optimized correlation threshold

# input includes abundance table of microbial entities (e.g., OTUs, ASVs),
# each row is a sample, each column is an OTU
table = read.table(
  'example_input/OTU.txt',sep = '\t',header = T,row.names = 1)
table <- as.matrix(table)
table <- table[which(rowSums(table) > 0),]
table <- table[,which(colSums(table) > 0)]


# 5.2
# RMT-based Co-occurrence network construction
# Reference:Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. 
# Taxonomic relatedness shapes bacterial assembly in activated sludge of 
# globally distributed wastewater treatment plants. Environmental Microbiology. 16(8):2421-2432

if (!requireNamespace("vegan", quietly=TRUE))
  install.packages("vegan")
if (!requireNamespace("igraph", quietly=TRUE))
  install.packages("igraph")
if (!requireNamespace("Hmisc", quietly=TRUE))
  install.packages("Hmisc")
if (!requireNamespace("RMThreshold", quietly=TRUE))
  install.packages("RMThreshold")
library(vegan)
library(igraph)
library(Hmisc)
library(RMThreshold)

# define function co_occurrence_network
# to construct co-occurrence network
co_occurrence_network <- 
  function(matrix, p.cutoff){
  
  # correlation analysis based on spearman's co-efficient
  matrix.dist<-rcorr(t(matrix),type="spearman")
  matrix.cor<-matrix.dist$r
  matrix.cor.p<-matrix.dist$P
  
  # using RMT theory to find the optimized correlation threshold or not
  cor.cutoff <- RMT_cutoff(matrix.cor)
  
  #Multiple testing correction using Benjamini-Hochberg standard false discovery rate correction ("FDR-BH")
  matrix.cor.p <- p.adjust(matrix.cor.p, method="BH")
  
  # Consider positive cooccurence at given coefficient (cor.cutoff) and p-value cutoffs
  matrix.cor1<-matrix.cor
  matrix.cor1.p<-matrix.cor.p
  matrix.cor1[which(matrix.cor1 <= cor.cutoff)]=0
  matrix.cor1[which(matrix.cor1.p>p.cutoff)]=0
  # delete those rows and columns with sum = 0
  matrix.cor1<-matrix.cor1[which(rowSums(matrix.cor1)!=1),]
  matrix.cor1<-matrix.cor1[,which(colSums(matrix.cor1)!=0)]
  
  # generate graph using igraph
  g1<-graph_from_adjacency_matrix(matrix.cor1,weight=T,mode="undirected")
  g1<-simplify(g1)
  V(g1)$label <- V(g1)$name
  V(g1)$degree <- degree(g1)
  
  # append the output into results
  result<-list()
  result$matrix.cor<-matrix.cor
  result$matrix.cor.p<-matrix.cor.p
  
  result$matrix.cor1<-matrix.cor1
  result$graph1<-g1
  result$cor.cutoff <- cor.cutoff
  
  return(result)
}

# define function RMT_cutoff
# using RMT theory to find the optimized correlation threshold
RMT_cutoff <- function(matrix.cor){
    
    res <-
      rm.get.threshold(matrix.cor, interactive = F,dist.method = 'LL',
                       plot.comp = F, save.fit = F, plot.spacing = F,
                       interval = c(0.6,0.99))
    ## click "ESC" on the keyboard twice, this will not affect the following analysis
    if (max(res$p.ks > 0.05)) {
      pks.id <- which(res$p.ks >= 0.05)[1]
    } else {
      pks.id <- which.max(res$p.ks)
    }
    
    RMT.cutoff <- res$tested.thresholds[pks.id]
    RMT.cutoff.round <- round(RMT.cutoff, 2)
    if (isTRUE(RMT.cutoff.round<RMT.cutoff)) {
      RMT.cutoff.round <- RMT.cutoff.round + 0.01
    }
    return(RMT.cutoff.round)
  }



# Construct co-occurrence network using defined function co_occurrence_network and output results
# Creating gml files of network (to be visulized in Gephi or Cytoscape)


## Notice: when the figure appears in Plots window, please click "ESC" on the keyboard twice. This will not affect the following analysis ##
pattern <- co_occurrence_network(
  matrix = table, p.cutoff = 0.05)  # cutoffs for correlation coefficient and P-value
pattern$cor.cutoff
View(pattern$matrix.cor1)

write_graph(pattern$graph1,'Network.gml',format='gml')    #network file for positive association

# Calculating network topological properties
g<-pattern$graph1   ###positive network
c <- cluster_walktrap(g)
# Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g, vids = NULL,
                   weights = NULL)
spl <- mean_distance(g, directed=FALSE, unconnected=TRUE)
gd  <- edge_density(g, loops=FALSE)
nd  <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)
node.degree <- degree(g, v = V(g), mode="all")
ad  <- mean(node.degree)
e <- ecount(g)
v <- vcount(g)
global.topology <- data.frame(e,v,cc,spl,md,gd,nd,ad)
write.csv(global.topology, file="Network.global.topology.csv")

# Node toplogical features
betweenness.centrality <- betweenness(g, v=V(g), 
                                      directed = FALSE, weights = NA,
                                      normalized = FALSE)
closeness.centrality <- closeness(g, vids = V(g),
                                  weights = NA, normalized = FALSE)
node.transitivity <- transitivity(g, type = c("local"), vids = NULL,
                                  weights = NA)
node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file="Network.node.topology.csv")
