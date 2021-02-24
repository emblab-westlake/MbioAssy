###MbioAssy 1.0
###This script integrated four modules including NST calculation, neutral model analysis, C-score variance analysis and co-occurrence network analysis,
###which aim to assess ecological stochasticity and determinism under microbial community assembly.
###email:yuanling@westlake.edu.cn

# input includes a) abundance table of microbial entities (e.g., OTUs, ASVs),
#                   each row is a sample, each column is an OTU
#                b) a one-column matrix indicating the group of each sample
table = read.table('example_input/AMB.txt',sep = '\t',header = T)
rownames(table) = table[,1]
sample_group = read.table('example_input/AMB.sample.group.txt',sep = '\t')
rownames(sample_group) = table[,1]
table = table[,-1]
table <- as.matrix(table)
table <- table[which(rowSums(table) > 0),]
table <- table[,which(colSums(table) > 0)]

# 1
# normalized stochasticity ratio (NST) calculation
# Reference: Daliang Ning, Ye Deng, James M. Tiedje, Jizhong Zhou. (2019) 
# A general framework for quantitatively assessing ecological stochasticity. 
# Proceedings of the National Academy of Sciences 116:34, 16892-16898. 

# NST calculation
if (!requireNamespace("NST", quietly=TRUE))
  install.packages("NST")
library("NST")
nst = tNST(comm = table, group = sample_group, 
           dist.method = "jaccard", abundance.weighted = TRUE, 
           rand = 20,null.model = "PF")
# for argument 'rand', 1000 is recommended; here set rand=20 to save test time
nst.sum=nst$index.grp
#View(nst.sum)
# output results
write.table(nst.sum,'NST.output.txt',sep="\t")

# --------------------------------------------------------------------------------------
# 2
# Neutral model
# This part of script was modified from a published research as below:
# Reference:https://www.nature.com/articles/ismej2015142 {Burns et al.} (2016)

if (!requireNamespace("minpack.lm", quietly=TRUE))
  install.packages("minpack.lm")
if (!requireNamespace("Hmisc", quietly=TRUE))
  install.packages("Hmisc")
if (!requireNamespace("stats4", quietly=TRUE))
  install.packages("stats4")
require(minpack.lm)
require(Hmisc)
require(stats4)

# Define function fit_sncm
# which returns several fitting statistics as well as predicted occurrence frequencies 
# for each ASV from an ASV table based on their abundance in the metacommunity
fit_sncm <- function(spp, pool=NULL, taxon=NULL){
  
  options(warn=-1)
  
  # Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  # Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  # Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  # Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  # Removes rows with any zero (absent in either source pool or local communities)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  # Calculate the limit of detection
  d = 1/N
  
  # Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.001))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  # Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  # Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  # Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  # Results
  fitstats <- data.frame(
    m=as.numeric(coef(m.fit)),
    m.ci=as.numeric(coef(m.fit)-m.ci[1]),
    poisLL=as.numeric(pois.mle@details$value),
    Rsqr=as.numeric(Rsqr), # measuring fit, # comparing fit differing datasets to the same model
    Rsqr.pois=as.numeric(Rsqr.pois),
    RMSE=as.numeric(RMSE), # measuring fit # comparing fit differing datasets to the same model
    RMSE.pois=as.numeric(RMSE.pois),
    AIC.pois=as.numeric(aic.pois),  # comparing differing models to the dataset
    BIC.pois=as.numeric(bic.pois), # comparing differing models to the dataset
    N=as.numeric(N),
    Samples=as.numeric(nrow(spp)),
    Richness=as.numeric(length(p)),
    Detect=as.numeric(d))
  
  A <- cbind(p, freq, freq.pred, pred.ci[,2:3])
  A <- as.data.frame(A)
  colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr')
  if(is.null(taxon)){
    B <- A[order(A[,1]),]
  } else {
    B <- merge(A, taxon, by=0, all=TRUE)
    row.names(B) <- B[,1]
    B <- B[,-1]
    B <- B[order(B[,1]),]
  }
  B <- B[!is.na(B$freq),]
  # fit_class for graphing
  B$fit_class <-"As predicted"
  B[which(B$freq < B$pred.lwr),"fit_class"]<- "Below prediction"
  B[which(B$freq > B$pred.upr),"fit_class"]<- "Above prediction"
  B[which(is.na(B$freq)),"fit_class"]<- "NA"
  
  # combine fit stats and predicitons into list
  i <- list(fitstats, B)
  names(i) <- c("fitstats", "predictions")
  return(i)
}

# define function plot_sncm_fit
# to plot the output from fit_sncm by ggpolt2
plot_sncm_fit <- function(spp.out, fill = NULL, title = NULL){
  
  tax_levels <- colnames(spp.out$predictions)[7:length(colnames(spp.out$predictions))-1]
  
  if(is.null(fill)){
    fill <- "fit_class"
  }
  
  r2_val <- paste("r^2 ==", round(spp.out$fitstats$Rsqr,4))
  m_val <- paste("m ==", round(spp.out$fitstats$m,4))
  df <- data.frame(t(table(spp.out$predictions$fit_class)))
  df <- df[,c(2,3)]
  colnames(df) <- c("Prediction", "AVS Abundance")
  
  p <- ggplot(data=spp.out$predictions)
  
  if(fill == "fit_class"){
    p <- p + geom_point(aes(x = log(p), y = freq, fill=eval(parse(text=fill))), shape =21, color="black", size =2, alpha=0.75)
    p <- p + scale_fill_manual(
      name = "Prediction",
      values = c("Above prediction" = "seagreen", "As predicted" = "black", "Below prediction" = "tan1", "NA" = "white"),
      breaks = c("Above prediction", "As predicted", "Below prediction", "NA"),
      labels = c(paste0("Above prediction (",round((df[1,2]/spp.out$fitstats$Richness)*100, 1),"%)"),
                 paste0("As predicted (",round((df[2,2]/spp.out$fitstats$Richness)*100, 1),"%)"),
                 paste0("Below Prediction (",round((df[3,2]/spp.out$fitstats$Richness)*100, 1),"%)"),
                 paste0("NA (",df[4,2],")")))
    
  }else if (fill %in% tax_levels){
    p <- p + geom_point(aes(x = log(p), y = freq, fill=eval(parse(text=fill))), shape =21, color="black", size =2, alpha=0.75)
    p <- p + scale_fill_discrete(name = "Taxon")
    
  } else{
    print(paste0("fill variable: ", fill, " is not a valid taxonomic level or fit_class"))
  }
  
  p <- p + geom_line(aes(x = log(p), y = freq.pred), color = "dodgerblue4", lwd=1.5)
  p <- p + geom_line(aes(x = log(p), y = pred.lwr), color = "dodgerblue4", linetype="dashed", lwd=1.5)
  p <- p + geom_line(aes(x = log(p), y = pred.upr), color = "dodgerblue4", linetype="dashed", lwd=1.5)
  p <- p + xlab("log(Mean Relative Abundance)")
  p <- p + ylab("Frequency")
  p <- p + ggtitle(title)
  p <- p + annotate("text", x=-5, y=0.65, size=5, label = r2_val, parse=TRUE)
  p <- p + annotate("text", x=-5, y=0.5, size=5, label = m_val, parse=TRUE)
  p <- p + theme_bw()
  p <- p + theme(panel.grid=element_blank(),element_line(size=1,colour="black"))
  return(p)
}

# Neutral model analysis and visualization using the example ASV table
nm.out <- fit_sncm(table)
p <- plot_sncm_fit(nm.out,title = 'AMB')
pdf('Neutral.model.plot.pdf',width = 6,height = 4)
p
dev.off()
write.table(nm.out$predictions,file = 'Neutral.model.details.txt',sep = '\t')

#---------------------------------------------------------------------------------------
# 3
# Checkerboard-score-var (C-score-var) analysis
# Reference: Stone L, Roberts A. 1990. The checkerboard score and species distributions.85(1):74-79. doi: 10.1007/BF00317345.
# Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. Taxonomic relatedness shapes bacterial assembly in activated sludge of globally distributed wastewater treatment plants. 
# Environmental Microbiology. 16(8):2421-2432

if (!requireNamespace("EcoSimR", quietly=TRUE))
  install.packages("EcoSimR")
if (!requireNamespace("devEMF", quietly=TRUE))
  install.packages("devEMF")
library("EcoSimR")
library("devEMF")
set.seed(56)        # for reproducible results

# Create presence-absence matrix
table01 <- t(table)
table01[table01 > 0] <- 1
# Filter out empty rows
table01.nonzerorow <- table01[which(rowSums(table01) > 0),]
table01 <- table01.nonzerorow

# C-score-var calculation  
csvarModel <- cooc_null_model(table01, algo = "sim9", metric = "c_score_var",
                              nReps = 500, saveSeed = FALSE, burn_in = 500, algoOpts = list(),
                              metricOpts = list(), suppressProg = FALSE)
# for argument 'nReps', 30000 is recommended; here set nReps = 500 to save test time

# output results
write.table('C-score-var summary',"c-score-var.summary.txt",append = TRUE)
sink("c-score-var.summary.txt", append = TRUE)
summary(csvarModel)
sink(NULL)

emf(file = "c-score-var.hist.500.emf", width = 7, height = 7,
    bg = "transparent", fg = "black", pointsize = 12,
    family = "Helvetica", custom.lty = FALSE);
plot(csvarModel,type = "hist");
dev.off()

#-------------------------------------------------------------------------
# 4
# Co-occurrence network construction
# Reference:Ju F, Xia Y, Guo F, Wang ZP, Zhang T. 2014. 
# Taxonomic relatedness shapes bacterial assembly in activated sludge of 
# globally distributed wastewater treatment plants. Environmental Microbiology. 16(8):2421-2432

if (!requireNamespace("vegan", quietly=TRUE))
  install.packages("vegan")
if (!requireNamespace("igraph", quietly=TRUE))
  install.packages("igraph")
if (!requireNamespace("Hmisc", quietly=TRUE))
  install.packages("Hmisc")
library(vegan)
library(igraph)
library(Hmisc)

# define function co_occurrence_network
# to construct co-occurrence network
co_occurrence_network<-function(matrix,cor.cutoff,p.cutoff){
  
  # correlation analysis based on spearman's co-efficient
  matrix.dist<-rcorr(t(matrix),type="spearman")
  matrix.cor<-matrix.dist$r
  matrix.cor.p<-matrix.dist$P
  
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
  g1<-graph.adjacency(matrix.cor1,weight=T,mode="undirected")
  g1<-simplify(g1)
  V(g1)$label <- V(g1)$name
  V(g1)$degree <- degree(g1)
  
  # append the output into results
  result<-list()
  result$matrix.cor<-matrix.cor
  result$matrix.cor.p<-matrix.cor.p
  
  result$matrix.cor1<-matrix.cor1
  result$graph1<-g1
  
  return(result)
}

# Construct co-occurrence network using defined function co_occurrence_network and output results
# Creating gml files of network (to be visulized in Gephi or Cytoscape)
pattern <- co_occurrence_network(t(table),0.8,0.05)  # cutoffs for correlation coefficient and P-value
write.graph(pattern$graph1,'AMB.Network.gml',format='gml')    #network file for positive association

# Calculating network topological properties
g <- pattern$graph1   ###positive network
c <- cluster_walktrap(g)
# Global toplogical features
modularity(c)
md <- modularity(g, membership(c), weights = NULL)
cc <- transitivity(g, vids = NULL,
                   weights = NULL)
spl <- average.path.length(g, directed=FALSE, unconnected=TRUE)
gd  <- graph.density(g, loops=FALSE)
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
                                      nobigint = TRUE, normalized = FALSE)
closeness.centrality <- closeness(g, vids = V(g),
                                  weights = NA, normalized = FALSE)
node.transitivity <- transitivity(g, type = c("local"), vids = NULL,
                                  weights = NA)
node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality, node.transitivity)
write.csv(node.topology, file="Network.node.topology.csv")