###MbioAssy 2.0
###This script integrated four modules including NST calculation, neutral model analysis, iCAMP null model analysis, C-score variance analysis and co-occurrence network analysis,
###which aim to assess ecological stochasticity and determinism under microbial community assembly.
###email:yuanling@westlake.edu.cn zhaoze@westlake.edu.cn
###Reference: Zhao, Z., et al. (2023) Hydrodynamic and anthropogenic disturbances co-shape microbiota rhythmicity and community assembly within intertidal groundwater-surface water continuum. Water Research, 120236.

# input files includes:
# a) OTU table of microbial entities (e.g., OTUs, ASVs), each row is an OTU, each column is a sample.
# b) Group table indicating the grouping infprmation of samples.
# c) Taxonomy table of microbia entities (e.g., OTUs, ASVs), each row is an OTU, seven columns corresponding to the taxonomy classification of kingdom, phylum, class, order, family, genus and species.
# d) A phylogenetic tree constructed using the OTUs or ASVs.
# e) Environment table indiacting the phygeochemical variables of samples.


table <- t(read.table(
  'example_input/OTU.txt',sep = '\t',header = T, row.names = 1))
sample_group = read.table(
  'example_input/Group.txt',sep = '\t', header = T, row.names = 1)
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
p <- plot_sncm_fit(nm.out,title = 'Test')
pdf('Neutral.model.plot.pdf',width = 6,height = 4)
p
dev.off()
write.table(nm.out$predictions,file = 'Neutral.model.details.txt',sep = '\t')



#---------------------------------------------------------------------------------------
# 3 
# iCAMP null model analysis
# This Script was modified from github: https://github.com/DaliangNing/iCAMP1
# Reference: https://www.nature.com/articles/s41467-020-18560-z (Ning D, et al., 2020)
# If you are dealing with a large dataset (e.g., > 20000 taxa), a server with enough CPU threads (e.g., >20) is preferred to finish the calculation in reasonable time.

if (!requireNamespace("iCAMP", quietly = TRUE))
  install.packages("iCAMP")
if (!requireNamespace("ape", quietly = TRUE))
  install.packages("ape")
require(iCAMP)
require(ape)


# 1 # set folder paths and file names, please change according to the folder paths and file names in your computer.

# create a folder for the iCAMP analysis
dir.create('D:/MbioAssy/iCAMP_analysis')

# create a subfolder and put all input files in that folder
data.wd <- 'D:/MbioAssy/iCAMP_analysis/Data'
if (!dir.exists(data.wd)) {dir.create(data.wd)}

# the OTU table file (Tab delimited txt file)
com.file <- "OTU.txt"

# the phylogenetic tree file
tree.file <- "Tree.nwk"

# the treatment informaiton table
treat.file <- "Group.txt"

# the environmental varialbes
env.file <- "Env.txt"

# the classification (taxonomy) information
clas.file <- "Taxonomy.txt"

# create a subfolder to save the output files
save.wd <- "D:/MbioAssy/iCAMP_analysis/Output"
if (!dir.exists(save.wd)) {dir.create(save.wd)}

# 2 # key parameter setting
prefix <- "Test"  # prefix of the output file names. usually use a project ID.
rand.time <- 1000  # randomization time, 1000 is usually enough. For example test, you may set as 100 or less to save time.
nworker <- 50 # nworker is thread number for parallel computing, which depends on the CPU core number of your computer.
memory.G <- 500 # to set the memory size as you need (but should be less than the available space in your hard disk), so that calculation of large tree will not be limited by physical memory. unit is Gb.


# 3 # load data
setwd(data.wd)
comm <- t(read.table(
  com.file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  as.is = TRUE,
  stringsAsFactors = FALSE,
  comment.char = "",
  check.names = FALSE
))
tree <- read.tree(file = tree.file)
treat <- read.table(
  treat.file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  as.is = TRUE,
  stringsAsFactors = FALSE,
  comment.char = "",
  check.names = FALSE)

env <- read.table(
  env.file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  as.is = TRUE,
  stringsAsFactors = FALSE,
  comment.char = "",
  check.names = FALSE) # skip this if you do not have env.file

clas <- read.table(
  clas.file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  as.is = TRUE,
  stringsAsFactors = FALSE,
  comment.char = "",
  check.names = FALSE)


# 4 # match sample IDs in OTU table and treatment information table
sampid.check <-
  match.name(rn.list = list(comm = comm, treat = treat, env = env))
# sampid.check=match.name(rn.list=list(comm=comm,treat=treat)) # if you do not have env.file
# for your data files, if you have not matched their IDs, the unmatched samples will be removed.
treat <- sampid.check$treat
comm <- sampid.check$comm
comm <- comm[, colSums(comm) > 0, drop = FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
# env=sampid.check$env # skip this if you do not have env.file

# 5 # match OTU IDs in OTU table and tree file
spid.check <- match.name(
  cn.list = list(comm = comm),
  rn.list = list(clas = clas),
  tree.list = list(tree = tree))
comm <- spid.check$comm
clas <- spid.check$clas
tree <- spid.check$tree

# 6 # calculate pairwise phylogenetic distance matrix.
setwd(save.wd)
if (!file.exists("pd.desc")) {
  pd.big <- iCAMP::pdist.big(
    tree = tree,
    wd = save.wd,
    nworker = nworker,
    memory.G = memory.G)
} else {
  # if you already calculated the phylogenetic distance matrix in a previous run
  pd.big <- list()
  pd.big$tip.label <- read.csv(
    paste0(save.wd, "/pd.taxon.name.csv"),
    row.names = 1,
    stringsAsFactors = FALSE)[, 1]
  pd.big$pd.wd <- save.wd
  pd.big$pd.file <- "pd.desc"
  pd.big$pd.name.file <- "pd.taxon.name.csv"
}


# you may skip step 7-8, if the "alternative way" based on stochasticity is applicable, as mentioned in the method part of iCAMP paper (Ning et al 2020 Nature Communications).

# 7 # assess niche preference difference between species
# env is required for this step.
setwd(save.wd)
niche.dif <- iCAMP::dniche(
  env = env,
  comm = comm,
  method = "niche.value",
  nworker = nworker,
  out.dist = FALSE,
  bigmemo = TRUE,
  nd.wd = save.wd)

# 8 # within-bin phylogenetic signal assessment.
# For real data, you may try several different settings of binning, and choose the one leading to the best within-bin phylogenetic signal.
# env is required for this step.
# 8.1 # try phylogenetic binning using current setttings.
ds <- 0.2
bin.size.limit <- 24 #For real data, usually try 12, 24 and 48 to explore the best choice.
phylobin <- taxa.binphy.big(
  tree = tree,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  ds = ds,
  bin.size.limit = bin.size.limit,
  nworker = nworker)

# 8.2 # test within-bin phylogenetic signal.
sp.bin <- phylobin$sp.bin[, 3, drop = FALSE]
sp.ra <- colMeans(comm / rowSums(comm))
abcut <- 3
commc <- comm[, colSums(comm) >= abcut, drop = FALSE]
dim(commc)
spname.use <- colnames(commc)
binps <- iCAMP::ps.bin(
  sp.bin = sp.bin,
  sp.ra = sp.ra,
  spname.use = spname.use,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  nd.list = niche.dif$nd,
  nd.spname = niche.dif$names,
  ndbig.wd = niche.dif$nd.wd,
  cor.method = "pearson",
  r.cut = 0.1,
  p.cut = 0.05,
  min.spn = 5)

if (file.exists(paste0(prefix, ".PhyloSignalSummary.csv"))) {
  appendy <- TRUE
  col.namesy <- FALSE
} else {
  appendy <- FALSE
  col.namesy <- TRUE
}

write.table(
  data.frame(ds = ds, n.min = bin.size.limit, binps$Index),
  file = paste0(prefix, ".PhyloSignalSummary.csv"),
  append = appendy,
  quote = FALSE,
  sep = ",",
  row.names = FALSE,
  col.names = col.namesy)

if (file.exists(paste0(prefix, ".PhyloSignalDetail.csv"))) {
  appendy2 <- TRUE
  col.namesy2 <- FALSE
} else{
  appendy2 <- FALSE
  col.namesy2 <- TRUE
}

write.table(
  data.frame(
    ds = ds,
    n.min = bin.size.limit,
    binID = rownames(binps$detail),
    binps$detail),
  file = paste0(prefix, ".PhyloSignalDetail.csv"),
  append = appendy2,
  quote = FALSE,
  sep = ",",
  row.names = FALSE,
  col.names = col.namesy2)

# usually, you are looking for a binning setting lead to higher RAsig.abj (relative abundance of bins with significant phylogenetic signal) and relative high meanR (mean correlation coefficient across bins).
# see help document of the function "ps.bin" for the meaning of output.


# 9 # iCAMP analysis
# 9.1 # without omitting small bins.
bin.size.limit <- 24 # For real data, usually use a proper number according to phylogenetic signal test or try some settings then choose the reasonable stochasticity level. our experience is 12, or 24, or 48.
sig.index <- "Confidence"
icres <- iCAMP::icamp.big(
  comm = comm,
  pd.desc = pd.big$pd.file,
  pd.spname = pd.big$tip.label,
  pd.wd = pd.big$pd.wd,
  rand = rand.time,
  tree = tree,
  prefix = prefix,
  ds = 0.2,
  pd.cut = NA,
  sp.check = TRUE,
  phylo.rand.scale = "within.bin",
  taxa.rand.scale = "across.all",
  phylo.metric = "bMPD",
  sig.index = sig.index,
  bin.size.limit = bin.size.limit,
  nworker = nworker,
  memory.G = memory.G,
  rtree.save = FALSE,
  detail.save = TRUE,
  qp.save = FALSE,
  detail.null = FALSE,
  ignore.zero = TRUE,
  output.wd = save.wd,
  correct.special = TRUE,
  unit.sum = rowSums(comm),
  special.method = "depend",
  ses.cut = 1.96,
  rc.cut = 0.95,
  conf.cut = 0.975,
  omit.option = "no",
  meta.ab = NULL)


# 10 # iCAMP bin level statistics
setwd(save.wd)

icbin <- iCAMP::icamp.bins(
  icamp.detail = icres$detail,
  treat = treat,
  clas = clas,
  silent = FALSE,
  boot = TRUE,
  rand.time = rand.time,
  between.group = TRUE)

save(icbin, file = paste0(prefix, ".iCAMP.Summary.rda")) # just to archive the result. rda file is automatically compressed, and easy to load into R.

write.csv(
  icbin$Pt,
  file = paste0(prefix, ".ProcessImportance_EachGroup.csv"),
  row.names = FALSE)

write.csv(
  icbin$Ptk,
  file = paste0(prefix, ".ProcessImportance_EachBin_EachGroup.csv"),
  row.names = FALSE)

write.csv(
  icbin$Ptuv,
  file = paste0(prefix, ".ProcessImportance_EachTurnover.csv"),
  row.names = FALSE)

write.csv(
  icbin$BPtk,
  file = paste0(prefix, ".BinContributeToProcess_EachGroup.csv"),
  row.names = FALSE)

write.csv(
  data.frame(
    ID = rownames(icbin$Class.Bin),
    icbin$Class.Bin,
    stringsAsFactors = FALSE),
  file = paste0(prefix, ".Taxon_Bin.csv"),
  row.names = FALSE)

write.csv(
  icbin$Bin.TopClass,
  file = paste0(prefix, ".Bin_TopTaxon.csv"),
  row.names = FALSE)

# output files:
# Test.iCAMP.Summary.rda: the object "icbin" saved in R data format. see help document of the function icamp.bins for description of each element in the object.
# Test.ProcessImportance_EachGroup.csv: Relative importance of each process in governing the turnovers in a group of samples.
# Test.ProcessImportance_EachBin_EachGroup.csv: Relative importance of each process in governing the turnovers of each bin among a group of samples.
# Test.ProcessImportance_EachTurnover.csv: Relative importance of each process in governing the turnovers between each pair of communities (samples).
# Test.BinContributeToProcess_EachGroup.csv: Bin contribution to each process, measuring the contribution of each bin to the relative importance of each process in the assembly of a group of communities.
# Test.Taxon_Bin.csv: a matrix showing the bin ID and classification information for each taxon.
# Test.Bin_TopTaxon.csv: a matrix showing the bin relative abundance; the top taxon ID, percentage in bin, and classification; the most abundant name at each phylogeny level in the bin.



#---------------------------------------------------------------------------------------
# 4
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
# 5.1
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

#-------------------------------------------------------------------------
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
  matrix = t(table), p.cutoff = 0.05)  # cutoffs for correlation coefficient and P-value
pattern$cor.cutoff
View(pattern$matrix.cor1)

write.graph(pattern$graph1,'Network.gml',format='gml')    #network file for positive association

# Calculating network topological properties
g<-pattern$graph1   ###positive network
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