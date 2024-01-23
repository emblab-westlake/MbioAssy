###@email zhaoze@westlake.edu.cn

# Scripts for iCAMP null model analysis for community assembly
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

# the treatment informaiton table (group information)
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
