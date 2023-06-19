###@email yuanling@westlake.edu.cn
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
