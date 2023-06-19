###@email yuanling@westlake.edu.cn
# input includes abundance table of microbial entities (e.g., OTUs, ASVs),
# each row is a sample, each column is an OTU
table = read.table('example_input/AMB.txt',sep = '\t',header = T)
rownames(table) = table[,1]
table = table[,-1]
table <- as.matrix(table)
table <- table[which(rowSums(table) > 0),]
table <- table[,which(colSums(table) > 0)]

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
