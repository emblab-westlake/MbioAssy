## MbioAssy1.0
### This folder contains R scripts for microbial community assembly analyses, including NST calculation, neutral model analysis, C-score variance analysis and co-occurrence network analysis.
input includes a) abundance table of microbial entities (e.g., OTUs, ASVs), each row is a sample, each column is an OTU.
               b) a one-column matrix indicating the group of each sample (only used in NST calculation)

Note: a minimum number of 20 to 25 samples is suggested for the analyses;
      Time-series data or local scale samples are recommended for the analyses, samples of large spatial scales may lead to strong bias.


### Functions:
### 1_NST.R
This script uses a developed R package NST (Reference doi:10.1073/pnas.1904623116) to calculate NST (normalized stochasticity ratio), which range from 0 to 100% and provides quantitative assessment of ecological stochasticity based on null model. A larger NST represents the community is more affected by stochasticity.
Reference:https://doi.org/10.1101/2020.06.03.131896


### 2_Neutral_model.R
This script was modified from (Reference doi:10.1038/ismej.2015.142) and used to compare observed community composition with that predicted by a neutral assembly model(Reference doi:10.1111/j.1462-2920.2005.00956.x),which assumes that community assembly is only driven by chance and dispersal.The results were output as a visualized figure and a detailed table which contains information of observed and predicted community composition.
Reference:https://doi.org/10.1101/2020.06.03.131896


### 3_C-score-var.R
This script was for the checkerboard score variance(C-score-var) analysis of non-random species co-occurrence patterns. Cite:

1. Niederdorfer R, Fragner L, Yuan L, Hausherr D, Wei J, Magyar P, Joss A, Lehmann MF, Ju F, Helmut Bürgmann. 2021. Distinct growth stages controlled by the interplay of deterministic and stochastic processes in functional anammox biofilms. Water Research. 200:117225
2. Hu AY, Ju F, Hou LY, Li JW, Yang XY, Wang HJ, Mulla SI, Sun Q, Bürgmann H, Yu CP. 2017. Strong impact of anthropogenic contamination on the co-occurrence patterns of a riverine microbial community. Environmental Microbiology (2017) 19(12), 4993–500

### 4_Co-occurrence_network.R
This script first defines the function "co_occurrence_network" to generate a filtered gml-formatted network file based on all pairwise Spearman's correlations and FDR-adjusted P-value calculation and filtration, and then uses input abundance table to generate co-occurrence network which can be visulized in Gephi (https://gephi.org/) and calculate network topological properties.

Reference:
1. Ju F, Zhang T. 2015. Bacterial assembly and temporal dynamics in activated sludge of a full-scale municipal wastewater treatment plant. The ISME Journal. 9: 683-695

3. Ju F, Xia Y, Guo F, Wang ZP, Zhang T*.2014. Taxonomic relatedness shapes bacterial assembly in activated sludge ofglobally distributed wastewater treatment plants. Environmental Microbiology.16(8):2421-2432


### MbioAssy1.0.R
This script integrates the above four modules which use the same abundance table of microbial entities as input. Pls cite the link and corresponding papers if you find the scripts useful:
https://github.com/emblab-westlake/MbioAssy1.0/
