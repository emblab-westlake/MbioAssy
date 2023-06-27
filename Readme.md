## MbioAssy2.0

### Update
The MbioAssy2.0 adds two function modules based on MbioAssy1.0:
* '3_iCAMP.R': quantify relative importance of basic community assembly processes at phylogenetic group level.
* '5.2_RMT_Co-occurrence_network.R': add the function to find the optimized correlation threshold instead of the 'cut-off' value based on the RMT theory.
The document was updated by Ze Zhao on June 27 2023
Reference: Zhao, Z., et al. (2023) Hydrodynamic and anthropogenic disturbances co-shape microbiota rhythmicity and community assembly within intertidal groundwater-surface water continuum. Water Research, 120236.

### This folder contains R scripts for microbial community assembly analyses, including NST calculation, neutral model analysis, iCAMP null model analysis, C-score variance analysis, co-occurrence network analysis, and generation of random newtorks.

input files includes:
a) OTU table of microbial entities (e.g., OTUs, ASVs), each row is an OTU, each column is a sample.
b) Group table indicating the grouping infprmation of samples.
c) Taxonomy table of microbia entities (e.g., OTUs, ASVs), each row is an OTU, seven columns corresponding to the taxonomy classification of kingdom, phylum, class, order, family, genus and species.
d) A phylogenetic tree constructed using the OTUs or ASVs.
e) Environment table indiacting the phygeochemical variables of samples.

Note: Time-series data or local scale samples are recommended for the analyses, but samples with large spatial scale or heterogeneities may lead to bias.
      

### Functions:
### 1_NST.R
This script uses a developed R package NST (Reference doi:10.1073/pnas.1904623116) to calculate NST (normalized stochasticity ratio), which range from 0 to 100% and provides quantitative assessment of ecological stochasticity based on null model. A larger NST represents the community is more affected by stochasticity.

Reference to cite: 
1. Zhang L, Yin W, Wang C, Zhang AJ, Zhang H, Ju F. 2021. Untangling Microbiota Diversity and Assembly Patterns in the World's Largest Water Diversion Canal. Water Research. 117617.  https://www.sciencedirect.com/science/article/pii/S0043135421008125

### 2_Neutral_model.R
This script was modified from (Reference doi:10.1038/ismej.2015.142) and used to compare observed community composition with that predicted by a neutral assembly model(Reference doi:10.1111/j.1462-2920.2005.00956.x),which assumes that community assembly is only driven by chance and dispersal.The results were output as a visualized figure and a detailed table which contains information of observed and predicted community composition.

Reference to cite: 
1. Zhang L, Yin W, Wang C, Zhang AJ, Zhang H, Ju F. 2021. Untangling Microbiota Diversity and Assembly Patterns in the World's Largest Water Diversion Canal. Water Research. 117617.  https://www.sciencedirect.com/science/article/pii/S0043135421008125
2. Niederdorfer R, Fragner L, Yuan L, Hausherr D, Wei J, Magyar P, Joss A, Lehmann MF, Ju F, Helmut Bürgmann. 2021. Distinct growth stages controlled by the interplay of deterministic and stochastic processes in functional anammox biofilms. Water Research. 200:117225

### 3_iCAMP.R
This script was modified from (Reference doi: 10.1038/s41467-020-18560-z) and used to quantify relative importance of basic community assembly processes at phylogenetic group ('bin') levels. The results were output as a detailed table which contains the relative importance of the five processes: homogeneous selection, heterogeneous selection, homoginizing dispersal, dispersal limitation, and drift.

Reference to cite:
1. Zhao Z, Zhang L, Zhang GQ, Gao H, Chen XG, Li L, Ju F. 2023. Hydrodynamic and anthropogenic disturbances co-shape microbiota rhythmicity and community assembly within intertidal groundwater-surface water continuum. Water Research, 120236.
2. Ning, D., Yuan, M., Wu, L. et al. A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nat Commun 11, 4717 (2020). https://doi.org/10.1038/s41467-020-18560-z

### 4_C-score-var.R

This script was for the checkerboard score variance(C-score-var) analysis of non-random species co-occurrence patterns. Cite:

Reference to cite:
1. Niederdorfer R, Fragner L, Yuan L, Hausherr D, Wei J, Magyar P, Joss A, Lehmann MF, Ju F, Helmut Bürgmann. 2021. Distinct growth stages controlled by the interplay of deterministic and stochastic processes in functional anammox biofilms. Water Research. 200:117225   https://www.sciencedirect.com/science/article/pii/S0043135421004231

### 5.1_Co-occurrence_network.R
This script first defines the function "co_occurrence_network" to generate a filtered gml-formatted network file based on all pairwise Spearman's correlations and FDR-adjusted P-value calculation and filtration, and then uses input abundance table to generate co-occurrence network which can be further visulized and explored in Gephi (https://gephi.org/) or Cytoscape (https://cytoscape.org/).

Reference to cite:
1. Ju F, Zhang T. 2015. Bacterial assembly and temporal dynamics in activated sludge of a full-scale municipal wastewater treatment plant. The ISME Journal. 9: 683-695
2. Ju F, Xia Y, Guo F, Wang ZP, Zhang T*.2014. Taxonomic relatedness shapes bacterial assembly in activated sludge ofglobally distributed wastewater treatment plants. Environmental Microbiology.16(8):2421-2432

### 5.2_RMT_Co-occurrence_network.R
Based on the '5.1_Co-occurrence_network.R' script, this script adds the function to find the optimized correlation threshold instead of the 'cut-off' value based on the RMT theory.

Reference to cite:
1. Zhao Z, Zhang L, Zhang GQ, Gao H, Chen XG, Li L, Ju F. 2023. Hydrodynamic and anthropogenic disturbances co-shape microbiota rhythmicity and community assembly within intertidal groundwater-surface water continuum. Water Research, 120236.
2. Ju F, Zhang T. 2015. Bacterial assembly and temporal dynamics in activated sludge of a full-scale municipal wastewater treatment plant. The ISME Journal. 9: 683-695
3. Ju F, Xia Y, Guo F, Wang ZP, Zhang T*.2014. Taxonomic relatedness shapes bacterial assembly in activated sludge ofglobally distributed wastewater treatment plants. Environmental Microbiology.16(8):2421-2432


### 6_Random_network.R
This script was for generating 10000 or other numbers of random network to topologically compared against an observed network

Reference to cite:
1. Ju F, Xia Y, Guo F, Wang ZP, Zhang T*.2014. Taxonomic relatedness shapes bacterial assembly in activated sludge ofglobally distributed wastewater treatment plants. Environmental Microbiology.16(8):2421-2432 
https://onlinelibrary.wiley.com/doi/10.1111/1462-2920.12355


### MbioAssy2.0
This script integrates the above modules which use the same abundance table of microbial entities as input. Pls cite the link and corresponding papers if you find the scripts useful. Pls cite the package as: 
MbioAssy2.0: https://github.com/emblab-westlake/MbioAssy2.0/

Publications:
?????
