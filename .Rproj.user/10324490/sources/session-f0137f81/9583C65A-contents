#######################################################
#####         Calculating Genetic distances       #####
#######################################################

### Libraries
library(here); library(ape); 
library(factoextra); library(writexl)

### Read in alignment data:
alignment <- ape::read.dna(format = 'fasta', file = here::here('1_LRV_analysis/4_LRV_MSA/LRV_WGS_w_REFS.fa'))

align.dist <- dist.dna(x = alignment, model = 'raw', as.matrix = T)
fviz_dist(as.dist(align.dist), lab_size = 8)

## Export distance matrix to excell:
writexl::write_xlsx(x = as.data.frame(align.dist), path = here::here('1_LRV_analysis/5_ML_Phylogenies/2_Genetic_Distance/Genetic_distance_clades.xlsx'))
