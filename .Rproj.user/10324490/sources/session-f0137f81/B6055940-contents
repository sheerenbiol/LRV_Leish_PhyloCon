################################################################
##### Comparing Leish - LRV genetic distances through MMRR #####
################################################################

### Libraries:
##############
library(here)
library(readxl); library(ape)
library(vegan); library(ggpubr)
library(ggplot2); library(data.table)
library(geodist); library(stringr); library(clusterSim)


### 1. Old World + New World Leish vs. LRV1 + LRV2
##################################################
## Binary matrix with all Leish-LRV links
link.data <- read.table(here::here("3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/Link_data.csv"), header = T, sep=';')
rownames(link.data) <- link.data$X; link.data <- link.data[,-c(1)]
colnames(link.data) <- gsub("\\.", "-", colnames(link.data))
colnames(link.data)[[3]] <- '013-2023-Lb-2018'
link.data <- as.matrix(link.data)

### LRV:
LRV.aln <- read.dna(here::here('3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/1_LRV/OW-NW_LRV_Pos_VG_Assembllies_REFS_Concat_genes.fa'), format = 'fasta')
LRV.aln.dist <- dist.dna(LRV.aln, model = 'K80', as.matrix = T)
LRV.aln.dist <- LRV.aln.dist[which(rownames(LRV.aln.dist) %in% colnames(link.data)),which(colnames(LRV.aln.dist) %in% colnames(link.data))]

## Fix order of row/colnames:
new_order <- sort(rownames(LRV.aln.dist))
LRV.aln.dist <- LRV.aln.dist[new_order,new_order]

## Make rownames/colnames a bit more similar to leish:
rownames(LRV.aln.dist)[25] <- colnames(LRV.aln.dist)[25] <- 'Lma-29ch'; rownames(LRV.aln.dist)[26] <- colnames(LRV.aln.dist)[26] <- 'Lma-37ch'; rownames(LRV.aln.dist)[27] <- colnames(LRV.aln.dist)[27] <- 'Lma-3T'; rownames(LRV.aln.dist)[28] <- colnames(LRV.aln.dist)[28] <- 'Lma-44T9'; rownames(LRV.aln.dist)[29] <- colnames(LRV.aln.dist)[29] <- 'Lb-CUM24'; rownames(LRV.aln.dist)[30] <- colnames(LRV.aln.dist)[30] <- 'Lb-CUM41'; rownames(LRV.aln.dist)[31] <- colnames(LRV.aln.dist)[31] <- 'Lb-LC2318'; rownames(LRV.aln.dist)[32] <- colnames(LRV.aln.dist)[32] <- 'Lb-PER010'; rownames(LRV.aln.dist)[33] <- colnames(LRV.aln.dist)[33] <- 'Lb-PER012'; rownames(LRV.aln.dist)[34] <- colnames(LRV.aln.dist)[34] <- 'Lb-PER065'; rownames(LRV.aln.dist)[35] <- colnames(LRV.aln.dist)[35] <- 'Lb-PER186'; rownames(LRV.aln.dist)[36] <- colnames(LRV.aln.dist)[36] <- 'Lb-PER130'; rownames(LRV.aln.dist)[37] <- colnames(LRV.aln.dist)[37] <- 'Lb-PER014'; rownames(LRV.aln.dist)[38] <- colnames(LRV.aln.dist)[38] <- 'Lb-PER067'; rownames(LRV.aln.dist)[39] <- colnames(LRV.aln.dist)[39] <- 'Lb-PER016'; rownames(LRV.aln.dist)[40] <- colnames(LRV.aln.dist)[40] <- 'Lb-LH825'; rownames(LRV.aln.dist)[41] <- colnames(LRV.aln.dist)[41] <- 'Lb-LC2319'; rownames(LRV.aln.dist)[42] <- colnames(LRV.aln.dist)[42] <- 'Lb-CUM68'; rownames(LRV.aln.dist)[43] <- colnames(LRV.aln.dist)[43] <- 'Lb-PER069'; rownames(LRV.aln.dist)[44] <- colnames(LRV.aln.dist)[44] <- 'Lb-LC2176'; rownames(LRV.aln.dist)[45] <- colnames(LRV.aln.dist)[45] <- 'Lb-LC2284'; rownames(LRV.aln.dist)[46] <- colnames(LRV.aln.dist)[46] <- 'Lb-PER002'; rownames(LRV.aln.dist)[47] <- colnames(LRV.aln.dist)[47] <- 'Lb-PER207'; rownames(LRV.aln.dist)[48] <- colnames(LRV.aln.dist)[48] <- 'Lb-PER212'; rownames(LRV.aln.dist)[49] <- colnames(LRV.aln.dist)[49] <- 'Lb-PER201'; rownames(LRV.aln.dist)[50] <- colnames(LRV.aln.dist)[50] <- 'Lad-LV30'; rownames(LRV.aln.dist)[51] <- colnames(LRV.aln.dist)[51] <- 'Lta-LV108'; 

## RE-DO Fix order of row/colnames:
new_order <- sort(rownames(LRV.aln.dist))
LRV.aln.dist <- LRV.aln.dist[new_order,new_order]

### Leish:
Leish.aln <- read.dna(here::here('3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/2_LEISH/4_Concatenated_alignment_MLtree/OW-NW_Concat_genes_Nomixes.fa'), format = 'fasta')
Leish.aln.dist <- dist.dna(Leish.aln, model = 'K80', as.matrix = T)
Leish.aln.dist <- Leish.aln.dist[which(rownames(Leish.aln.dist) %in% rownames(link.data)),which(colnames(Leish.aln.dist) %in% rownames(link.data))]

## Fix order of row/colnames:
new_order <- sort(rownames(Leish.aln.dist))
Leish.aln.dist <- Leish.aln.dist[new_order,new_order]

##### MMRR on standardized distance matrices:
#MMRR function: Wang, Ian J. (2013), Data from: Examining the full effects of landscape heterogeneity on spatial genetic variation: a multiple matrix regression approach for quantifying geographic and ecological isolation, Dryad, Dataset, https://doi.org/10.5061/dryad.kt71r
source(here::here('3_Cophylo_analysis/3_Genetic_Dist_Regression/R_Scripts/MMRR_Wang2013.R'))

##Creating a list of independent matrices:
X <- list()
X$Leish <- Leish.aln.dist

## MMRR with genetic distance as response variable
lrv.leish.mmrr <- MMRR(LRV.aln.dist, X, nperm = 1000) 
lrv.leish.mmrr

### Scatter plot for pairwise genetic distance:
###############################################
## Regression plots:
distances <- data.frame(lrv=LRV.aln.dist[lower.tri(LRV.aln.dist)], 
                        leish=Leish.aln.dist[lower.tri(Leish.aln.dist)])

ggplot(distances, aes(x=leish, y=lrv)) + 
  #geom_point(size=3, aes(color="#AAAC95", ))+ 
  scale_color_manual(values=alpha("#DFDF8D", 0.7))+
  geom_smooth(method = 'lm', formula = y~x, color = "#E57A61") + xlab('Leishmania genetic distance') + ylab('LRV genetic distance') +
  stat_summary_bin(fun.data='mean_cl_boot', bins = 200,colour=alpha("#DFDF8D", 0.8), size=1)+
  theme(
    panel.background = element_rect(fill = "#05598C" ,
                                    colour = "#05598C" ,
                                    size = 0.6, linetype = "solid"),
    panel.grid.major = element_line(size = 0.6, linetype = 'solid',
                                    colour =  "#929C96"), 
    panel.grid.minor = element_line(size = 0.35, linetype = 'solid',
                                    colour =  "#929C96")
  )



### 2. New World Leish vs. LRV1 (SNP-Based alignments)
######################################################

## Binary matrix with all Leish-LRV links
link.data <- read.table(here::here("3_Cophylo_analysis/1_Preparing_MSA/2_Species/Link_data.csv"), header = T, sep=';')
rownames(link.data) <- link.data$X; link.data <- link.data[,-c(1)]
colnames(link.data) <- gsub("\\.", "-", colnames(link.data))
colnames(link.data)[[3]] <- '013-2023-Lb-2018'
link.data <- as.matrix(link.data)

### LRV:
LRV.aln <- read.dna(here::here('1_LRV_analysis/4_LRV_MSA/LRV_WGS_w_REFS.fa'), format = 'fasta')
LRV.aln.dist <- dist.dna(LRV.aln, model = 'K80', as.matrix = T)
LRV.aln.dist <- LRV.aln.dist[which(rownames(LRV.aln.dist) %in% colnames(link.data)),which(colnames(LRV.aln.dist) %in% colnames(link.data))]

## Fix order of row/colnames:
new_order <- sort(rownames(LRV.aln.dist))
LRV.aln.dist <- LRV.aln.dist[new_order,new_order]

## Make rownames/colnames a bit more similar to leish:
rownames(LRV.aln.dist)[25] <- colnames(LRV.aln.dist)[25] <- 'CUM24'; rownames(LRV.aln.dist)[26] <- colnames(LRV.aln.dist)[26] <- 'CUM41'; rownames(LRV.aln.dist)[27] <- colnames(LRV.aln.dist)[27] <- 'LC2318'; rownames(LRV.aln.dist)[28] <- colnames(LRV.aln.dist)[28] <- 'PER010'; rownames(LRV.aln.dist)[29] <- colnames(LRV.aln.dist)[29] <- 'PER012'; rownames(LRV.aln.dist)[30] <- colnames(LRV.aln.dist)[30] <- 'PER065'; rownames(LRV.aln.dist)[31] <- colnames(LRV.aln.dist)[31] <- 'PER186'; rownames(LRV.aln.dist)[32] <- colnames(LRV.aln.dist)[32] <- 'PER130'; rownames(LRV.aln.dist)[33] <- colnames(LRV.aln.dist)[33] <- 'PER014'; rownames(LRV.aln.dist)[34] <- colnames(LRV.aln.dist)[34] <- 'PER067'; rownames(LRV.aln.dist)[35] <- colnames(LRV.aln.dist)[35] <- 'PER016'; rownames(LRV.aln.dist)[36] <- colnames(LRV.aln.dist)[36] <- 'LH825'; rownames(LRV.aln.dist)[37] <- colnames(LRV.aln.dist)[37] <- 'LC2319'; rownames(LRV.aln.dist)[38] <- colnames(LRV.aln.dist)[38] <- 'CUM68'; rownames(LRV.aln.dist)[39] <- colnames(LRV.aln.dist)[39] <- 'PER069'; rownames(LRV.aln.dist)[40] <- colnames(LRV.aln.dist)[40] <- 'LC2176'; rownames(LRV.aln.dist)[41] <- colnames(LRV.aln.dist)[41] <- 'LC2284'; rownames(LRV.aln.dist)[42] <- colnames(LRV.aln.dist)[42] <- 'PER002'; rownames(LRV.aln.dist)[43] <- colnames(LRV.aln.dist)[43] <- 'PER207'; rownames(LRV.aln.dist)[44] <- colnames(LRV.aln.dist)[44] <- 'PER212'; rownames(LRV.aln.dist)[45] <- colnames(LRV.aln.dist)[45] <- 'PER201';

## RE-DO Fix order of row/colnames:
new_order <- sort(rownames(LRV.aln.dist))
LRV.aln.dist <- LRV.aln.dist[new_order,new_order]

### Leish:
Leish.aln <- read.dna(here::here('3_Cophylo_analysis/1_Preparing_MSA/2_Species/1_Viannia_SNP_alignment/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.Noart.NoNA.vcf.gz.filtered.fa'), format = 'fasta')
Leish.aln.dist <- dist.dna(Leish.aln, model = 'K80', as.matrix = T)
Leish.aln.dist <- Leish.aln.dist[which(rownames(Leish.aln.dist) %in% rownames(link.data)),which(colnames(Leish.aln.dist) %in% rownames(link.data))]

## Fix order of row/colnames:
new_order <- sort(rownames(Leish.aln.dist))
Leish.aln.dist <- Leish.aln.dist[new_order,new_order]

### Adding geographic data:
info <- read_xlsx(here::here('3_Cophylo_analysis/3_Genetic_Dist_Regression/coordinates_for_MMRR.xlsx'), sheet=1)
info <- info[info$`Analysis Code` %in% rownames(Leish.aln.dist),]
#info$StateLat <- jitter(as.numeric(info$StateLat), 10)
#info$StateLon <- jitter(as.numeric(info$StateLon), 10)

### Creating geographic distance matrix:
Coords <- info[,c(6,5)]; Coords <- as.data.frame(Coords)
Coords$Longitude <- as.numeric(unlist(Coords$Longitude)); Coords$Latitude <- as.numeric(unlist(Coords$Latitude))
rownames(Coords) <- info$`Analysis Code`
geograph.dist <- geodist(Coords, measure = 'haversine')
geograph.mtrx <- as.matrix(geograph.dist)
rownames(geograph.mtrx) <- rownames(Coords); colnames(geograph.mtrx) <- rownames(Coords)
geograph.mtrx <- geograph.mtrx/1000 #convert to km's
#fviz_dist(as.dist(geograph.mtrx), lab_size = 8)

### MMRR on standardized distance matrices:

### 2.1 Viannia vs. LRV
#######################
## Creating a list of independent matrices:
X <- list()
X$std.Leish <- Leish.aln.dist

## MMRR with genetic distance as response variable
lrv.leish.mmrr <- MMRR(LRV.aln.dist, X, nperm = 1000) #R2 = 0.65
lrv.leish.mmrr


### 2.2 Geographic distance vs. LRV1 genetic distance
#####################################################
##Creating a list of independent matrices:
X <- list()
X$geo <- geograph.mtrx

## MMRR with genetic distance as response variable
lrv.geo.mmrr <- MMRR(LRV.aln.dist, X, nperm = 1000) 
lrv.geo.mmrr


### 3. Population level: 
### L. braziliensis / L. guyanensis vs. LRV1
############################################

### 3.1 L. braziliensis vs. LRV1
################################
X <- list()
X$std.Lbraz <- Leish.aln.dist[c(1:4,8,9,16,22,23,25,28:45),c(1:4,8,9,16,22,23,25,28:45)]
lrv.leish.mmrr <- MMRR(LRV.aln.dist[c(1:4,8,9,16,22,23,25,28:45),c(1:4,8,9,16,22,23,25,28:45)], X, nperm = 1000) #R2 = 0.65
lrv.leish.mmrr

### 3.2 L. guyanensis vs. LRV1
##############################
X <- list()
X$std.Lguy <- Leish.aln.dist[c(5:7,10:15,17,18,24,26,27),c(5:7,10:15,17,18,24,26,27)]
lrv.leish.mmrr <- MMRR(LRV.aln.dist[c(5:7,10:15,17,18,24,26,27),c(5:7,10:15,17,18,24,26,27)], X, nperm = 1000) 
lrv.leish.mmrr

### EOF ------------------------------------------------------------------------