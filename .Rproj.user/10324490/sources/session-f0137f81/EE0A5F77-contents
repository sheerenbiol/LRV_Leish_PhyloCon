#####           Co-phylogenetic analysis          #####
#####            Global-fit statistics:           #####
#######################################################

### Libraries
library(ape); library(vegan); library(phytools)
library(paco); library(tidyverse); library(here)


#############################################
### 1. Old & New World Leishmania & LRV  ####
#############################################

### Data Input:
###############
## Binary matrix with all Leish-LRV links
link.data <- read.table(here::here("3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/Link_data.csv"), header = T, sep=';')
rownames(link.data) <- link.data$X; link.data <- link.data[,-c(1)]
colnames(link.data) <- gsub("\\.", "-", colnames(link.data))
colnames(link.data)[[3]] <- '013-2023-Lb-2018'
link.data <- as.matrix(link.data)


## ML tree of Leish
Leish.tree <- read.tree(here::here("3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/2_LEISH/4_Concatenated_alignment_MLtree/OW-NW_Concat_genes_Nomixes_Aln.fa.treefile"))
tips2drop <- c(Leish.tree$tip.label[!(Leish.tree$tip.label %in% rownames(link.data))])
Leish.tree <- drop.tip(Leish.tree, tip = tips2drop)
plot(Leish.tree)

## ML tree of LRV
LRV.tree <- read.tree(here::here("3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/1_LRV/OW-NW_LRV_Pos_VG_Assembllies_REFS_Concat_genes.fa.treefile"))
LRV.tree$tip.label <- gsub("_", "-", LRV.tree$tip.label) 
LRV.tree$tip.label <- gsub("\\.", "-", LRV.tree$tip.label)
LRV.tree$tip.label <- gsub("\\'", "", LRV.tree$tip.label)
tips2drop <- c(LRV.tree$tip.label[!(LRV.tree$tip.label %in% colnames(link.data))])
LRV.tree <- drop.tip(LRV.tree, tip = tips2drop)
plot(LRV.tree)


### RF - Distance:
##################
## Tips:
tips <- read.table(here::here('3_Cophylo_analysis/1_Preparing_MSA/1_Subgenus/Tipnames.txt'), header = T, sep = '\t')
colnames(tips) <- c('Leish','LRV')
tips$LRV <- as.factor(tips$LRV)
tips$Leish <- as.factor(tips$Leish)
tips <- tips[order(tips$LRV),]
tips <- as.matrix(tips)

## Permutation test:
assoc.tips <- tips
cophyl.stat.perm <- cospeciation(Leish.tree, LRV.tree, distance = 'RF', assoc = assoc.tips, nsim = 1000, method = 'permutation')
cophyl.stat.perm; plot(cophyl.stat.perm) 

### Tanglegram:
###############
par(mar=c(0.5,0.5,0.5,0.5))
cophyl <- cophylo(midpoint_root(Leish.tree), midpoint_root(LRV.tree), assoc = tips, use.edge.length=T, rotate = T)
plot(cophyl, link.type="curved", fsize=0.5, link.lwd=4, link.lty="solid", link.col=make.transparent("red",0.25),use.edge.length=T, scale.bar=c(0.1,0.1))

### ParaFit:
############
H <- cophenetic.phylo(cophyl$trees[[1]])
P <- cophenetic.phylo(cophyl$trees[[2]])
par.fit.1 <- parafit(host.D = H, para.D = P, HP=link.data, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.1

### Export to trees and Link data for eMPRess:
##############################################
write.tree(cophyl$trees[[1]], file = here:here('3_Cophylo_analysis/5_Event_based_cophylo/1_Subgenus/Leish.nwk'))
write.tree(cophyl$trees[[2]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/1_Subgenus/LRV.nwk'))
write.table(cophyl$assoc, here::here('3_Cophylo_analysis/5_Event_based_cophylo/1_Subgenus/tipnames.mapping'), row.names = F, col.names = F, quote = F, sep = ':')
write.table(cophyl$assoc[,c(2,1)], here::here('3_Cophylo_analysis/5_Event_based_cophylo/1_Subgenus/tipnames_rev.mapping'), row.names = F, col.names = F, quote = F, sep = ':')

### PACo:
#########
HP <- link.data; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)

paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000)
paco.fit$gof

paco.r2 <- 1- paco.fit$gof$ss
paco.r2



########################################
### 2. New World Leishmania & LRV1  ####
########################################

### Data Input:
###############
## Binary matrix with all Leish-LRV links
link.data <- read.table(here::here("3_Cophylo_analysis/1_Preparing_MSA/2_Species/Link_data.csv"), header = T, sep=';')
rownames(link.data) <- link.data$X; link.data <- link.data[,-c(1)]
colnames(link.data) <- gsub("\\.", "-", colnames(link.data))
colnames(link.data)[[3]] <- '013-2023-Lb-2018'
link.data <- as.matrix(link.data)

## ML tree of Leish
Leish.tree <- read.tree(here::here("3_Cophylo_analysis/1_Preparing_MSA/2_Species/2_Viannia_phylo/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.Noart.NoNA.vcf.gz.filtered.fa.treefile"))
# Drop tips:
tips2drop <- c(Leish.tree$tip.label[!(Leish.tree$tip.label %in% rownames(link.data))])
Leish.tree <- drop.tip(Leish.tree, tip = tips2drop)
plot(Leish.tree)

## ML tree of LRV
LRV.tree <- read.tree(here::here("1_LRV_analysis/5_ML_Phylogenies/1_IQtree/LRV_WGS_w_REFS/LRV_Pos_VG_Assembllies_REFS.fa.treefile"))
LRV.tree$tip.label <- gsub("_", "-", LRV.tree$tip.label)
LRV.tree$tip.label <- gsub("\\'", "", LRV.tree$tip.label)
LRV.tree$tip.label <- gsub("\\.", "-", LRV.tree$tip.label)
tips2drop <- c(LRV.tree$tip.label[!(LRV.tree$tip.label %in% colnames(link.data))], 'IOCL3574.Ln.2015', 'IOCL3539.Lg.2014')
LRV.tree <- drop.tip(LRV.tree, tip = tips2drop)
plot(LRV.tree)

### RF - Distance:
##################
## Tips:
tips <- read.table(here::here('3_Cophylo_analysis/1_Preparing_MSA/2_Species/Tipnames.txt'), header = T, sep = '\t')
colnames(tips) <- c('Leish','LRV')
tips$LRV <- as.factor(tips$LRV)
tips$Leish <- as.factor(tips$Leish)
tips <- tips[order(tips$LRV),]
tips <- as.matrix(tips)

## Permutation test:
assoc.tips <- tips
cophyl.stat.perm <- cospeciation(Leish.tree, LRV.tree, distance = 'RF', assoc = assoc.tips, nsim = 1000, method = 'permutation')
cophyl.stat.perm; plot(cophyl.stat.perm) 

### Tanglegram:
###############
par(mar=c(0.5,0.5,0.5,0.5))
cophyl <- cophylo(midpoint_root(Leish.tree), midpoint_root(LRV.tree), assoc = tips, use.edge.length=T, rotate = T)
plot(cophyl, link.type="curved", fsize=0.5, link.lwd=4, link.lty="solid", link.col=make.transparent("red",0.25),use.edge.length=T, scale.bar=c(0.1,0.1))

### ParaFit:
############
H <- cophenetic.phylo(cophyl$trees[[1]])
P <- cophenetic.phylo(cophyl$trees[[2]])
par.fit.1 <- parafit(host.D = H, para.D = P, HP=link.data, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.1

### PACo:
#########
HP <- link.data; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)
#paco.input <- prepare_paco_data(H=P, P=H, HP=t(HP))

paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000)
paco.fit$gof

paco.r2 <- 1- paco.fit$gof$ss
paco.r2

### Export to trees for eMPRess:
###########################################################
write.tree(cophyl$trees[[1]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/2_Species/Leish.nwk'))
write.tree(cophyl$trees[[2]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/2_Species/LRV.nwk'))
write.table(cophyl$assoc, here::here('3_Cophylo_analysis/5_Event_based_cophylo/2_Species/tipnames.mapping'), row.names = F, col.names = F, quote = F, sep = ':')
write.table(cophyl$assoc[,c(2,1)], here::here('3_Cophylo_analysis/5_Event_based_cophylo/2_Species/tipnames_rev.mapping'), row.names = F, col.names = F, quote = F, sep = ':')


###################################
### 3. L. braziliensis & LRV1  ####
###################################

### Tip link:
link.data.Lb <- link.data[c(3,4,16,19,22,23,24:45), c(3,4,16,19,22,23,24:45)]
assoc.tips.Lb <- assoc.tips[c(1,5,6,13,20,21,24,27:47),]

### ML tree of Leish
Leish.tree.Lb <- read.tree(here::here("3_Cophylo_analysis/1_Preparing_MSA/2_Species/2_Viannia_phylo/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.Noart.NoNA.vcf.gz.filtered.fa.treefile"))
tips2drop <- c(Leish.tree.Lb$tip.label[!(Leish.tree.Lb$tip.label %in% rownames(link.data.Lb))])
Leish.tree.Lb <- drop.tip(Leish.tree.Lb, tip = tips2drop)
plot(Leish.tree.Lb)

## ML tree of LRV
LRV.tree.Lb <- read.tree(here::here("1_LRV_analysis/5_ML_Phylogenies/1_IQtree/LRV_WGS_w_REFS/LRV_Pos_VG_Assembllies_REFS.fa.treefile"))
LRV.tree.Lb$tip.label <- gsub("_", "-", LRV.tree.Lb$tip.label)
LRV.tree.Lb$tip.label <- gsub("\\'", "", LRV.tree.Lb$tip.label)
LRV.tree.Lb$tip.label <- gsub("\\.", "-", LRV.tree.Lb$tip.label)
tips2drop <- c(LRV.tree.Lb$tip.label[!(LRV.tree.Lb$tip.label %in% colnames(link.data.Lb))])
LRV.tree.Lb <- drop.tip(LRV.tree.Lb, tip = tips2drop)
plot(LRV.tree.Lb)

## RF-distance:
###############
cophyl.Lb.stat.perm <- cospeciation(Leish.tree.Lb, LRV.tree.Lb, distance = 'RF', assoc = assoc.tips.Lb, nsim = 1000, method = 'permutation')
cophyl.Lb.stat.perm; plot(cophyl.Lb.stat.perm) 

## ParaFit:
###########
cophyl.Lb <- cophylo(midpoint_root(Leish.tree.Lb), midpoint_root(LRV.tree.Lb), assoc = assoc.tips.Lb, use.edge.length=T, rotate = T)
plot(cophyl.Lb, link.type="curved", fsize=0.5, link.lwd=4, link.lty="solid", link.col=make.transparent("red",0.25),use.edge.length=T, scale.bar=c(0.01,0.1))

H <- cophenetic.phylo(cophyl.Lb$trees[[1]])
P <- cophenetic.phylo(cophyl.Lb$trees[[2]])
par.fit.1 <- parafit(host.D = H, para.D = P, HP=link.data.Lb, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.1

## PACo:
########
HP <- link.data.Lb; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)
paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000)
paco.fit$gof
paco.r2 <- 1- paco.fit$gof$ss
paco.r2

## Export to trees for eMPRess:
###############################
write.tree(cophyl.Lb$trees[[1]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/1_Lbraziliensis/Leish.nwk'))
write.tree(cophyl.Lb$trees[[2]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/1_Lbraziliensis/LRV.nwk'))
write.table(cophyl.Lb$assoc, here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/1_Lbraziliensis/tipnames.mapping'), row.names = F, col.names = F, quote = F, sep = ':')
write.table(cophyl.Lb$assoc[,c(2,1)], here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/1_Lbraziliensis/tipnames_rev.mapping'), row.names = F, col.names = F, quote = F, sep = ':')


#################################
### 4. L. guyanensis & LRV1  ####
#################################
## Tip link:
link.data.Lg <- link.data[c(2,5:9,11:15,17,18,21), c(2,5:9,11:15,17,18,21)]
assoc.tips.Lg <- assoc.tips[c(2:4,7:12,14,15,19,23,25,26),]

## ML tree of Leish
Leish.tree.Lg <- read.tree(here::here("3_Cophylo_analysis/1_Preparing_MSA/2_Species/2_Viannia_phylo/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.Noart.NoNA.vcf.gz.filtered.fa.treefile"))
tips2drop <- c(Leish.tree.Lg$tip.label[!(Leish.tree.Lg$tip.label %in% rownames(link.data.Lg))])
Leish.tree.Lg <- drop.tip(Leish.tree.Lg, tip = tips2drop)
plot(Leish.tree.Lg)

## ML tree of LRV
LRV.tree.Lg <- read.tree(here::here("1_LRV_analysis/5_ML_Phylogenies/1_IQtree/LRV_WGS_w_REFS/LRV_Pos_VG_Assembllies_REFS.fa.treefile"))
LRV.tree.Lg$tip.label <- gsub("_", "-", LRV.tree.Lg$tip.label)
LRV.tree.Lg$tip.label <- gsub("\\'", "", LRV.tree.Lg$tip.label)
LRV.tree.Lg$tip.label <- gsub("\\.", "-", LRV.tree.Lg$tip.label)
tips2drop <- c(LRV.tree.Lg$tip.label[!(LRV.tree.Lg$tip.label %in% colnames(link.data.Lg))])
LRV.tree.Lg <- drop.tip(LRV.tree.Lg, tip = tips2drop)
plot(LRV.tree.Lg)


## RF-distance:
###############
cophyl.Lg.stat.perm <- cospeciation(Leish.tree.Lg, LRV.tree.Lg, distance = 'RF', assoc = assoc.tips.Lg, nsim = 1000, method = 'permutation')
cophyl.Lg.stat.perm; plot(cophyl.Lg.stat.perm) 

## ParaFit:
###########
cophyl.Lg <- cophylo(midpoint_root(Leish.tree.Lg), midpoint_root(LRV.tree.Lg), assoc = assoc.tips.Lg, use.edge.length=T, rotate = T)
plot(cophyl.Lg, link.type="curved", fsize=0.5, link.lwd=4, link.lty="solid", link.col=make.transparent("red",0.25),use.edge.length=T, scale.bar=c(0.01,0.1))
H <- cophenetic.phylo(cophyl.Lg$trees[[1]])
P <- cophenetic.phylo(cophyl.Lg$trees[[2]])
par.fit.1 <- parafit(host.D = H, para.D = P, HP=link.data.Lg, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.1

## PACo:
########
HP <- link.data.Lg; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)
paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000)
paco.fit$gof
paco.r2 <- 1- paco.fit$gof$ss
paco.r2

## Export to trees for eMPRess:
###############################
write.tree(cophyl.Lg$trees[[1]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/2_Lguyanensis/Leish.nwk'))
write.tree(cophyl.Lg$trees[[2]], file = here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/2_Lguyanensis/LRV.nwk'))
write.table(cophyl.Lg$assoc, here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/2_Lguyanensis/tipnames.mapping'), row.names = F, col.names = F, quote = F, sep = ':')
write.table(cophyl.Lg$assoc[,c(2,1)], here::here('3_Cophylo_analysis/5_Event_based_cophylo/3_Population/2_Lguyanensis/tipnames_rev.mapping'), row.names = F, col.names = F, quote = F, sep = ':')

### EOF ------------------------------------------------------------------------