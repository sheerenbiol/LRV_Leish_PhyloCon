#########################################################################
#####  Global fit Cophylogenetic assessment WGS/CP/RDRP alignments  #####
#########################################################################

### IMPORTANT:
### MAKE SURE THAT TREE FILES FROM IQTREE ARE IN NEXUS FORMAT!
### CAN BE DONE BY OPENING .TREEFILE FROM IQTREE IN FIGTREE AND EXPORTING AS NEXUS FILE ANS SAVE AS .TREE

### Libraries
library(here)
library(ape); library(paco); library(phytools)

### Read in data:
###################
### Link data (linkage between tips of trees) -- For Parafit and PACo
link.data <- read.table(here::here('1_LRV_analysis/5_ML_Phylogenies/3_Gene_Cophylo/Link_data_ALL.csv'), header = T, sep=';')

### Link data for Robinson-Foulds distance
tips <- read.table(here::here('1_LRV_analysis/5_ML_Phylogenies/3_Gene_Cophylo/Tipnames_ALL.txt'), header = T, sep = '\t')
colnames(tips) <- c('CP','RDRP')
tips$CP <- as.factor(tips$CP)
tips$RDRP <- as.factor(tips$RDRP)
tips <- tips[order(tips$CP),]
tips <- as.matrix(tips)

## Make column and rownames identical:
rownames(link.data) <- link.data$X; link.data <- link.data[,-c(1)]
link.data <- as.matrix(link.data)
colnames(link.data) <- rownames(link.data)

### ML tree - WGS (NEXUS format)
WGS.tree <- read.nexus(here::here("1_LRV_analysis/5_ML_Phylogenies/1_IQtree/LRV_WGS_w_REFS/LRV_Pos_VG_Assembllies_REFS.fa.tree"))
## Dropping out References sequences (NC_003601.1, NC_002063.1) and IOCL1398-Lg-1989 (i.e.,Not present in RDRP).
WGS.tree$tip.label <- gsub("'","",WGS.tree$tip.label); WGS.tree <- drop.tip(WGS.tree, c('NC_003601.1','NC_002063.1','IOCL1398-Lg-1989'))

### ML tree - CP (NEXUS format)
CP.tree <- read.nexus(here::here("1_LRV_analysis/5_ML_Phylogenies/1_IQtree/LRV_CP_w_REFS/LRV_Pos_VG_Assembllies_REFS_CP.fa.tree"))
## Dropping out References sequences (NC_003601.1, NC_002063.1) and IOCL1398-Lg-1989 (i.e.,Not present in RDRP).
CP.tree$tip.label <- gsub("'","",CP.tree$tip.label); CP.tree<- drop.tip(CP.tree, c('NC_003601.1','NC_002063.1','IOCL1398-Lg-1989'))

### ML tree of RDRP (NEXUS format)
RDRP.tree <- read.nexus(here::here("1_LRV_analysis/5_ML_Phylogenies/1_IQtree/LRV_RDRP_w_REFS/LRV_Pos_VG_Assembllies_REFS_RDRP.fa.tree"))
## Dropping out References sequences (NC_003601.1, NC_002063.1) and IOCL1398-Lg-1989 (i.e.,Not present in RDRP).
RDRP.tree$tip.label <- gsub("'","",RDRP.tree$tip.label); RDRP.tree<- drop.tip(RDRP.tree, c('NC_003601.1','NC_002063.1','IOCL1398-Lg-1989'))


### Robinson - Foulds Distance:
###############################
### Permutation test:
assoc.tips <- tips[,c(2,1)] ## Change order of columns.

### CP - RDRP
cophyl1.stat.perm <- cospeciation(RDRP.tree, CP.tree, distance = 'RF', assoc = assoc.tips, nsim = 1000, method = 'permutation')

### WGS - CP
cophyl2.stat.perm <- cospeciation(CP.tree, WGS.tree, distance = 'RF', assoc = assoc.tips, nsim = 1000, method = 'permutation')

### WGS - RDRP
cophyl3.stat.perm <- cospeciation(RDRP.tree, WGS.tree, distance = 'RF', assoc = assoc.tips, nsim = 1000, method = 'permutation')

### Return + plot permutation results:
cophyl1.stat.perm; plot(cophyl1.stat.perm)
cophyl2.stat.perm; plot(cophyl2.stat.perm)
cophyl3.stat.perm; plot(cophyl3.stat.perm)


### ParaFit:
############
### CP - RDRP
H <- cophenetic.phylo(cophyl$trees[[1]])
P <- cophenetic.phylo(cophyl$trees[[2]])
par.fit.1 <- parafit(host.D = H, para.D = P, HP=link.data, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.1

### WGS - CP
H <- cophenetic.phylo(cophyl2$trees[[1]])
P <- cophenetic.phylo(cophyl2$trees[[2]])
par.fit.2 <- parafit(host.D = H, para.D = P, HP=link.data, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.2

### WGS - RDRP
H <- cophenetic.phylo(cophyl3$trees[[1]])
P <- cophenetic.phylo(cophyl3$trees[[2]])
par.fit.3 <- parafit(host.D = H, para.D = P, HP=link.data, nperm = 1000, test.links = T, correction = 'cailliez')
par.fit.3


### PACo:
#########
### CP - RDRP
H <- cophenetic.phylo(cophyl$trees[[1]])
P <- cophenetic.phylo(cophyl$trees[[2]])
HP <- link.data; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)

paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000, shuffled = T)
paco.fit$gof

paco.r2 <- 1- paco.fit$gof$ss
paco.r2

### WGS - CP
H <- cophenetic.phylo(cophyl2$trees[[1]])
P <- cophenetic.phylo(cophyl2$trees[[2]])
HP <- link.data; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)

paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000, shuffled = T)
paco.fit$gof

paco.r2 <- 1- paco.fit$gof$ss
paco.r2

### WGS - RDRP
H <- cophenetic.phylo(cophyl3$trees[[1]])
P <- cophenetic.phylo(cophyl3$trees[[2]])
HP <- link.data; rownames(HP) <-  rownames(H); colnames(HP) <- colnames(P)
paco.input <- prepare_paco_data(H=H, P=P, HP=HP)

paco.Pc <- add_pcoord(paco.input, correction = 'cailliez')
paco.fit <- paco::PACo(D=paco.Pc, symmetric = F, nperm = 1000, shuffled = T)
paco.fit$gof

paco.r2 <- 1- paco.fit$gof$ss
paco.r2


### Tanglegrams:
###############
par(mar=c(0.5,0.5,0.5,0.5))
cophyl <- cophylo(midpoint_root(CP.tree), midpoint_root(RDRP.tree), assoc = tips[,c(2,1)], use.edge.length=T, rotate = T)
plot(cophyl, link.type="curved", fsize=0.5, link.lwd=4, link.lty="solid", link.col=make.transparent("red",0.25),use.edge.length=T, scale.bar=c(0.1,0.1))

cophyl2 <- cophylo(midpoint_root(WGS.tree), midpoint_root(CP.tree), assoc = tips[,c(2,1)], use.edge.length=T, rotate = T)
cophyl3 <- cophylo(midpoint_root(WGS.tree), midpoint_root(RDRP.tree), assoc = tips[,c(2,1)], use.edge.length=T, rotate = T)