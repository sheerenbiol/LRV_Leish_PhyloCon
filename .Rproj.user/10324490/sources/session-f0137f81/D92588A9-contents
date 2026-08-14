##############################################
##### Functions to read in genotype data #####
##############################################

### Libraries:
##############
library(readxl);library(data.table); library(stringr)

### Functions:
##############
##Function to read in genotype files:
read.geno <- function(file) {
  geno <- fread(file, data.table = F, header = F)[,-1]
  rownames(geno) <- as.character(read.table(paste(file, 'indv', sep='.'))[,1]) #indv names as rownames
  genopos <- read.table(paste(file, 'pos', sep='.'))
  colnames(geno) <- as.character(paste(genopos[,1], genopos[,2], sep=';')) #Chr.pos as colnames
  return(as.data.frame(geno))
}

##Function to convert genotype file to genlight format:
geno2gl <- function(geno) {
  list <- as.list(as.data.frame(t(geno)))
  loci <- colnames(geno)
  positions <- as.character(lapply(str_split(loci,';'), function(x) x[2]))
  chromosomes <- as.character(lapply(str_split(loci,';'), function(x) x[1]))
  gl <- new('genlight', as.list(as.data.frame(t(geno))))
  gl@chromosome <- as.factor(chromosomes)
  gl@position <- as.factor(positions)
  gl@loc.names <- loci
  return(gl)
}

### EOF ----------------------------------------------------------------------