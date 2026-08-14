#####################################################
#####               SNP Counts                  #####
#####################################################

### Libraries:
##############
library(here)
library(FSA); library(ggplot2)
library(tidyr); library(ggforce)
library(ggpubr); library(scico)

### Calling functions:
######################
source(here::here('2_Leish_analysis/5_SNP_Counts/2_R_Scripts/1_Functions.R'))

### ========================================================================

### Read in genotype data:
##########################
### Genotype index data
chromosomes <- read.table(here::here('1_LRV_analysis/1_Reference_Mapping/3_REFERENCE/TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome.fasta.fai'), colClasses = c('character', 'numeric', 'numeric', 'numeric', 'numeric'))
chromosomes <- chromosomes[c(1:35),] # only need for the 35 nuclear chromosomes!

### Genotype data 
geno <- read.geno(here::here('2_Leish_analysis/5_SNP_Counts/1_Conversion_012_Format/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.NoArt.NoNA.012'))
geno <- geno[order(rownames(geno)),]
dim(geno)
## removing -1 (= NA) values (if any):
col_rm <- apply(geno, 2, function(x) all(x !=-1 )); geno <- geno[,col_rm]

## list of individuals:
All.inds <- rownames(geno)
Study.inds <- All.inds[c(1:40)] ## samples from this study
Ref.inds <- All.inds[c(41:88)]  ## reference samples from public data

## Genotype subsets:
study.geno <- geno[c(1:40),]    ## samples from this study
ref.geno <- geno[c(41:88),]     ## reference samples from public data


### Read in metadata for samples:
#################################
info <- read_excel(here::here('1_LRV_analysis/7_Maps/Table.xlsx'), sheet = 1)
info <- info[info$`Analysis Code` %in%  rownames(geno),]
info <- info[order(info$`Analysis Code`),]

### Changing column properties:
info$`Species/Strain (species typing - method?)` <- as.factor(info$`Species/Strain (species typing - method?)`)
levels(info$`Species/Strain (species typing - method?)`) <- c("Lb", "Lb1","Lb2","Lb3","Lb4","Lg","Lla","Lli","Ln","Lpa","Lpe","Lsh","Lspp.")
info$Species <- as.factor(info$Species)

### Metadata subsets:
study.info <- info[info$`Analysis Code` %in%  Study.inds,] ## samples from this study
ref.info <- info[info$`Analysis Code` %in%  Ref.inds,]     ## reference samples from public data

### Setting species colours:
scico::scico(n=12, palette = 'batlowK')
info$sp_cols <- info$Species
levels(info$sp_cols)
levels(info$sp_cols) <- c("#04050A", "#223A50", "#38545F", "#4B6358", "#62714B", "#80803E", "#A99139",  "#EFA575", "#FAAFA1", "#FDBCCB", "#F9CCF9")
###  Order                   Mixture.       Lb1           Lb2          Lb3         Lb4          Lg           Lla            Ln          Lpa          Lpe.         Lsh
 

### Extract the number of heterozygous & homozygous sites:
#########################################################
### All samples:
Nhetero <- apply(geno, 1, function(x) sum(x==1)) ## Number of heterozygous sites
Nhomo <- apply(geno, 1, function(x) sum(x==2))   ## Number of homozygous sites (for the ALT allele)

### Creating a full data frame with SNP counts and metadata
datacounts <- cbind(Nhomo, Nhetero, info)
## Adding sum of heterozygous and homozygous SNPs:
datacounts$NSum <- datacounts$Nhetero + datacounts$Nhomo
table(datacounts$Species)

## Setting species colours as character type:
datacounts$sp_cols <- as.character(datacounts$sp_cols)
## Setting point shapes variable for plotting (For discrimination between reference and study samples):
datacounts$pch <- as.integer(c(rep(16,22),18,rep(16,6),18,rep(16,6),rep(10,48)))


### PLOTS:
##########
### Quick & Dirty:
plot(datacounts$Nhomo/1000, datacounts$Nhetero/1000,  las = 1, col = datacounts$sp_cols, pch = datacounts$pch, cex = 2, xlab = 'Number of homozygous SNPs (x1000)', ylab= 'Number of heterozygous SNPs (x1000)')

### Nice plot:
SNPcount.plot <- ggplot(data = datacounts, aes(x=Nhomo/1000, y=Nhetero/1000, group = Species, colour = Species))+
  geom_point(shape=datacounts$pch, size=5, alpha=0.8) + scale_color_manual(values = levels(as.factor(datacounts$sp_cols))) +
  xlab('Number of homozygous SNPs (x1000)') + ylab('Number of heterozygous SNPs (x1000)')+
  geom_text(data = datacounts[c(23,30),], aes(label = datacounts[c(23,30),]$`Analysis Code`), nudge_x=4)+
  theme_classic()
SNPcount.plot
SNPcount.plot + facet_zoom(xlim=c(0,10),  ylim = c(0,5))

### Source data for figure:
#write.table(datacounts[,c(1,2,5,7,18,20)], file = './SNPCounts_SourceData.csv', dec = ',', sep = ';', quote = F, row.names = F, col.names = T)


### Comparing the number of SNPs between lineages/species:
##########################################################
### Creating subsets per lineage/species:
L1.snps <- subset(datacounts, datacounts$Species == 'Lb1')
L2.snps <- subset(datacounts, datacounts$Species == 'Lb2')
L3.snps <- subset(datacounts, datacounts$Species == 'Lb3')
L4.snps <- subset(datacounts, datacounts$Species == 'Lb4')
Lg.snps <- subset(datacounts, datacounts$Species == 'Lg')
Lla.snps <- subset(datacounts, datacounts$Species == 'Lla')
Lli.snps <- subset(datacounts, datacounts$Species == 'Lli')
Ln.snps <- subset(datacounts, datacounts$Species == 'Ln')
Lpa.snps <- subset(datacounts, datacounts$Species == 'Lpa')
Lpe.snps <- subset(datacounts, datacounts$Species == 'Lpe')
Lsh.snps <- subset(datacounts, datacounts$Species == 'Lsh')
Mix.snps <- subset(datacounts, datacounts$Species == 'Mixed') 

### Median of total No. of SNPs:
median(L1.snps$NSum)
median(L2.snps$NSum)
median(L3.snps$NSum)
median(L4.snps$NSum)
median(Lg.snps$NSum)
median(Lla.snps$NSum)
median(Ln.snps$NSum)
median(Lpa.snps$NSum)
median(Lpe.snps$NSum)
median(Lsh.snps$NSum)
median(Mix.snps$NSum)
## Kruskal-Wallis test and Dunn's post-hoc test for assessing difference in total No. of SNPs:
KW.test <- kruskal.test(datacounts$NSum ~ datacounts$Species); KW.test
Dunn.tst <- dunnTest(NSum ~ Species, dat=datacounts, method = "bh"); Dunn.tst

### Median of heterozygous SNPs:
median(L1.snps$Nhetero)
median(L2.snps$Nhetero)
median(L3.snps$Nhetero)
median(L4.snps$Nhetero)
median(Lg.snps$Nhetero)
median(Lla.snps$Nhetero)
median(Ln.snps$Nhetero)
median(Lpa.snps$Nhetero)
median(Lpe.snps$Nhetero)
median(Lsh.snps$Nhetero)
median(Mix.snps$Nhetero)
## Kruskal-Wallis test and Dunn's post-hoc test for assessing difference in No. of heterozygous SNPs:
KW.test <- kruskal.test(datacounts$Nhetero ~ datacounts$Species); KW.test
Dunn.tst <- dunnTest(Nhetero ~ Species, dat=datacounts, method = "bh"); Dunn.tst # "bonferroni" OR "bh" method

### Median of homozygous SNPs:
median(L1.snps$Nhomo)
median(L2.snps$Nhomo)
median(L3.snps$Nhomo)
median(L4.snps$Nhomo)
median(Lg.snps$Nhomo)
median(Lla.snps$Nhomo)
median(Ln.snps$Nhomo)
median(Lpa.snps$Nhomo)
median(Lpe.snps$Nhomo)
median(Lsh.snps$Nhomo)
median(Mix.snps$Nhomo)
### Kruskal-Wallis test and Dunn's post-hoc test for assessing difference in No. of homozygous SNPs:
KW.test <- kruskal.test(datacounts$Nhomo ~ datacounts$Species); KW.test
Dunn.tst <- dunnTest(Nhomo ~ Species, dat=datacounts, method = "bh"); Dunn.tst # "bonferroni" OR "bh" method


### Alternate Allele Frequency spectra per lineage/species:
###########################################################
### export individual names per subgroup into a txt file:
write.table(rownames(L1.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lb1.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(L2.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lb2.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(L3.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lb3.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(L4.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lb4.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(Lg.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lg.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(Lla.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lla.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(Ln.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Ln.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(Lpa.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lpa.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(Lpe.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lpe.txt'), row.names = F, col.names = F, quote = F)
write.table(rownames(Lsh.snps), file = here::here('2_Leish_analysis/5_SNP_Counts/Lsh.txt'), row.names = F, col.names = F, quote = F)

### text file is made through the following in bash: 
## $ for i in *.txt; do bcftools view -S $i -Oz /PATH/TO/VCF/FILE/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.NoArt.NoNA.vcf.gz | bcftools query -f '%CHROM %POS %AN %AC{0}\n' | awk '{printf "%s %s %f\n",$1,$2,$4/$3}' > Alt.All.Freq.${i}; done
system(for i in *.txt; do bcftools view -S $i -Oz /PATH/TO/VCF/FILE/Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.NoArt.NoNA.vcf.gz | bcftools query -f '%CHROM %POS %AN %AC{0}\n' | awk '{printf "%s %s %f\n",$1,$2,$4/$3}' > Alt.All.Freq.${i}; done)

### Read in the Alternate Allele Freq data per species:
L1.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lb1.txt'), header = F); colnames(L1.AF) <- c('chrom','pos','AF')
L1.AF <- subset(L1.AF, L1.AF$AF > 0.000000)
L2.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lb2.txt'), header = F); colnames(L2.AF) <- c('chrom','pos','AF')
L2.AF <- subset(L2.AF, L2.AF$AF > 0.000000)
L3.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lb3.txt'), header = F); colnames(L3.AF) <- c('chrom','pos','AF')
L3.AF <- subset(L3.AF, L3.AF$AF > 0.000000)
L4.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lb4.txt'), header = F); colnames(L4.AF) <- c('chrom','pos','AF')
L4.AF <- subset(L4.AF, L4.AF$AF > 0.000000)
Lg.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lg.txt'), header = F); colnames(Lg.AF) <- c('chrom','pos','AF')
Lg.AF <- subset(Lg.AF, Lg.AF$AF > 0.000000)
Lla.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lla.txt'), header = F); colnames(Lla.AF) <- c('chrom','pos','AF')
Lla.AF <- subset(Lla.AF, Lla.AF$AF > 0.000000)
Ln.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Ln.txt'), header = F); colnames(Ln.AF) <- c('chrom','pos','AF')
Ln.AF <- subset(Ln.AF, Ln.AF$AF > 0.000000)
Lpa.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lpa.txt'), header = F); colnames(Lpa.AF) <- c('chrom','pos','AF')
Lpa.AF <- subset(Lpa.AF, Lpa.AF$AF > 0.000000)
Lpe.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lpe.txt'), header = F); colnames(Lpe.AF) <- c('chrom','pos','AF')
Lpe.AF <- subset(Lpe.AF, Lpe.AF$AF > 0.000000)
Lsh.AF <- read.table(here::here('2_Leish_analysis/5_SNP_Counts/Alt.All.Freq.Lsh.txt'), header = F); colnames(Lsh.AF) <- c('chrom','pos','AF')
Lsh.AF <- subset(Lsh.AF, Lsh.AF$AF > 0.000000)

## Plot Alternate Allele Frequency per species:
plot.AF.L1  <- ggplot(data = L1.AF,  aes(x=AF, fill="#223A50")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#223A50") + scale_y_sqrt() #+ scale_y_log10()
plot.AF.L2  <- ggplot(data = L2.AF,  aes(x=AF, fill="#38545F")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#38545F")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.L3  <- ggplot(data = L3.AF,  aes(x=AF, fill="#4B6358")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#4B6358")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.L4  <- ggplot(data = L4.AF,  aes(x=AF, fill="#62714B")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#62714B")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.Lg  <- ggplot(data = Lg.AF,  aes(x=AF, fill="#80803E")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#80803E")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.Lla <- ggplot(data = Lla.AF, aes(x=AF, fill="#A99139")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#A99139")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.Ln  <- ggplot(data = Ln.AF,  aes(x=AF, fill="#EFA575")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#EFA575")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.Lpa <- ggplot(data = Lpa.AF, aes(x=AF, fill="#FAAFA1")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#FAAFA1")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.Lpe <- ggplot(data = Lpe.AF, aes(x=AF, fill="#FDBCCB")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#FDBCCB")+ scale_y_sqrt() #+ scale_y_log10()
plot.AF.Lsh <- ggplot(data = Lsh.AF, aes(x=AF, fill="#F9CCF9")) + geom_histogram(binwidth = 0.05, show.legend = FALSE) + scale_fill_manual(values = "#F9CCF9")+ scale_y_sqrt() #+ scale_y_log10()

ggarrange(nrow = 5, ncol = 2, plot.AF.L1, plot.AF.L2, plot.AF.L3, plot.AF.L4, plot.AF.Lg, plot.AF.Lla, plot.AF.Ln, plot.AF.Lpa, plot.AF.Lpe, plot.AF.Lsh)

## Percantage of Alt alleles fixed in the lineage/species:
L1.fixed <- subset(L1.AF, L1.AF$AF == 1)
L2.fixed <- subset(L2.AF, L2.AF$AF == 1)
L3.fixed <- subset(L3.AF, L3.AF$AF == 1)
L4.fixed <- subset(L4.AF, L4.AF$AF == 1) 
Lg.fixed <- subset(Lg.AF, Lg.AF$AF == 1)
Lla.fixed <- subset(Lla.AF, Lla.AF$AF == 1)
Ln.fixed <- subset(Ln.AF, Ln.AF$AF == 1) 
Lpa.fixed <- subset(Lpa.AF, Lpa.AF$AF == 1)
Lpe.fixed <- subset(Lpe.AF, Lpe.AF$AF == 1)
Lsh.fixed <- subset(Lsh.AF, Lsh.AF$AF == 1)

## Percentage of alt alleles that are rare (AF < 0.1)
L1.rare <- subset(L1.AF, L1.AF$AF < 0.1) 
L2.rare <- subset(L2.AF, L2.AF$AF < 0.1)
L3.rare <- subset(L3.AF, L3.AF$AF < 0.1) 
L4.rare <- subset(L4.AF, L4.AF$AF < 0.1) 
Lg.rare <- subset(Lg.AF, Lg.AF$AF < 0.1) 
Lla.rare <- subset(Lla.AF, Lla.AF$AF < 0.1) 
Ln.rare <- subset(Ln.AF, Ln.AF$AF < 0.1) 
Lpa.rare <- subset(Lpa.AF, Lpa.AF$AF < 0.1)
Lpe.rare <- subset(Lpe.AF, Lpe.AF$AF < 0.1) 
Lsh.rare <- subset(Lsh.AF, Lsh.AF$AF < 0.1)

### EOF ----------------------------------------------------------------------