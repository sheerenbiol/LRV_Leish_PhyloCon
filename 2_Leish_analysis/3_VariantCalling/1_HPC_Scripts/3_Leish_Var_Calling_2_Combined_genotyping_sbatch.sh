#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=30g
#SBATCH --time=24:00:00
#SBATCH -o slurm-%j.out.file
#SBATCH -e slurn-%j.err.file
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mailto@me.be

###  RE-Setting the environment:
module --force purge
module load calcua/2020a
## loading my modules:
module load BioTools/2020a.00-intel-2020a
module load BWA
module load SMALT
module load GATK
module load Java


#########################################
### Setting the REFERENCE GENOME & INDEX: SHOULD ALREADY BE DONE IN MAPPING STEP. REMOVE '#' IF NOT.
#########################################
REF=/PATH/TO/REFERENCE/TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome.fasta

### Indexing the reference genome (HAS TO BE DONE ONLY ONCE):
#cd /PATH/TO/REFERENCE/

### Smalt indexing --> Normally this has already been done!!
#smalt index -k 13 -s 2 TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome $REF

### BWA indexing
#bwa index $REF

### Samtools indexing
#samtools faidx $REF

### INDEX:
IDX=/PATH/TO/REFERENCE/TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome

### Create Dictionary of Reference Genome for GATK
#gatk --java-options "-Xmx30G" CreateSequenceDictionary -R $REF -O ${IDX}.dict

### Indexing gVCF files:
cd /PATH/TO/GVCFs
for file in /PATH/TO/GVCFs/*g.vcf.gz ; do bcftools index -t $file; done

########################################
### Variant Calling - COMBINING SAMPLES:
########################################
### Change to VCF directory:
mkdir -p /PATH/TO/VCFs
cd /PATH/TO/VCFs

### Combining all gVCF files into one VCF file:
gatk CombineGVCFs -R $REF -V /PATH/TO/10X_gvcfs.LRV-Pos.list -O Leish.Viannia.combined.vcf.gz

## "10X_gvcfs.LRV-Pos.list" represents a list of all gVCF files from samples with
## a median CDS coverage of at least 10X + only covering CDS regions with at least 10X coverage.
###
## IMPORTANT: This list contains PATHS to the filtered gVCF files and thus should be update to your needs.

### Joint Genotyping of the combined gVCF file:
gatk GenotypeGVCFs -R $REF -V Leish.Viannia.combined.vcf.gz -O Leish.Viannia.GENO.vcf.gz

### Separating SNPs and INDELs from each other:
gatk SelectVariants -R $REF -V Leish.Viannia.GENO.vcf.gz -select-type SNP -O Leish.Viannia.GENO.SNP.vcf.gz
gatk SelectVariants -R $REF -V Leish.Viannia.GENO.vcf.gz -select-type INDEL -O Leish.Viannia.GENO.INDEL.vcf.gz

### -----
### Optional: extract a table with variant statistics (e.g., QUAL, QD, MQ, DP, GQ, FS, SOR, MQRankSum, ReadPosRankSum)
gatk VariantsToTable -V Leish.Viannia.GENO.SNP.vcf.gz -F CHROM -F POS -F QUAL -F QD -F MQ -F DP -GF DP -GF GQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O Leish.Viannia.GENO.SNP.UNfiltered.txt
gatk VariantsToTable -V Leish.Viannia.GENO.INDEL.vcf.gz -F CHROM -F POS -F QUAL -F QD -F MQ -F DP -GF DP -GF GQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O Leish.Viannia.GENO.INDEL.UNfiltered.txt
### -----


### AFTER THIS, PROCEED WITH VARIANT FILTRATION --> Leish_Var_Calling_3_Hard_Filtering_sbatch.sh


### EOF ----------------------------------------------------------------------------------------------