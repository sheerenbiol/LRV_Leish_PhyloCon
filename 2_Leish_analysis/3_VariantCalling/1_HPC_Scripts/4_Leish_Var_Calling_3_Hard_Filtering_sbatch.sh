#!/usr/bin/env bash

#SBATCH --account=ap_itg_mpu
#SBATCH --partition=broadwell
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=30g
#SBATCH --time=01:00:00
#SBATCH -o slurm-%j.out.file
#SBATCH -e slurn-%j.err.file
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sheeren@itg.be

###  RE-Setting the environment:
module --force purge
module load calcua/2020a
## loading my modules:
module load BioTools/2020a.00-intel-2020a
#module load BWA
#module load SMALT
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


##########################################
### Hard Filtering of Variant Annotations:
##########################################
### Change directory to the VCF files:
cd /PATH/TO/VCFs


### Filtering based on GATK Recommendations:
############################################
### SNPs:
gatk VariantFiltration -R $REF -V Leish.Viannia.GENO.SNP.vcf.gz --filter-expression 'QD<2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'GATKrecom' -O Leish.Viannia.GENO.SNP.GATK.vcf.gz
## Optional: Variants table
#gatk VariantsToTable -V Leish.Viannia.GENO.SNP.GATK.vcf.gz -F CHROM -F POS -F QUAL -F QD -F MQ -F DP -GF DP -GF GQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O Leish.Viannia.GENO.SNP.GATK.txt

### INDELs:
#gatk VariantFiltration -R $REF -V Leish.Viannia.GENO.INDEL.vcf.gz --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0' --filter-name 'GATKrecom' -O Leish.Viannia.GENO.INDEL.GATK.vcf.gz
## Optional: Variants table
#gatk VariantsToTable -V Leish.Viannia.GENO.INDEL.GATK.vcf.gz -F CHROM -F POS -F QUAL -F QD -F MQ -F DP -GF DP -GF GQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O Leish.Viannia.GENO.INDEL.GATK.txt


### Filtering on Format GQ  --> 10 (individual (per sample) Genotype Quality):
##############################################################################
bcftools view -i  'MIN(FMT/GQ)>10' Leish.Viannia.GENO.SNP.GATK.vcf.gz -Oz -o Leish.Viannia.GENO.SNP.GATK.fGQ10.vcf.gz
gatk IndexFeatureFile -I Leish.Viannia.GENO.SNP.GATK.fGQ10.vcf.gz ## to create vcf index file after using bcftools
## Optional: Variants table
#gatk VariantsToTable -V Leish.Viannia.GENO.SNP.GATK.fGQ10.vcf.gz -F CHROM -F POS -F QUAL -F QD -F MQ -F DP -GF DP -GF GQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O Leish.Viannia.GENO.SNP.GATK.fGQ10.txt

### Filtering on Quality --> 100:
#################################
gatk VariantFiltration -R $REF -V Leish.Viannia.GENO.SNP.GATK.fGQ10.vcf.gz --filter-expression 'QUAL < 100 ' --filter-name 'QUAL100' -O Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.vcf.gz
## Optional: Variants table
#gatk VariantsToTable -V Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.vcf.gz -F CHROM -F POS -F QUAL -F QD -F MQ -F DP -GF DP -GF GQ -F FS -F SOR -F MQRankSum -F ReadPosRankSum -O Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.txt


##########################################################################
### Additional filtering of SNPs present on artificial scaffolds, 
### removing missing genotypes, and making sure only sites that passed all 
### filters are retained
##########################################################################

### Extracting only SNPs that passed all filters::
zcat Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.vcf.gz | grep '^#\|PASS' | grep -v 'LbrM.00' | grep -v '\./\.' > Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.vcf
bgzip Leish.Viannia.GENO.SNP.GATK.fGQ10.Q100.PASS.vcf ## BGZIP comes from the HTSLIB package (acompanied with SAMtools and BCFtools) and is used to compress VCF files. 


## EOF ----------------------------------------------------------------------------------------------
