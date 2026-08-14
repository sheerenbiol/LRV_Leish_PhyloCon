#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=30g
#SBATCH --time=15:00:00
#SBATCH -o slurm-%j.out.file
#SBATCH -e slurm-%j.err.file
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mailto@me.be

###  RE-Setting the environment:
module --force purge
module load calcua/2020a
## loading my modules:
module load atools
module load BioTools/2020a.00-intel-2020a
#module load BWA
#module load SMALT
module load GATK
module load Java


#######################
### atools environment:
#######################
## Start logging:
alog --state start

source <(aenv --data /PATH/TO/SAMPLE/LIST/Updated_Sample_List.csv --no_sniffer)


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


#################################
### Variant Calling - PER SAMPLE:
#################################
### Go to bam file directory:
cd /PATH/TO/BAMS

### solve multi-sample file error:
gatk AddOrReplaceReadGroups -I ${IN}.sort.md.bam -O ${IN}.sort.md.rg.bam -SORT_ORDER coordinate -RGLB ${IN} -RGPL illumina -RGSM ${IN} -CREATE_INDEX True -RGPU unit1

### HAPLOTYPECALLER: Variant calling with local re-de-novo assembly
gatk HaplotypeCaller -R $REF --min-base-quality-score 25 -I ${IN}.sort.md.rg.bam -O ${IN}.g.vcf.gz -ERC GVCF

### Remove the RRRRG corrected bam file to safe space:
rm ${IN}.sort.md.rg.bam
rm ${IN}.sort.md.rg.bai

### Move g.VCF file to g.vcf folder:
mkdir -p /PATH/TO/GVCFs
mv ${IN}.g.vcf.gz /PATH/TO/GVCFs
mv ${IN}.g.vcf.gz.tbi /PATH/TO/GVCFs

### AFTER THIS, PROCEED WITH Filtering GVCF files by BED files with sufficiently high read depths --> filter_GVCF_by_BED.slrm


## End logging of atools:
alog --state end  --exit $?

### EOF ----------------------------------------------------------------------------------------------
