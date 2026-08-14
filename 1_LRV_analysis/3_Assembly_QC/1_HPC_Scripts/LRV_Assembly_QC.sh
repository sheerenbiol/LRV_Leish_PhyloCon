#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=20g
#SBATCH --time=00:10:00
#SBATCH --partition=HPC_PARTITION
#SBATCH -o slurm-%j.out.file
#SBATCH -e slurm-%j.err.file
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mailto@me.be

###  RE-Setting the environment:
module --force purge
module load HPC/2020a
## loading my modules:
module load atools
module load BLAST+/2.10.0-intel-2020a
module load BioTools/2020a.00-intel-2020a
module load SMALT
module load Pilon/1.23-Java-11

#######################
### atools environment:/scratch/antwerpen/207/vsc20717/LRV_Leish_CoDiv/1_GenData/1_Reads_simlinks
#######################
## Start logging:
alog --state start

source <(aenv --data /PATH/TO/SAMPLE/LIST/LRV_assembly_QC.csv --no_sniffer)


#####################################################################
### Mapping Reads to GENOME ASSEMBLY & Performing Pilon Polishing ###
#####################################################################
### Change directory to where the genome assemblies are located:
cd /PATH/TO/ASSEMBLIES/

### Set assembled genome as reference:
REF=/PATH/TO/ASSEMBLIES/${GENOME}.fa

### Create Index files of REF for SMALT:
smalt index -k 13 -s 2 ${GENOME} $REF

### Set index variable:
IDX=/PATH/TO/ASSEMBLIES/${GENOME}

### Mapping reads (unmapped reads from previous step) to genome assembly:
smalt map -f sam ${IDX} ${READ_PATH}/${SAMPLE}_refmap.sort.unmapped.1.fq.gz ${READ_PATH}/${SAMPLE}_refmap.sort.unmapped.2.fq.gz | samtools sort -o ${GENOME}.sort.bam
samtools index -b ${GENOME}.sort.bam

############################
### Polishing with PILON ###
############################
java -Xmx16G -jar /PATH/TO/pilon-1.23.jar --genome $REF --bam ${GENOME}.sort.bam --output ${GENOME}.pilon

### Rename output file:
mv ${GENOME}.pilon.fasta ${GENOME}.pilon.fa

### Remove intermediate files:
rm ${GENOME}.sort.bam
rm ${GENOME}.sort.bam.bai

####################################################
### RE-Mapping Reads to POLISHED GENOME ASSEMBLY ###
####################################################
### Set assembled genome as reference:
REF=/PATH/TO/ASSEMBLIES/${GENOME}.pilon.fa

### Create Index files of REF for SMALT:
smalt index -k 13 -s 2 ${GENOME}.pilon $REF

### Set index variable:
IDX=/PATH/TO/ASSEMBLIES/${GENOME}.pilon

### RE-Mapping:
smalt map -f sam ${IDX} ${READ_PATH}/${SAMPLE}_refmap.sort.unmapped.1.fq.gz ${READ_PATH}/${SAMPLE}_refmap.sort.unmapped.2.fq.gz | samtools sort -o ${GENOME}.pilon.sort.bam
samtools index -b ${GENOME}.pilon.sort.bam

#################################
### Extract alignment statistics:
#################################
### Extracting read depths per position in the genome
samtools depth -a ${GENOME}.pilon.nogap.sort.bam > ${GENOME}.pilon.nogap.sort.depths

### Extracting mapping statistics:
samtools flagstat ${GENOME}.pilon.nogap.sort.bam > ${GENOME}.pilon.nogap.sort.flagstat ## flagstat file


##############################
### SNP Calling of assemblies:
##############################
### indexing the fasta/reference file
samtools faidx $REF 

### SNP Calling
bcftools mpileup -a INFO/AD -f $REF ${GENOME}.pilon.sort.bam | bcftools call -V indels -v -m -O z -o ${GENOME}.pilon.sort.SNP.vcf.gz


## End logging of atools:
alog --state end  --exit $?

### EOF ----------------------------------------------------------------------------------------------