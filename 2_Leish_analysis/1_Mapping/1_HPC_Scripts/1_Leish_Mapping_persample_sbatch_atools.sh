#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=30g
#SBATCH --time=36:00:00
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
module load BioTools/2020a.00-intel-2020a
module load BWA
module load SMALT
module load GATK
module load Java


###############################################
### atools environment -- array job per sample
###############################################
## Start logging:
alog --state start

source <(aenv --data /PATH/TO/SAMPLE/LIST/Example_list.csv --no_sniffer)


#########################################
### Setting the REFERENCE GENOME & INDEX:
#########################################
REF=/PATH/TO/REFERENCE/TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome.fasta

### Indexing the reference genome (HAS TO BE DONE ONLY ONCE):
cd /PATH/TO/REFERENCE/

### Smalt indexing --> Normally this has already been done!!
#smalt index -k 13 -s 2 TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome $REF

### BWA indexing
bwa index $REF

### Samtools indexing
samtools faidx $REF

### INDEX:
IDX=/PATH/TO/REFERENCE/TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome

### Create Dictionary of Reference Genome for GATK
gatk --java-options "-Xmx30G" CreateSequenceDictionary -R $REF -O ${IDX}.dict


#######################################################
### Mapping paired end reads against reference genomes:
#######################################################
### Change to Reads directory
cd /PATH/TO/FASTQs

### input files: forward read NAME_R1.fastq.gz ; reversed read NAME_R2.fastq.gz
smalt map -f sam ${IDX} /PATH/TO/FASTQs/${IN}_R1.fastq.gz /PATH/TO/FASTQs/${IN}_R2.fastq.gz | samtools sort -o ${IN}.sort.bam

### indexing bam file
samtools index ${IN}.sort.bam

### Marking duplicates & Filter mapping quality (Keep MQ > 25) in bam file
gatk --java-options "-Xmx30G"  MarkDuplicates -I ${IN}.sort.bam -O ${IN}.sort.m.bam -M ${IN}.sort.md.metrics.txt
samtools view -q 25 -b -o ${IN}.sort.md.bam ${IN}.sort.m.bam

### Indexing the final bam file
samtools index ${IN}.sort.md.bam

### Removing intermediate files
rm ${IN}.sort.bam
rm ${IN}.sort.bam.bai
rm ${IN}.sort.m.bam

### Calculate read depths & move to separate folder
samtools depth -a ${IN}.sort.md.bam > /PATH/TO/2_DEPTHS/${IN}.sort.md.depths

### Moving bam files and associated files to separate folder:
mkdir -p Bams
mv ${IN}.sort.md.bam /PATH/TO/BAMS/
mv ${IN}.sort.md.bam.bai /PATH/TO/BAMS/
mv ${IN}.sort.md.metrics.txt /PATH/TO/BAMS/


### AFTER THIS, PROCEED TO VARIANT CALLING --> Leish_Var_Calling_1_persample_sbatch_atools.sh


## End logging of atools:
alog --state end  --exit $?

