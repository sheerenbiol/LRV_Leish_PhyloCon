#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=01:00:00
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


###############################################
### atools environment -- array job per sample
###############################################
## Start logging:
alog --state start

source <(aenv --data /PATH/TO/SAMPLE/LIST/Sample_list.csv)


###############################
### Trimming raw reads - fastp:  NEED TO BE DOWNLOADED FIRST!!!
###############################
## Change to reads directory
cd /PATH/TO/FASTQs

## Trimming of raw reads with fastp:
fastp -i ${IN}_R1.fastq.gz -I ${IN}_R2.fastq.gz -o ${IN}_trim1_fastq.gz -O ${IN}_trim2_fastq.gz -q 30 -u 10 -5 -3 -W 1 -M 30 --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 -l 100 -b 150


########################################################################
### Extracting unmapped reads to Leishmania Ref genome for LRV assembly:
########################################################################
### Indexing Reference genome (HAS TO BE DONE ONLY ONCE):
smalt index -k 13 -s 2 TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome /PATH/TO/REFERENCE/TriTrypDB-65_LbraziliensisMHOMBR75M2904_2019_Genome.fasta

### Mapping trimmed reads to Leishmania Ref genome:
smalt map -f sam ${IDX} /PATH/TO/FASTQs/${IN}_R1.fastq.gz /PATH/TO/FASTQs/${IN}_R2.fastq.gz | samtools sort -o ${IN}.sort.bam

### Indexing BAM file:
samtools index -b ${IN}.sort.bam

### Extract unmapped reads:
samtools view -b -f 4 ${IN}.sort.bam > ${IN}.unmapped.bam
samtools sort -n -o ${IN}.unmapped.sort.bam ${IN}.unmapped.bam
samtools index -b ${IN}.unmapped.sort.bam

### Convert BAM to FASTQ files:
bedtools bamtofastq -i /PATH/TO/FASTQs/${IN}.unmapped.sort.bam -fq ${IN}_refmap.sort.unmapped.1.fq -fq2 ${IN}_refmap.sort.unmapped.2.fq

### gzip FASTQ files:
gzip ./*fq

### Remove intermediate files:
rm ${IN}.sort.bam
rm ${IN}.sort.bam.bai
rm ${IN}.unmapped.bam
rm ${IN}.unmapped.sort.bam
rm ${IN}.unmapped.sort.bam.bai


## End logging of atools:
alog --state end  --exit $?

### EOF ----------------------------------------------------------------------------------------------



