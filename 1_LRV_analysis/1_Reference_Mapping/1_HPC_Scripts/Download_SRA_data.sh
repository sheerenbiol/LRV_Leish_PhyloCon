#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=8
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
module load SRA-Toolkit/3.0.5-gompi-2021a


###############################################
### atools environment -- array job per sample
###############################################
## Start logging:
alog --state start

source <(aenv --data /PATH/TO/SAMPLE/LIST/Sample_list.csv)


########################################
### Downloading FASTQ files from SRA ###
########################################

### Downloading/ Pre-fetching the data's SRA-formatted file
prefetch ${ACCESSION}

# 2. Convert to SRA file to FASTQ
fasterq-dump ${ACCESSION} --split-files -e 8 -O ./2_FASTQs

### IMPORTANT:
echo Change fastq file name from ${ACCESSION} to ${IN}, the actual sample name before proceeding with the analyses.


## End logging of atools:
alog --state end  --exit $?

### EOF ----------------------------------------------------------------------------------------------
