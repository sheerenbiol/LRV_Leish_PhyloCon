#!/usr/bin/env bash

#SBATCH --account=CREDIT_ACCOUNT
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=20g
#SBATCH --time=00:30:00
#SBATCH --partition=HPC_PARTITION
#SBATCH -o stout_orthol-%j.out.file
#SBATCH -e sterr_orthol-%j.err.file
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


###############################################
### atools environment -- array job per sample
###############################################
## Start logging:
alog --state start

source <(aenv --data /PATH/TO/OW_NW_Sample_list.csv --no_sniffer)


###############################
### Trimming raw reads - fastp:
###############################
## Change to reads directory
cd /PATH/TO/FASTQs

### Trimming of the raw reads
fastp -i ${IN}_R1.fastq.gz -I ${IN}_R2.fastq.gz -o ${IN}_trim1_fastq.gz -O ${IN}_trim2_fastq.gz -q 30 -u 10 -5 -3 -W 1 -M 30 --cut_right --cut_right_window_size 10 --cut_right_mean_quality 30 -l 100 -b 150


###############################
### De novo assembly - MEGAHit:
###############################
### Make Assembly folder:
mkdir -p /PATH/TO/ASSEMBLIES

### Perform de novo assembly:
megahit -1 ${IN}_trim1_fastq.gz -2 ${IN}_trim2_fastq.gz --k-list 53 -o /PATH/TO/ASSEMBLIES/${IN}_trim_MH



##########################################
### Extracting Orthologous Genes - BLASTn:
##########################################
### Change to sample specific assembly directory
cd /PATH/TO/ASSEMBLIES/${IN}_trim_MH

### Set variables for Reference genomes of orthologs
NH=/PATH/TO/2_Reference_orthologs/NH.fa
PEIF2AS=/PATH/TO/2_Reference_orthologs/PEIF2AS.fa
PSS1=/PATH/TO/2_Reference_orthologs/PSS1.fa
PTIFAS=/PATH/TO/2_Reference_orthologs/PTIFAS.fa
RPII=/PATH/TO/2_Reference_orthologs/RPII.fa
ZBDH=/PATH/TO/2_Reference_orthologs/ZBDH.fa

### Set variable forQuery genome assemblies - Sample
QUER=/PATH/TO/ASSEMBLIES/${IN}_trim_MH/final.contigs.fa


### Perform BLASTn search
blastn -query $QUER -subject $NH -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tee ${IN}.NH.blast.txt
blastn -query $QUER -subject $PEIF2AS -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tee ${IN}.PEIF2AS.blast.txt
blastn -query $QUER -subject $PSS1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tee ${IN}.PSS1.blast.txt
blastn -query $QUER -subject $PTIFAS -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tee ${IN}.PTIFAS.blast.txt
blastn -query $QUER -subject $RPII -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tee ${IN}.RPII.blast.txt
blastn -query $QUER -subject $ZBDH -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tee ${IN}.ZBDH.blast.txt


### Extract contigs with matches
for i in *.blast.txt; do
    awk '{print $1}' $i | tee ${IN}.$i.contigs_extracted
    for contig in $(cat ${IN}.$i.contigs_extracted); do
        echo $contig
        grep -A 1 $contig /PATH/TO/ASSEMBLIES/${IN}_trim_MH/final.contigs.fa | tee $IN.$i.$contig.fa
    done
done


## End logging:
alog --state end  --exit $?

### EOF ----------------------------------------------------------------------------------------------
