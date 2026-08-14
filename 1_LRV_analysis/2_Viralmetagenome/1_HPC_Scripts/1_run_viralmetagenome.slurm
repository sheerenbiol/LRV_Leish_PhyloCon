#!/bin/bash -l

#SBATCH --partition=HPC_PARTITION
#SBATCH --account="CREDIT_ACCOUNT"
#SBATCH --chdir="/PATH/TO/SOFTWARE/viralgenie"
#SBATCH --mail-type="END,FAIL,TIME_LIMIT"
#SBATCH --mail-user="mailto@me.be"
#SBATCH --nodes="1"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task="1"
#SBATCH --mem-per-cpu="10000M"
#SBATCH --time="72:00:00"

module load Nextflow/23.04.2
module load git

##### DISCLAIMER #####
### At the time of using this pipeline the name was Viralgenie, which changed at a later stage to ViralMetagenome. 
### Available at: https://nf-co.re/viralmetagenome/1.1.1/


# Set Apptainer/Singularity environment variables to define caching and tmp
# directories. These are used during the conversion of Docker images to
# Apptainer/Singularity ones.
export APPTAINER_CACHEDIR="${HPC_SCRATCH}/apptainer/cache"
export APPTAINER_TMPDIR="${HPC_SCRATCH}/apptainer/tmp"
export NXF_APPTAINER_CACHEDIR="${HPC_SCRATCH}/apptainer/nextflow_cache"

####################################
##### Run a nextflow analysis ######
####################################

#### How to run: ####
# sbatch --export=ALL,ARGS='--input ./Samplesheet.csv --outdir SET/PATH/TO/RESULTS' run_viralmetagenome.slurm

# write a function that catches any errors
abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

# activate the appropriate conda environment
#conda activate nf-core

# write some cluster metadata to the error and output files
hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

# switching to leishamania_if_host_annotaiton branch
# 23-02-2023 16:45 still not merged with dev branch
##git checkout 65-select-best-ref
git checkout dev
## git checkout leishamania_if_host_annotaiton
echo '
**** current hash ****
'"$(git log --pretty=format:'%h' -n 1)"'
**********************
'
# run the nextflow analysis, when running with tower make sure to have /bin/bash -l to load the .bashrc
nextflow run main.nf \
    -ansi-log false \
    -resume \
    -profile vsc_calcua \
    -params-file /PATH/TO/params.yaml \
    --project CREDIT_ACCOUNT \
    $ARGS

### vsc_calcua profile is available at: https://nf-co.re/configs/vsc_calcua/
### CREDIT_ACCOUNT should match the --account flag in the header of the script; i.e. the credit account to make computations on your HPC cluster.

# write some cluster metadata to the error and output files 
printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE ***
************
'
