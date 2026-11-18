#!/bin/bash

## Note - Slurm script comments require two hash symbols (##).  A single
## hash symbol immediately followed by SBATCH indicates an SBATCH
## directive.  "##SBATCH" indicates the SBATCH command is commented
## out and is inactive.

## NTasks is not thread count, be generally it's fine to leave it set at 1
##SBATCH --ntasks=1

## CPUS = threads for most purposes. 
#SBATCH --cpus-per-task=1

## The job does not require much memory
#SBATCH --mem=4G

## Time - make sure that the executor job for nextflow
## runs for at least as long as it takes to finish the whole analysis.
#SBATCH --time=24:00:00
## Jobs can restart automatically if the cluster experiences
## an unforeseen issue.  This may not be desired if you want to retain
## the work that's been performed by your script so far. Default is no-requeue  
## --no-requeue # don't restart if it fails
## --requeue # restart if it fails

## Normal Slurm options
## SBATCH -p shared
#SBATCH --job-name="AmoebaME"
#SBATCH --output=AmoebaME_%j.out

## Partition
## usually can be 'normal' but might want large in compute_zone2 or longjobs.
#SBATCH --partition=compute_zone2,normal 

## Load the appropriate modules or software first
module purge # clears any existing software

## For singularity
module load singularity

## For nextflow:
## Run the Nextflow controller as a small Slurm job.  It submits the workflow's
## processes as further Slurm jobs (see nextflow.config) and run longer than any
## of its daughter-jobs; because it's just orchestrating the other slurm jobs, it
## does not need a large allocation of resources.
eval "$(${HOME}/miniconda3/bin/conda shell.bash hook)"
conda activate nextflow

# Set system variables
export GENOME_HOME="$HOME/Genome-Pipeline"
export NXF_CONDA_CACHEDIR="$GENOME_HOME/conda_cache"


# Change me for each input and output you want to run
INPUT="$GENOME_HOME/data/incendi.fastq.gz"
OUTPUT="$HOME/incendi/Genome_results"
mkdir -p $OUTPUT

## Instruct your program to make use of the number of desired threads.
## As your job will be allocated an entire node, this should normally
## be 24.
cd ${GENOME_HOME}
#nextflow run main.nf \
#  -profile slurm \
#  --input_fgz $INPUT \
#  --outdir    $OUTPUT \
#  --db_dir    refdbs \
#  -resume

# Two or more samples
nextflow run main.nf \
    -profile slurm \
    -c "$GENOME_HOME/conf/zest.config" \
    --samplesheet "$HOME/incendi/genome_samplesheet.csv" \
    --outdir "$HOME/incendi/Genome_results" \
    --db_dir "$HOME/Genome-Pipeline/refdbs" \
    -resume
