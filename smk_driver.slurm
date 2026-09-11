#!/bin/bash
#SBATCH --job-name=snakemake_driver
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/snakemake_driver_%j.out
#SBATCH --error=logs/snakemake_driver_%j.err

module load GCC/12.3.0 OpenMPI/4.1.5 snakemake/8.4.2 Anaconda3/2024.02-1

# Ensure the log directory for child jobs exists
mkdir -p logs

snakemake \
  --use-conda \
  --use-singularity \
  --conda-frontend conda \
  --jobs 11 \
  --executor cluster-generic \
  --allowed-rules \
    all  \
    dotplot_curated_all \
    dotplot_curated \
    scaffold_with_daedalus \
    mummer_scaffolded \
    dotplot_scaffolded \
  --cluster-generic-submit-cmd \
    "sbatch -A 132714815093 \
            -p medium \
            -t {resources.runtime} \
            --ntasks=1 \
            --mem={resources.mem_mb} \
            --cpus-per-task={threads} \
            --job-name=smk-{rule} \
            --output=logs/{rule}-%j.out \
            --error=logs/{rule}-%j.err" \
  --latency-wait 60 \
  -p
  -n
