#!/bin/bash

#SBATCH -J snakemake

#SBATCH --account=mdbf
#SBATCH --qos=short_mdbf
#SBATCH --partition=short_mdbf

#SBATCH --output=/cta/users/cazgari/pipelines/g4/logs/cluster/snakemake-%j.out
#SBATCH --error=/cta/users/cazgari/pipelines/g4/logs/cluster/snakemake-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

eval "$(conda shell.bash hook)"
conda activate snakemake

snakemake --unlock
snakemake --profile config/slurm
