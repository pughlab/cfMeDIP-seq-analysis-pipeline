#! /bin/bash -login
#SBATCH -J snakemake-submission
#SBATCH -t 5-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eric.zhao@rmp.uhn.ca

echo 'Running on H4H cluster'
snakemake -p --profile slurm
