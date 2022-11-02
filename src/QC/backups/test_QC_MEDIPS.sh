#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./
#SBATCH --mem=32G
#SBATCH -J run_QC_MEDIPS
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

source /cluster/home/t110409uhn/bin/miniconda3/bin/activate /cluster/home/t110409uhn/git/cfmedip_medremix_bedpe_git/workflow/.snakemake/conda/6f4c279e409bd806d2e93df306a0b542

echo "Job started at "$(date) 
time1=$(date +%s)

#BAM_PATH="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis003_2toys_input_ftype/some_cohort_ftype_fastq/results/bam_markdup/toy01.aligned.sorted.markdup.bam"
BAM_PATH="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/data/CMP-01-02-cfDNA-03.chr22.bam"
OUT_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis005_toy_QC/QC_MEDIPS"
BSGENOME="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/reference/bsgenome/BSgenome.Hsapiens.UCSC.hg38"

Rscript ./QC_MEDIPS.R \
    --bamFile ${BAM_PATH} \
    --outputDir ${OUT_DIR} \
    --windowSize 300

#    --genome ${BSGENOME}






time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
