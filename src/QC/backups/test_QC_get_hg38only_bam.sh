#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./
#SBATCH --mem=32G
#SBATCH -J get_hg38only_bam
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

source /cluster/home/t110409uhn/bin/miniconda3/bin/activate /cluster/home/t110409uhn/git/cfmedip_medremix_bedpe_git/workflow/.snakemake/conda/8ee40725fd397237febfcd4d814d9f27

echo "Job started at "$(date) 
time1=$(date +%s)

SAMPLE_NAME="toy01"
BAM_PATH="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis003_2toys_input_ftype/some_cohort_ftype_fastq.bak/results/bam_markdup/toy01.aligned.sorted.markdup.bam"
OUT_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis005_toy_QC/QC_picard"

samtools view ${BAM_PATH} | grep -v "F19K16\|F24B22" \
    | cat <(samtools view -H ${BAM_PATH} | grep -v "F19K16\|F24B22") - \
    | samtools view -b - \
    > ${OUT_DIR}/${SAMPLE_NAME}_hg38only.bam



time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
