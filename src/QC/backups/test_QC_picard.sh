#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./
#SBATCH --mem=32G
#SBATCH -J run_QC_picard
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

source /cluster/home/t110409uhn/bin/miniconda3/bin/activate /cluster/home/t110409uhn/git/cfmedip_medremix_bedpe_git/workflow/.snakemake/conda/8ee40725fd397237febfcd4d814d9f27
PICARD_DIR="/cluster/tools/software/picard/2.10.9"

echo "Job started at "$(date) 
time1=$(date +%s)

SAMPLE_NAME="toy01"
BAM_PATH="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis005_toy_QC/QC_picard/toy01_hg38only.bam"

OUT_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis005_toy_QC/QC_picard"
#hg38_BAC_FA="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/reference/hg38_F19K16_F24B22/hg38_F19K16_F24B22.fa"
hg38_FA="/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa"

### CollectGcBiasMetrics
#java -jar ${PICARD_DIR}/picard.jar CollectGcBiasMetrics \
#    I=${BAM_PATH} \
#    O=${OUT_DIR}/${SAMPLE_NAME}.gc_bias_metrics.txt \
#    CHART=${OUT_DIR}/${SAMPLE_NAME}.gc_bias_metrics.pdf \
#    S=${OUT_DIR}/${SAMPLE_NAME}.gc_bias_summary_metrics.txt \
#    R=${hg38_FA}
##    R=${hg38_BAC_FA}
#
### MarkDuplicates
#java -jar ${PICARD_DIR}/picard.jar MarkDuplicates \
#    I=${BAM_PATH} \
#    O=${OUT_DIR}/${SAMPLE_NAME}.picardMarkDup.bam \
#    M=${OUT_DIR}/${SAMPLE_NAME}.marked_dup_metrics.txt
#
### CollectAlignmentSummaryMetrics
#java -jar ${PICARD_DIR}/picard.jar CollectAlignmentSummaryMetrics \
#    I=${BAM_PATH} \
#    O=${OUT_DIR}/${SAMPLE_NAME}.alignmentMetrics.txt \
#    R=${hg38_FA} 
#
### CollectQualityYieldMetrics
#java -jar ${PICARD_DIR}/picard.jar CollectQualityYieldMetrics \
#    I=${BAM_PATH} \
#    O=${OUT_DIR}/${SAMPLE_NAME}.quality_yield_metrics.txt
#
### CollectRrbsMetrics
#java -jar ${PICARD_DIR}/picard.jar CollectRrbsMetrics \
#    I=${BAM_PATH} \
#    M=${OUT_DIR}/${SAMPLE_NAME}.Rrbs_metrics \
#    R=${hg38_FA}






time2=$(date +%s)
echo "Job ended at "$(date)
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"

## EOF
