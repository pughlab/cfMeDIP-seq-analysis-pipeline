#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./logs_slurm/
#SBATCH --mem=60G
#SBATCH -J step7_getFilterMetrics
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

# getopts ###################################################
usage(){
    echo 
    echo "Usage: bash getFilterMetrics.sh -r BAM_RAW -a BAM_FILT1 -b BAM_FILT2 -c BAM_FILT3 -d BAM_DEDUPD -o FILT_METRICS" 
    echo 
}
no_args="true"

## Help 
Help()
{   
    # Display Help
    echo 
    echo "Get how many reads were removed during each filtering step."
    echo
    echo "Usage: bash getFilterMetrics.sh -r BAM_RAW -a BAM_FILT1 -b BAM_FILT2 -c BAM_FILT3 -d BAM_DEDUPD -o FILT_METRICS" 
    echo "options:"
    echo "-h   [HELP]      print help"
    echo "-r   [REQUIRED]  raw/unfiltered bam (full path)"
    echo "-a   [REQUIRED]  bam with unmapped and secondary reads removed (full path)"
    echo "-b   [REQUIRED]  bam with reads belonging to inserts shorter than 119nt or greater than 501nt removed (full path)"
    echo "-c   [REQUIRED]  bam with reads with edit distance > 7 from reference removed (full path)"
    echo "-d   [REQUIRED]  bam UMI deduplicated (full path)"
    echo "-o   [REQUIRED]  output directory (full path)"
    echo
}

## Get the options
while getopts ":hr:a:b:c:d:o:" option; do
    case "${option}" in
        h) Help
           exit;;
        r) BAM_RAW=${OPTARG};;
        a) BAM_FILT1=${OPTARG};;
        b) BAM_FILT2=${OPTARG};;
        c) BAM_FILT3=${OPTARG};;
        d) BAM_DEDUPD=${OPTARG};;
        o) FILT_METRICS=${OPTARG};;
       \?) echo "Error: Invalid option"
           exit;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }


# Main program ##############################################

echo "Processing getFilterMetrics..." 
echo "Job started at "$(date) 
time1=$(date +%s)

bam_raw=$(samtools view -c ${BAM_RAW})
bam_filter1=$(samtools view -c ${BAM_FILT1})
bam_filter2=$(samtools view -c ${BAM_FILT2})
bam_filter3=$(samtools view -c ${BAM_FILT3})
bam_dedupd=$(samtools view -c ${BAM_DEDUPD})

echo -e "bam_stage\tread_count" > ${FILT_METRICS}
echo -e "bam_raw\t$bam_raw" >> ${FILT_METRICS}
echo -e "bam_filter1\t$bam_filter1" >> ${FILT_METRICS}
echo -e "bam_filter2\t$bam_filter2" >> ${FILT_METRICS}
echo -e "bam_filter3\t$bam_filter3" >> ${FILT_METRICS}
echo -e "bam_dedupd\t$bam_dedupd" >> ${FILT_METRICS}


echo "Finished processing getting filtering metrics."

time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"
echo ""

## EOF
