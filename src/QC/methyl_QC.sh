#!/bin/bash

SEQMETH="F19K16"
SEQUMETH="F24B22"

BAM_PATH=$1
METH_COUNT=$2
METH_QC=$3

samtools view ${BAM_PATH} | \
    cut -f 3 | sort | uniq -c | sort -nr | \
    sed -e 's/^ *//;s/ /\t/' | \
    awk 'OFS="\t" {print $2,$1}' | \
    sort -n -k1,1 \
    > ${METH_COUNT}

total=$(samtools view -c ${BAM_PATH})
unmap=$(cat ${METH_COUNT} | grep '^\*' | cut -f2); if [[ -z $unmap ]]; then unmap="0"; fi
methyl=$(cat ${METH_COUNT} | grep ${SEQMETH} | cut -f2); if [[ -z $methyl ]]; then methyl="0"; fi
unmeth=$(cat ${METH_COUNT} | grep ${SEQUMETH} | cut -f2); if [[ -z $unmeth ]]; then unmeth="0"; fi
pct_meth_ctrl=$(echo "scale=5; ($methyl + $unmeth)/$total * 100" | bc -l); if [[ -z $pct_meth_ctrl ]]; then pct_meth_ctrl="0"; fi
pct_bac_meth=$(echo "scale=5; $methyl/($methyl + $unmeth) * 100" | bc -l); if [[ -z $pct_bac_meth ]]; then pct_bac_meth="0"; fi
pct_bac_unmeth=$(echo "scale=5; $unmeth/($methyl + $unmeth) * 100" | bc -l); if [[ -z $pct_bac_unmeth ]]; then pct_bac_unmeth="0"; fi
methyl_specificity=$(echo "scale=5; 100 - $pct_bac_unmeth/$pct_bac_meth" | bc -l); if [[ -z $methyl_specificity ]]; then methyl_specificity="0"; fi
bet_meth_ctrl=$(echo "scale=5; $methyl/($methyl + $unmeth)" | bc -l); if [[ -z $bet_meth_ctrl ]]; then bet_meth_ctrl="0"; fi

echo -e "Total_reads\tUnmapped_reads\tNum_reads_align_methyl_F19K16\tNum_reads_align_unmethyl_F24B22\tPctg_reads_aligned_to_BACs\tPctg_BACs_is_methyl_F19K16\tPctg_BACs_is_unmethyl_F24B22\tMethylation_specificity\tMethylation_beta" > ${METH_QC}
echo -e "$total\t$unmap\t$methyl\t$unmeth\t$pct_meth_ctrl\t$pct_bac_meth\t$pct_bac_unmeth\t$methyl_specificity\t$bet_meth_ctrl" >> ${METH_QC}

cat ${METH_QC} \
    | awk -F "\t" '{for (i=1; i<=NF; i++) a[i,NR]=$i
                    max=(max<NF?NF:max)}
                    END {for (i=1; i<=max; i++)
                            {for (j=1; j<=NR; j++) 
                             printf "%s%s", a[i,j], (j==NR?RS:FS)
                            }
                        }' \
    | awk 'BEGIN {OFS="\t"} $1="F19K16_F24B22_methyl_qc\t"$1' \
    > ${METH_QC}.parsed

