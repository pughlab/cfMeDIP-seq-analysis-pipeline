'BAM QC
Uses MEDIPS to output QC metrics

Usage:
    bam_qc.R -b BAM -r RDATA -s SUMMARY

Options:
    -b --bam BAM                Path to input BAM file
    -r --rdata RDATA            Path to output RDATA file
    -s --summary SUMMARY        Path to output QC summary file
' -> doc

library(docopt)
args <- docopt(doc, version='BAM QC 1.1')
print(args)

library(MEDIPS)
packageVersion('MEDIPS')
library(BSgenome.Hsapiens.UCSC.hg19)

BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=0
extend=0
shift=0
ws=300

major_chromosomes <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')

bam_file=args[['bam']]

cr = MEDIPS.seqCoverage(
    file = bam_file,
    pattern="CG",
    BSgenome=BSgenome,
    extend=extend,
    shift=shift,
    uniq=uniq,
    paired=TRUE,
    chr.select=major_chromosomes
)

er = MEDIPS.CpGenrich(
    file=bam_file,
    BSgenome=BSgenome,
    extend=extend,
    shift=shift,
    uniq=uniq,
    paired=TRUE,
    chr.select=major_chromosomes
)

bam_QC_df = data.frame(
    Sample = args[['bam']], 
    numberReads = cr$numberReads,
    percentReadsW0 = (cr$numberReadsWO / cr$numberReads),
    PercentCpGw0 = (table(cr$cov.res)[1] / sum(table(cr$cov.res))),
    PercentCpGgt5 = (sum(table(cr$cov.res)[6:length(table(cr$cov.res))]) / sum(table(cr$cov.res))),
    enrichment.score.relH = er$enrichment.score.relH,
    enrichment.score.GoGe = er$enrichment.score.GoGe
)

save(bam_QC_df, args[['rdata']])

write.table(bam_QC_df, file=args[['summary']], sep="\t", row.names=FALSE)
