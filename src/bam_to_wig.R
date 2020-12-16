'BAM to WIG file
Uses MEDIPS to export a WIG file from BAM for MeDIP-seq analysis.

Usage:
    bam_to_wig.R -b BAM -o OUTPUT

Options:
    -b --bam BAM        Path to input BAM file
    -o --output OUTPUT  Path to output WIG file
' -> doc

library(docopt)
args <- docopt(doc, version='BAM to WIG 1.1')
print(args)

library(MEDIPS)
packageVersion('MEDIPS')
library(BSgenome.Hsapiens.UCSC.hg19)

BSgenome='BSgenome.Hsapiens.UCSC.hg19'
uniq=0
extend=0
shift=0
ws=300

major_chromosomes <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
Raligned.hg19.sorted.markdup.bam_MSet <- MEDIPS.createSet(file=args[['bam']], BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=major_chromosomes, paired=TRUE)

MEDIPS.exportWIG(Set=Raligned.hg19.sorted.markdup.bam_MSet, file=args[['output']], format='count', descr=args[['bam']])
