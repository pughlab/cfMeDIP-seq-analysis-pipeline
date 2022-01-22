'BAM to WIG file
Uses MEDIPS to export a WIG file from BAM for MeDIP-seq analysis.

Usage:
    bam_to_wig.R -b BAM -o OUTPUT [ -g GENOME ]

Options:
    -b --bam BAM        Path to input BAM file
    -o --output OUTPUT  Path to output WIG file
    -g --genome GENOME  BSgenome ID. If not specified, will use BSgenome.Hsapiens.UCSC.hg38
' -> doc

library(docopt)
args <- docopt(doc, version='BAM to WIG 1.1')
print(args)

library(MEDIPS)
packageVersion('MEDIPS')

if (is.null(args[['genome']])) {
    BSgenome='BSgenome.Hsapiens.UCSC.hg38'
} else {
    BSgenome = args[['genome']]
}

uniq=0
extend=0
shift=0
ws=300

major_chromosomes <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')
bam_MSet <- MEDIPS.createSet(file=args[['bam']], BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=major_chromosomes, paired=TRUE)

MEDIPS.exportWIG(Set=bam_MSet, file=args[['output']], format='count', descr=args[['bam']])
