'CpG Content Correction
Uses MeDEStrand to correct for CpG content and calculate absolute methylation.

Usage:
    correct_cpg.R -b BAM -o OUTPUT

Options:
    -b --bam BAM        Path to input BAM file
    -o --output OUTPUT  Path to output WIG file
' -> doc

library(docopt)
args <- docopt(doc, version='BAM to WIG 1.1')
print(args)

library(MeDEStrand)
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

BSgenome="BSgenome.Hsapiens.UCSC.hg19"
uniq=1
extend=0
shift=0
ws=300
chr.select=paste0("chr", c(1:22,"X","Y"))
bam_file = args[['bam']]

MeDIP_set <- MeDEStrand.createSet(file=bam_file, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws, chr.select=chr.select, paired=T)
CS <- MeDEStrand.countCG(pattern="CG", refObj=MeDIP_set)
MeDIP_binmethyl <- MeDEStrand.binMethyl(MeDIP_set, CS)

write.table(MeDIP_binmethyl, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file = args[['output']])
