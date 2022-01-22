' run_medestrand.R
Run MeDEStrand.

Usage:
    run_medestrand.R -b BAM -o OUTPUT [ -m MEDESTRAND --bsgenome BSGENOME ]

Options:
    -b --bam BAM                Path to input BAM file
    -o --output OUTPUT          Output path (feather file format)
    -m --medestrand MEDESTRAND  Path to MeDEStrand Package
    --bsgenome BSGENOME         Name of BSgenome package, or path to BSgenome package files
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Running MeDEStrand')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

if (is.null(args[['medestrand']])) {
    library(MeDEStrand)
} else {
    devtools::load_all(args[['medestrand']])
}

if (file.exists(paste(args[['--bsgenome']], 'DESCRIPTION', sep='/'))) {
    devtools::load_all(args[['--bsgenome']])
    BSgenome_name = basename(args[['--bsgenome']])
} else {
    BSgenome_name = args[['--bsgenome']]
}

library(GenomicRanges)
library(tidyverse)
library(arrow)

methylset <- MeDEStrand.createSet(
    file = args[['bam']],
    BSgenome='BSgenome.Hsapiens.UCSC.hg38',
    uniq = 1,
    extend = 200,
    shift = 0,
    window_size = 300,
    chr.select=paste0('chr', c(1:22)),
    paired = T
)

CS = MeDEStrand.countCG(pattern='CG', refObj=methylset)

absolute_methylation = MeDEStrand.binMethyl(MSetInput = methylset, CSet = CS, Granges = TRUE, BSgenome=getBSgenome(BSgenome_name))

absolute_methylation %>%
    as_tibble %>%
    dplyr::select(
        bin_chr = seqnames,
        bin_start = start,
        bin_end = end,
        CF,
        binMethyl
    ) %>%
    write_feather(args[['output']])
