' bin_stats v1.0
Compute BAM stats by bin

Usage:
    bin_stats.R -b BAM -g GENOME -c CHROM -o OUTPUT [ --bsgchr BSGCHR --filtered FILTERED]

Options:
    -b --bam BAM        Path to input BAM file
    -g --genome GENOME  Either the name of BSgenome (usually BSgenome.Hsapiens.UCSC.hg38
                            or BSgenome.Athaliana.TAIR.TAIR9), or the path to a folder
                            containing a custom BSgenome as a package, which will be loaded
                            using devtools::load_all()

    -c --chrom CHROM    Chromosome
    -o --output OUTPUT  Output path

    --bsgchr BSGCHR     If CHROM does not match the corresponding chromosome name in BSgenome,
                            provide here the actual name of the corresponding BSgenome chromosome
    --filtered FILTERED Path to dump reads that were filted out of the analysis

' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Bin Stats vs. 1.0')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

library(arrow)
library(GenomicRanges)
library(BSgenome)
library(Rsamtools)
library(tidyverse)

BIN_WIDTH = 300
FRAGMENT_LENGTH_LIMIT = 500

bam_file <- BamFile(args[['bam']], asMates=TRUE)
chrom_length <- bam_file %>% seqinfo %>% seqlengths %>% .[args[['chrom']]]

if (file.exists(paste(args[['genome']], 'DESCRIPTION', sep='/'))) {
    devtools::load_all(args[['genome']])
    bsgenome <- getBSgenome(basename(args[['genome']]))
} else {
    bsgenome <- getBSgenome(args[['genome']])
}

if (!is.null(args[['bsgchr']])) {
    stopifnot(args[['bsgchr']] %in% seqnames(bsgenome))
    bsgenome_chr = args[['bsgchr']]
} else {
    bsgenome_chr = args[['chrom']]
}

bins = GRanges(
    seqnames = args[['chrom']],
    ranges = IRanges(
        start = seq(1, chrom_length, BIN_WIDTH),
        end = c(seq(BIN_WIDTH, chrom_length, BIN_WIDTH), chrom_length)
    )
)

bam_reads <- scanBam(
    bam_file,
    param = ScanBamParam(
        which = GRanges(sprintf('%s:1-%s', args[['chrom']], chrom_length)),
        what = scanBamWhat()
    )
)

bam_reads_tibble <- tibble(
    seqnames = args[['chrom']],
    start = bam_reads[[1]]$pos,
    width = bam_reads[[1]]$qwidth,
    mate_status = bam_reads[[1]]$mate_status,
    mateid = bam_reads[[1]]$groupid,
    strand = bam_reads[[1]]$strand,
    mapq = bam_reads[[1]]$mapq
) %>%
    mutate(end = start + width)

fragments_tibble <- bam_reads_tibble %>%
    filter(mapq >= 20) %>%
    group_by(seqnames, mateid) %>%
    filter(n() == 2, '-' %in% strand, '+' %in% strand) %>% 
    filter(mate_status == 'mated') %>%
    summarise(
        paired_end_reads = sprintf(
            '+:%s-%s,-:%s-%s',
            start[strand=='+'],
            end[strand=='+'],
            start[strand=='-'],
            end[strand=='-']
        ),
        start = start[strand == '+'],
        end = end[strand == '-'],
        width = end-start,
        mean_mapq = mean(mapq)
    ) %>% ungroup

fragments_tibble_limited <- fragments_tibble %>% filter(width < FRAGMENT_LENGTH_LIMIT, width > 0)

filtered_reads <- bam_reads_tibble %>% filter(! mateid %in% fragments_tibble_limited$mateid)

message(
    sprintf('Filtered out %s fragments with insufficient quality, length <= 0, or length > %s',
    nrow(fragments_tibble) - nrow(fragments_tibble_limited),
    FRAGMENT_LENGTH_LIMIT
    )
)

fragment_bin_overlaps <- findOverlaps(GRanges(fragments_tibble_limited), bins)

bin_coverage <- bind_cols(
        fragments_tibble_limited[queryHits(fragment_bin_overlaps), ],
        bins[subjectHits(fragment_bin_overlaps)] %>%
            as_tibble %>%
            rename(bin_chr = seqnames, bin_start = start, bin_end = end) %>%
            select(-width)
    ) %>%
    mutate(
        overlap_length = pmin(end, bin_end) - pmax(start, bin_start) + 1
    ) %>%
    group_by(bin_chr, bin_start, bin_end) %>%
    summarise(
        n_fragments = n(),
        coverage_bp = sum(overlap_length),
        mean_fragment_length = mean(width),
        mean_fragment_mapq = mean(mean_mapq)
    ) %>%
    ungroup() %>%
    right_join(
        bins %>%
            as_tibble %>%
            select(seqnames, start, end),
        by = c('bin_chr' = 'seqnames', 'bin_start' = 'start', 'bin_end' = 'end')
    ) %>%
    replace_na(list(
        n_fragments = 0,
        coverage_bp = 0
    )) %>%
    arrange(bin_chr, bin_start) %>%
    mutate(mean_coverage = coverage_bp / BIN_WIDTH) %>%
    mutate(
        seq = getSeq(
            bsgenome,
            names = rep(bsgenome_chr, n()),
            start = bin_start,
            end = bin_end
        ) %>% as.character
    ) %>%
    mutate(
        known_bps = str_count(seq, '[TCGA]'),
        gc_content = str_count(seq, '[GC]') / ifelse(known_bps == 0.0, 1.0, known_bps),
        cpg_count = str_count(seq, 'CG'),
    ) %>%
    select(-known_bps, -seq)

bin_coverage %>%
    write_tsv(
        args[['output']],
        append=TRUE,
        col_names = TRUE
    )

if (!is.null(args[['filtered']])) {
    filtered_reads %>% write_tsv(args[['filtered']])
}
