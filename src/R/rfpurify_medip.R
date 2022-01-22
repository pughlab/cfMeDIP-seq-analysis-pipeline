' rfpurify_medip.R

Apply RFpurify to cfMeDIP-seq data.

Usage:
    rfpurify_medip.R -s SIG ( -g GLOB | -p PATHS ) -o OUTPUT [ --clip CLIP ]

Options:
    -s --sig SIG                Path to signature file
                                    (see https://github.com/pughlab/methylation-signatures-resource).
    -g --glob GLOB              Globstring pointing to data files.
    -p --paths PATHS            A file containing paths of input
                                    one path per line.
    -o --output OUTPUT          Path to output (parquet file format)
    --clip CLIP                 Text to clip from the file basename.
                                    Will be substituted with blank.
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Extract signatures v1.0')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

library(tidyverse)
library(arrow)
library(data.table)

message('Reading signature')
signature <- read_tsv(
    args[['--sig']],
    col_types = cols(
        chrom=col_character(),
        chromStart=col_integer(),
        chromEnd=col_integer()
    )
)
signature_dt <- as.data.table(signature)
setkey(signature_dt, chrom, chromStart, chromEnd)

if (!is.null(args[['paths']])) {
    paths <- read_lines(args[['--paths']])
} else {
    paths <- Sys.glob(args[['--glob']])
}

merged_data <- NULL
for (p in paths) {
    message(sprintf('Loading data: %s', p))
    data <- read_feather(p)
    overlap_result <- foverlaps(
        as.data.table(data),
        as.data.table(signature_dt),
        by.x = c('bin_chr', 'bin_start', 'bin_end')
    ) %>% 
        na.omit(cols=c('chromStart', 'bin_start')) %>%
        as_tibble %>%
        .[, c('bin_chr', 'bin_start', 'bin_end', colnames(signature)[4:ncol(signature)], 'methylated_posterior')] %>%
        distinct()

    colnames(overlap_result)[ncol(overlap_result)] <- basename(p)
    if (!is.null(args[['--clip']])) {
        colnames(overlap_result)[ncol(overlap_result)] <- gsub(args[['--clip']], '', colnames(overlap_result)[ncol(overlap_result)])
    }

    if (is.null(merged_data)) {
        merged_data <- overlap_result
    } else {
        merged_data <- full_join(merged_data, overlap_result, by=c('bin_chr', 'bin_start', 'bin_end', colnames(signature)[4:ncol(signature)]))
    }
}

write_parquet(merged_data, args[['--output']])
message(sprintf('Output written to %s', args[['--output']]))
