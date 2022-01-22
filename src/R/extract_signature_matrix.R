' extract_signature_matrix.R

Extract specific regions specified in a signature into a matrix
of scores.

Usage:
    extract_signature_matrix.R -s SIG ( -g GLOB | -p PATHS ) -c COL -f FILETYPE -o OUTPUT [ --clip CLIP ]

Options:
    -s --sig SIG                Path to signature file
                                    (see https://github.com/pughlab/methylation-signatures-resource).
    -g --glob GLOB              Globstring pointing to data files.
    -p --paths PATHS            A file containing paths of input
                                    one path per line.
    -c --col COL                Name of the column to extract values from
                                    to construct the matrix. e.g. coverage
    -f --filetype FILETYPE      Filetype, e.g. feather, parquet, csv, or tsv
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
additional_signature_cols = colnames(signature)[! colnames(signature) %in% c('chrom', 'chromStart', 'chromEnd')]

if (!is.null(args[['--paths']])) {
    paths <- read_lines(args[['--paths']])
} else {
    paths <- Sys.glob(args[['--glob']])
}

merged_data <- NULL
for (p in paths) {
    message(sprintf('Loading data: %s', p))

    if (args[['--filetype']] == 'tsv') {
        data <- suppressMessages(read_tsv(p))
    } else if (args[['--filetype']] == 'csv') {
        data <- suppressMessages(read_csv(p))
    } else if (args[['--filetype']] == 'feather') {
        data <- read_feather(p)
    } else if (args[['--filetype']] == 'parquet') {
        data <- read_parquet(p)
    } else {
        stop('Input file format (--filetype) must be one of tsv, csv, feather, or parquet.')
    }

    overlap_result <- foverlaps(
        as.data.table(data),
        as.data.table(signature_dt),
        by.x = c('bin_chr', 'bin_start', 'bin_end')
    ) %>% 
        na.omit(cols=c('chromStart', 'bin_start')) %>%
        as_tibble %>%
        .[, c('bin_chr', 'bin_start', 'bin_end', 'chromStart', 'chromEnd', additional_signature_cols, args[['--col']])] %>%
        distinct()

    colnames(overlap_result)[ncol(overlap_result)] <- basename(p)
    if (!is.null(args[['--clip']])) {
        colnames(overlap_result)[ncol(overlap_result)] <- gsub(args[['--clip']], '', colnames(overlap_result)[ncol(overlap_result)])
    }

    if (is.null(merged_data)) {
        merged_data <- overlap_result
    } else {
        merged_data <- full_join(merged_data, overlap_result, by=c('bin_chr', 'bin_start', 'bin_end', additional_signature_cols))
    }
}

write_parquet(merged_data, args[['--output']])
message(sprintf('Output written to %s', args[['--output']]))
