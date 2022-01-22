' merge_range.R
Merge results into matrix.

Usage:
    merge_range.R -b BINS ( -g GLOB | -p PATHS ) -c COL -o OUTPUT [ -s START -e END --clip CLIP ]

Options:
    -b --bins BINS              File containing all bins.
    -g --glob GLOB              Globstring pointing to data files.
    -p --paths PATHS            A file containing paths of input
                                    one path per line.
    -c --col COL                Column of the data to merge on
    -o --output OUTPUT          Path to output (CSV file)

    -s --start START            If subsetting, the first row in BINS file
                                    to include
    -e --end END                If subsetting, the final row in BINS file
                                    to include
    --clip CLIP                 Text to clip from the file basename.
                                    Will be substituted with blank.
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Bin Stats vs. 1.0')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

library(readr)
library(data.table)
library(arrow)

message('Reading bins')
bins <- read_feather(args[['bins']])
setDT(bins)

if (is.null(args[['start']])) {
    args[['start']] = 1
}
if (is.null(args[['end']])) {
    args[['end']] = nrow(bins)
}
bins <- bins[as.integer(args[['start']]) : as.integer(args[['end']])]
bins <- bins[, c(1,2,3)] # remove cpg count data for now

message('Extracting rows')

if (!is.null(args[['paths']])) {
    paths <- read_lines(args[['paths']])
} else {
    paths <- Sys.glob(args[['glob']])
}

paths <- unique(paths)

for (p in paths) {
    message(p)
    data <- read_feather(p, select=c(1,2,3,as.integer(args[['col']])))
    setDT(data)
    colnames(data)[4] = basename(p)
    if (!is.null(args[['clip']])) {
        colnames(data)[4] <- gsub(args[['clip']], '', colnames(data)[4])
    }
    bins <- merge(bins, data, by=c('bin_chr', 'bin_start', 'bin_end'))
}

write_feather(as.data.frame(bins), args[['output']])
message(sprintf('Output written to %s', args[['output']]))
