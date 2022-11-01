' row_bind_tables.R

Usage: row_bind_tables.R ( -p PATHS | -i INPUT | -G GLOB ) ( --in-csv | --in-tsv | --in-feather | --in-parquet) (--out-csv | --out-tsv | --out-feather | --out-parquet ) -o OUTPUT [ --omit-paths --trim-paths --index-col-name ICN ]

Options:
    -p --paths PATHS            Comma separated list of paths to TSV files
    -i --input INPUT            Path to a file containing paths to tables, one per line
    -G --glob GLOB              Globstring matching paths

    -o --output OUTPUT          Path to output

    --omit-paths                Exclude column listing the path to file
    --trim-paths                Includes only the file name instead of full path
    --index-col-name ICN        Name of the index column containing the paths. Defaults to "paths"
' -> doc

library(docopt)

args <- docopt(doc)

library(plyr)
library(readr)
library(doParallel)
library(arrow)

if ( ! is.null(args[['paths']]) ) {
    paths = strsplit(args[['paths']], ',')[[1]]
} else if (! is.null(args[['input']])) {
    paths = readLines(args[['input']])
} else if ( ! is.null(args[['glob']]) ) {
    paths = Sys.glob(args[['glob']])
} else {
    stop('Must provide one of --paths, --input, or --glob.')
}

paths <- unique(paths)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 5)

message('Merging files')

registerDoParallel()

options(readr.show_progress = FALSE)

output <- ddply(data.frame(paths), 'paths', function(z) {
      path = as.character(z$paths)
      message(paste0('Reading file: ', path))

      if (args[['--in-tsv']]) {
          infile=suppressMessages(read_tsv(path))
      } else if (args[['--in-csv']]) {
          infile=suppressMessages(read_csv(path))
      } else if (args[['--in-feather']]) {
          infile=read_feather(path)
      } else if (args[['--in-parquet']]) {
          infile=read_parquet(path)
      } else {
          stop('Must provide an input file format')
      }
      return(infile)
}, .parallel = TRUE)

if (args[['--trim-paths']]) {
    output <- dplyr::mutate(output,
        paths = sapply(as.character(paths), function(p) { tail(strsplit(p, '/')[[1]], 1) })
    )
}

if (args[['--omit-paths']]) {
    output <- dplyr::select(output, -paths)
}

if (is.null(args[['--index-col-name']])) {
    icn = 'paths'
} else {
    icn = as.character(args[['--index-col-name']])
}
colnames(output)[colnames(output) == 'paths'] <- icn

message('Done merging')

if (args[['--out-tsv']]) {
    write_tsv(output, args[['output']])
} else if (args[['--out-csv']]) {
    write_csv(output, args[['output']])
} else if (args[['--out-feather']]) {
    write_feather(output, args[['output']])
} else if (args[['--out-parquet']]) {
    write_feather(output, args[['output']])
} else {
    stop('Must provide an output file format')
}

message(paste0('Output written to ', args[['output']]))
