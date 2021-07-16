' cfmedip_to_array_sites.R
Strip down cfMeDIP-seq data to only the array sites specified in a manifest.

Usage:
    cfmedip_to_array_sites.R -c CFMEDIP -m MANIFEST -o OUTPUT

Options:
    -c --cfmedip CFMEDIP        Path to cfMeDIP-seq data output from cfmedip_nbglm.R
    -m --manifest MANIFEST      Manifest files for the appropriate array, such as can be obtained
                                    from https://zwdzwd.github.io/InfiniumAnnotation in txt format
    -o --output OUTPUT          Path to output cfMeDIP-seq data filtered to array coordinates
                                    with probe data annotated.
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='cfMeDIP-seq to array sites')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

library(tidyverse)
library(data.table)

cfmedip_data <- fread(args[['cfmedip']])
manifest <- fread(args[['manifest']]) %>% filter(!is.na(CpG_beg))

setkey(cfmedip_data, bin_chr, bin_start, bin_end)
setkey(manifest, CpG_chrm, CpG_beg, CpG_end)

overlaps <- foverlaps(
    cfmedip_data,
    manifest %>%
        mutate(
            probe_position = CpG_beg,
            CpG_beg = CpG_beg-50,
            CpG_end = CpG_end+50
        )
    )

overlaps_out <- unique(overlaps[
    !is.na(CpG_beg),
    c('bin_chr',
      'bin_start',
      'bin_end',
      'probeID',
      'probe_position',
      'methylation_status', 
      'coverage', 
      'cpg_count', 
      'gc_content', 
      'methylated_posterior', 
      'methylated_mu', 
      'mean_fragment_length'
    )])

fwrite(overlaps_out, args[['output']])
