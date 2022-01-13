' tcga_signature.R

Usage:
    tcga_signature.R -t TCGA -o OUTPUT [ -c CORES ]

Options:
    -t --tcga TCGA              Glob for capturing relevant TCGA files
    -o --output OUTPUT          Path to output
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='TCGA-based cancer signature')
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

library(tidyverse)
library(data.table)
library(purrr)
library(ChAMP)

output_path <- args[['output']]
tcga_globstring = args[['tcga']]
tcga_sample_types <- read_tsv('/cluster/home/zhaoe/git/circulating_immunome/analysis/data/manual/tcga_codes/sampleType.tsv')
tcga_paths <- Sys.glob(tcga_globstring)

meta <- tibble(path = tcga_paths) %>%
    mutate(
        sample_barcode = gsub('.*?(TCGA-..-....-...-...-....-..).*', '\\1', path),
        sample_type = gsub('TCGA-..-....-(..).-...-....-..', '\\1', sample_barcode)
    ) %>% 
    left_join(tcga_sample_types, by = c('sample_type' = 'Code'))

data <- meta$path %>% reduce(function(x, y) {
    if (typeof(x) == 'character') {
        sample_barcode = gsub('.*?(TCGA-..-....-...-...-....-..).*', '\\1', x)
        message(sprintf('Reading sample %s', sample_barcode))
        x = fread(x)[, c('Composite Element REF', 'Beta_value')]
        colnames(x) <- c('probe', sample_barcode)
    }
    sample_barcode = gsub('.*?(TCGA-..-....-...-...-....-..).*', '\\1', y)
    message(sprintf('Reading sample %s', sample_barcode))
    y = fread(y)[, c('Composite Element REF', 'Beta_value')]
    colnames(y) <- c('probe', sample_barcode)
    full_join(x, y, by='probe')
})

data_matrix <- data %>% column_to_rownames('probe') %>% as.matrix
# Filter out probes with all NAs
data_matrix <- data_matrix[apply(data_matrix, 1, function(z) {! all(is.na(z))}), ]
data_imputed <- champ.impute(beta = data_matrix, pd = NULL)

dmr_metadata <- tibble(sample_barcode = colnames(data_imputed)) %>%
    left_join(meta, by = 'sample_barcode') %>%
    mutate(
        include = Definition %in% c('Primary Solid Tumor', 'Solid Tissue Normal')
    ) %>%
    filter(include) %>%
    mutate(Definition = factor(Definition))
levels(dmr_metadata$Definition) <- gsub(' ', '.', levels(dmr_metadata$Definition))

dmr_input <- data_imputed[, colnames(data_imputed) %in% dmr_metadata$sample_barcode]

result_dmrcate <- champ.DMR(
    beta = dmr_input,
    pheno = dmr_metadata$Definition,
    compare.group = c('Primary.Solid.Tumor', 'Solid.Tissue.Normal'),
    arraytype = '450K',
    method = 'DMRcate',
    cores = detectCores() - 1
)

result_probelasso <- champ.DMR(
    beta = dmr_input,
    pheno = dmr_metadata$Definition,
    compare.group = c('Primary.Solid.Tumor', 'Solid.Tissue.Normal'),
    arraytype = '450K',
    method = 'ProbeLasso',
    cores = detectCores() - 1
)


tibble(sample_barcode = colnames(data_imputed)) %>%
    left_join(meta, by='sample_barcode')

saveRDS(data, output_path)
