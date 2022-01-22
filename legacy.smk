# ---------------------------------- #
#  Assemble matrices of the results  #
# ---------------------------------- #
#
# This is legacy code gradually being retired and superceded
# by the above "signatures" section.

def get_matrix_targets(wildcards, path_template):
    """Get all targets for constructing a matrix, breaking files into segments.
    """
    col = name_to_col[wildcards.colname]
    total_bins = sum(1 for line in open(get_cohort_config(wildcards.cohort)['bins']))
    start = list(range(1, total_bins, 3000 * 1000))
    end = [n-1 for n in start[1:]] + [total_bins]
    start_end = [str(s) + '_' + str(e) for s, e in zip(start, end)]
    samples = get_cohort_data(wildcards.cohort).sample_name
    targets = expand(
        path_template.format(cohort = wildcards.cohort, colname = wildcards.colname),
        start_end = start_end,
    )
    return(targets)

def matrix_targets(wildcards):
    """Calls get_matrix_targets and assembles the paths"""
    template = path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp_filtered/matrix_{colname}_{{start_end}}.feather'
    return(get_matrix_targets(wildcards, template))

def matrix_stats_targets(wildcards):
    """Calls get_matrix_targets and assembles the paths for matrix summary stats output"""
    template = path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp_stats/matrix_{colname}_{{start_end}}.feather'
    return(get_matrix_targets(wildcards, template))

# Merges together the data from paths for matrix of results.
rule cfmedip_nbglm_matrix_segment:
    input:
        bins = lambda wildcards: get_cohort_config(wildcards.cohort)['bins'],
        paths = path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp/matrix_{colname}.paths.txt'
    output:
        temp(path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp/matrix_{colname}_{start}_{end}.feather')
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    params:
        colnumber = lambda wildcards: name_to_col[wildcards.colname]
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/merge_range.R -b {input.bins} -p {input.paths} -c {params.colnumber} -o {output} -s {wildcards.start} -e {wildcards.end} --clip "_fit_nbglm.feather"'

# Assembles a matrix subsetting to specific bins
rule matrix_binset:
    input:
        bins = lambda wildcards: get_cohort_config(wildcards.cohort)['binsets'][wildcards.binset],
        paths = path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp/matrix_{colname}.paths.txt'
    output:
        path_to_data + '/{cohort}/results/matrix_binset/{cohort}_{binset}_{colname}.feather'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    params:
        colnumber = lambda wildcards: name_to_col[wildcards.colname]
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/merge_range.R -b {input.bins} -p {input.paths} -c {params.colnumber} -o {output} --clip "_fit_nbglm.feather"'

# Filtered out rows (coverage)
rule cfmedip_nbglm_matrix_segment_coverage_filter:
    input:
        path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp/matrix_coverage_{start}_{end}.feather'
    output:
        temp(path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp_filtered/matrix_coverage_{start}_{end}.feather')
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    run:
        with open(output[0], 'w') as outfile:
            with open(input[0], 'r') as infile:
                outfile.write(next(infile))
                for line in infile:
                    row = [int(n) == 0 for n in line.rstrip('\n').split(',')[3:]] # skip chr, start, and end cols
                    if not all(row):
                        outfile.write(line)

# Filtered out rows (posterior)
rule cfmedip_nbglm_matrix_segment_posterior_filter:
    input:
        path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp/matrix_posterior_{start}_{end}.feather'
    output:
        temp(path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp_filtered/matrix_posterior_{start}_{end}.feather')
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    run:
        with open(output[0], 'w') as outfile:
            with open(input[0], 'r') as infile:
                outfile.write(next(infile))
                for line in infile:
                    row = [float(n) < 0.2 for n in line.rstrip('\n').split(',')[3:]]
                    if not all(row):
                        outfile.write(line)

# Summary stats of NBGLM data
rule cfmedip_nbglm_matrix_segment_stats:
    input:
        path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp/matrix_{metric}_{start}_{end}.csv'
    output:
        temp(path_to_data + '/{cohort}/tmp/cfmedip_nbglm_matrix_temp_stats/matrix_{metric}_{start}_{end}.csv')
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    run:
        with open(output[0], 'w') as outfile:
            outfile.write(','.join([
                'bin_chr', 'bin_start', 'bin_end', 'mean', 'median', 'stdev', 'n_nonzero'
            ]) + '\n')
            with open(input[0], 'r') as infile:
                next(infile)
                for line in infile:
                    row = line.rstrip('\n').split(',')
                    locus = row[0:3]
                    values = [float(n) for n in row[3:]]

                    if not all([n == 0 for n in values]):
                        mean = statistics.mean(values)
                        median = statistics.median(values)
                        stdev = statistics.stdev(values)
                        n_nonzero = sum([n > 0 for n in values])
                        outfile.write(','.join(locus + [str(round(n, 5)) for n in [mean, median, stdev, n_nonzero]]) + '\n')

def row_bind_files(pathlist, outpath):
    """ Bind rows of files based on matching column names
    """
    # check all colnames are the same
    colname_list = [pd.read_csv(i, nrows=10).columns for i in pathlist]
    all_same = all([all(x == colname_list[0]) for x in colname_list])
    if not all_same:
        raise Exception('Header rows are not all the same')
    with open(outpath,"wb") as fout:
        # first file:
        with open(pathlist[0], "rb") as f:
            fout.write(f.read())
        # now the rest:    
        for i in range(1,len(pathlist)):
            with open(pathlist[i], "rb") as f:
                next(f) # skip the header
                fout.write(f.read())

# Assembles the final merged matrix using multiple submatrices
rule cfmedip_nbglm_matrix:
    input:
        matrix_targets
    output:
        path_to_data + '/{cohort}/results/results_matrix/{cohort}_{colname}.csv'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    run:
        row_bind_files(input, output[0])

# Assembles the final merged bin stats
rule cfmedip_nbglm_matrix_stats:
    input:
        matrix_stats_targets
    output:
        path_to_data + '/{cohort}/results/results_matrix/{cohort}_{colname}_stats.csv'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    run:
        row_bind_files(input, output[0])

rule tcga_tissue_signature:
    input:
        ancient('/cluster/projects/pughlab/projects/ezhao/data/TCGA/download/methylation/450k/TCGA-COAD')
    output:
        '/cluster/projects/pughlab/projects/ezhao/projects/COMPARISON/output/signatures/tcga_coad_matrix.Rds'
    resources: cpus=1, mem_mb=16000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/methylation_signatures/tcga_colorectal_signature.R -t {input} -o {output}'

##################
## NGSCheckMate ##
##################
## Provides quality control check to ensure all samples are appropriately paired.


# ---------------------------------------- #
#  Generate Wig and BigWig files from Bam  #
# ---------------------------------------- #
#
# This code is not currently part of the core pipeline
# but is an optional add on if you wish to compare the
# core pipeline results against results from MeDIPs.
        
rule bam_to_wig:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        path_to_data + '/{cohort}/results/bam_to_wig/{sample}_counts.wig',
    resources: mem_mb=30000, time_min='72:00:00'
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome']['human'],
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        clean(r'''
        Rscript src/bam_to_wig.R -b {input} -o {output} -g {params.bsgenome}
        ''')

rule bam_to_bigwig:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        path_to_data + '/{cohort}/results/bigwig/{sample}.aligned.sorted.markdup.bam.{metric}.bigWig',
    conda: 'conda_env/deeptools.yml'
    params:
        blacklist = lambda wildcards: get_cohort_config(wildcards.cohort)['blacklist'],
        normalize = lambda wildcards: ' --normalizeUsing RPKM --ignoreForNormalization chrX chrM Arabidopsis1 Arabidopsis3 ' if wildcards.metric == 'fpkm' else ''
    shell:
        clean('''
        bamCoverage
            --bam {input}
            --outFileName {output}
            --outFileFormat bigwig
            --binSize 300
            --minMappingQuality 30
            --blackListFileName {params.blacklist}
            --numberOfProcessors "max/2"
            --maxFragmentLength 1000
            --samFlagInclude 64
            --extendReads
            {params.normalize}
        ''')

rule bigwig_to_matrix:
    input:
        lambda wildcards: expand(
                path_to_data + '/{{cohort}}/results/bigwig/{sample}.aligned.sorted.markdup.bam.{{metric}}.bigWig',
                sample = get_all_samples(wildcards.cohort).sample_name.unique()
            )
    output:
        path_to_data + '/{cohort}/results/results_matrix/{cohort}_{metric}.npz'
    conda: 'conda_env/deeptools.yml'
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    shell:
        'multiBigwigSummary bins -b {input} -o {output} --smartLabels'


