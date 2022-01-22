import yaml
import pandas as pd
import gzip
import statistics

configfile: "config.yml"

path_to_data = config['data']['base_path']

excluded_cases = config['data']['excluded_cases']
if excluded_cases:
    print('Excluding cases:')
    for case in excluded_cases:
        print(case)
else:
    excluded_cases = []

# ---------------------------- #
#  Important Shared Functions  #
# ---------------------------- #

def get_active_cohorts():
    """Returns a list of active cohorts.
    Active cohorts are defined in config.yml data > cohorts > [cohortname] > active: True
    """
    return(c for c in config['data']['cohorts'] if config['data']['cohorts'][c]['active'])

def get_cohort_config(cohort):
    """Returns the configuration details for specific cohort.
    Importantly, retrieving a cohort using this function (rather than directly
    accessing it via the config dictionary) will also inject the default settings
    (config > data > defaults) wherever a specific setting is not specified.
    """
    config_data = dict(config['data']['defaults'])
    if 'settings' in config['data']['cohorts'][cohort]:
        for key in config['data']['cohorts'][cohort]['settings']:
            if key in config_data:
                config_data[key] = config['data']['cohorts'][cohort]['settings'][key]
            else:
                raise(Exception('Setting {} does not exist for cohort {}. Available settings: {}'.format(
                    key, cohort, ', '.join(config_data.keys())
                )))
    return(config_data)

def get_cohort_data(cohort):
    """Parses the samplesheet for a specific cohort.
    Also removes any excluded_cases from the samplesheet before it is returned.
    """
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort]['samplesheet'], comment='#').drop_duplicates()
    samplesheet = samplesheet[~samplesheet.sample_name.isin(excluded_cases)]
    return(samplesheet)

def get_fastq_path(cohort, sample, library, read_in_pair=1):
    """Retrieves the path to the fastq file.

    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
        library -- integer representing the library index as specified in the samplesheet.
        read_in_pair -- 1 or 2 - representing read 1 or read 2 in paired end data.
    """
    library = int(library)
    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample) &
        (cohort_data.library_index == library) &
        (cohort_data.read_in_pair == read_in_pair)
    ]
    return(sample_line.path.to_list()[0])

def get_all_samples(cohort=None):
    """Retrieves all samples to be processed.
    Does so by calling get_cohort_data, and therefore filters out excluded_cases.

    Keyword arguments:
        cohort -- Name of a cohort, OPTIONAL. If not specified, returns all samples
                  across all cohorts.
    """
    all_samples = pd.concat([
        get_cohort_data(cohort_name).assign(cohort_name = cohort_name)
        for cohort_name
        in config['data']['cohorts']
        if config['data']['cohorts'][cohort_name]['active']
        ])

    if cohort is None:
        return(all_samples)
    else:
        return(all_samples[all_samples.cohort_name == cohort])

def get_all_samples_list():
    """Returns all samples in list format.
    By calling get_all_samples()
    """
    return(get_all_samples().sample_name.unique().tolist())

def get_all_samples_with_cohorts():
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples()[['cohort_name', 'sample_name']]
    return(zip(
        samples.drop_duplicates().cohort_name.tolist(),
        samples.drop_duplicates().sample_name.tolist()
    ))

def clean(command):
    """Cleans a snakemake command string by replacing whitespace with a single space.
    Useful for stripping down multiline shell commands so that they read more
    cleanly on snakemake's printed output.
    """
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return(command.strip())

def get_cohort_signature_pairs():
    """Iterator returning every signature for every active cohort."""
    for active_cohort in get_active_cohorts():
        for signature in get_cohort_config(active_cohort)['signatures']:
            yield (active_cohort, signature)

# ------------------------------ #
#  Beginning of Snakemake Rules  #
# ------------------------------ #

rule all:
    input:
        #expand(
        #    path_to_data + '/{cohort}/results/results_matrix/{cohort}_{colname}{ext}',
        #    cohort = get_active_cohorts(),
        #    colname = ['coverage', 'posterior'],
        #    ext = ['.csv']
        #),
        expand([
            path_to_data + '/inspire/results/extract_signature_matrix/{}-{}-{{valuetype}}.parquet'.format(cohort, signature)
                for cohort, signature in get_cohort_signature_pairs()
            ],
            valuetype = ['coverage', 'posterior', 'medestrand'])
        #[path_to_data + '/{cohort}/results/bam_to_wig/{sample}_counts.wig'.format(
        #        cohort=v[0],
        #        sample=v[1]
        #    ) for v in get_all_samples_with_cohorts()
        #],

# ------------------------------------------ #
#  Install conda dependencies automatically  #
# ------------------------------------------ #
# This section is only necessary to support air-gapped environments,
# where the cluster nodes have no direct access to internet.
# In this setting, you must first load all necessary packages from
# a build node with internet access. This can be done simply
# by running 
# bash update_conda.sh
# from a node with internet access, which loads install_dependencies
# rule below to install all of the necessary conda dependencies
# in advance. Any time the git repo is pulled or if you modify
# any of the conda environments in conda_env, this needs to be
# re-run.

rule install_dependencies:
    input:
        'conda_env/biopython_installed',
        'conda_env/cfmedip_r_installed',
        'conda_env/fastqc_installed',
        'conda_env/qualimap_installed',
        'conda_env/samtools_installed',
        'conda_env/trimgalore_installed'
    output:
        'conda_env/dependencies'
    shell:
        'touch {output}'

rule biopython:
    output:
        'conda_env/biopython_installed'
    conda: 'conda_env/biopython.yml'
    shell:
        'touch {output}'

rule cfmedip_r:
    output:
        'conda_env/cfmedip_r_installed'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'touch {output}'

rule fastqc:
    output:
        'conda_env/fastqc_installed'
    conda: 'conda_env/fastqc.yml'
    shell:
        'touch {output}'

rule qualimap:
    output:
        'conda_env/qualimap_installed'
    conda: 'conda_env/qualimap.yml'
    shell:
        'touch {output}'

rule samtools:
    output:
        'conda_env/samtools_installed'
    conda: 'conda_env/samtools.yml'
    shell:
        'touch {output}'

rule trimgalore:
    output:
        'conda_env/trimgalore_installed'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'touch {output}'

# -------------------------- #
#  Pre-process input FASTQs  #
# -------------------------- #

rule gunzip_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(wildcards.read)),
    output:
        temp(path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'gunzip -dc {input} > {output}'

# QC of input FASTQs using FASTQC
rule fastqc_fastq:
    input:
        path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq'
    output:
        html=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.html',
        zipfile=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Extract Barcodes using ConsensusCruncher
# Pulls the path to extract_barcodes.py from config > paths > dependencies > extract_barcodes_path
rule extract_barcodes:
    input:
        R1_qc = path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
        R2_qc = path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html',
        R1 = path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R1.fastq',
        R2 = path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R2.fastq'
    output:
        R1 = temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R1.fastq'),
        R2 = temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R2.fastq')
    params:
        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0],
        barcodes = lambda wildcards: get_cohort_config(wildcards.cohort)['barcodes']
    resources: cpus=1, mem_mb=16000, time_min='5-00:00:00'
    conda: 'conda_env/biopython.yml'
    shell:
        clean(r'''
        python {extract_barcodes}
            --read1 {{input.R1}}
            --read2 {{input.R2}}
            --outfile {{params.outprefix}}
            {{params.barcodes}}
        '''.format(extract_barcodes = config['paths']['dependencies']['extract_barcodes_path']))

# Trims FASTQ using trimgalore to remove barcode sequences
# By default, trims 10 base pairs from the 5' end, which seems to be correct for OICR cfMeDIP-seq output.
# This can be configured in the config.yml under data > cohorts > settings > trimgalore.
rule trim_fastq:
    input:
        R1 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R1.fastq',
        R2 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R2.fastq'
    output:
        trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq'),
        trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq'),
        report_1 = path_to_data + '/{cohort}/results/qc/{sample}_lib{lib}_extract_barcode_R1.fastq_trimming_report.txt',
        report_2 = path_to_data + '/{cohort}/results/qc/{sample}_lib{lib}_extract_barcode_R2.fastq_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=4, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 4 --dont_gzip --paired {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} {input.R2} && cp ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq/{wildcards.sample}_lib{wildcards.lib}_extract_barcode_R*.fastq_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc'

# --------------------------- #
#  Align FASTQs to reference  #
# --------------------------- #

def get_read_group_from_fastq(fastq_file, sample_name):
    """Extracts the read group data from FASTQs automatically.
    This information is used in the shell script of rule bwa_mem.
    """
    with gzip.open(fastq_file, 'rt') as fastq:
        header = next(fastq)
        (instrument, run_number, flowcell, lane, tile, xpos, ypos) = header.split(' ')[0].split(':')
        lib_value = header.strip().split(' ')[1].split(':')[3]
        rg_line = r"@RG\tID:{flowcell}_{lane}\tSM:{sample}\tPL:Illumina\tPU:.\tLB:{lib_value}".format(
            flowcell = flowcell,
            lane = lane,
            sample = sample_name,
            lib_value = lib_value
        )
        return(rg_line)

# Run BWA mem on FASTQs after extracting barcodes.
rule bwa_mem:
    input:
        path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq',
        path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq',
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam')
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r"""
        bwa mem -M -t4
        -R'{params.read_group}' 
        {params.bwa_index}
        {input} > {output}""")

# Converts SAM to BAM with sorting
rule sam_to_sorted_bam:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam'),
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        # Try it without fixmate for replication purposes
        #"samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"
        clean(r'''
        samtools view -buS -f 2 -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam}
        ''')

def get_libraries_of_sample(sample):
    """Returns all library indices of a sample based on samplesheet."""
    filtered_table = get_all_samples()[get_all_samples().sample_name == sample]
    return(list(set(filtered_table.library_index.to_list())))

# If there are multiple libraries for a given sample, as specified in samplesheet,
# these libraries are automatically merged at this step into a single unified BAM.
rule merge_bam:
    input:
        lambda wildcards: expand(
                path_to_data + '/{{cohort}}/tmp/bwa_mem/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        temp(path_to_data + '/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output} {input} && samtools index {output}'

# Bam markdup and create index. 
# This step finalizes the definitive BAM file.
rule bam_markdup:
    input:
        path_to_data + '/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam'
    output:
        bam = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
        index = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

# ----------------- #
#  QC of BAM files  #
# ----------------- #

# Run FASTQC on final BAM files.
rule fastqc_bam:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        html = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.html',
        zipfile = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.zip',
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Run Qualimap QC metrics on final BAM files.
rule qualimap_bam:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        html = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.html',

# Run Flagstat to get basic stats on final BAM files.
rule bam_flagstat:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup.bam.flagstat',
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools flagstat {input} > {output}'

# Unified QC rule which runs all of the above QCs
rule bam_qc:
    input:
        fastqc = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.html',
        flagstat = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup.bam.flagstat',
    output:
        path_to_data + '/{cohort}/results/qc/{sample}_qc_complete',
    resources: cpus=1, mem_mb=1000, time_min='00:01:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'touch {output}'

# ------------------------------------------ #
#  Compute BAM bin stats including coverage  #
# ------------------------------------------ #

# Generate bin stats for each chromosome individually.
rule bam_bin_stats:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        binstat = temp(path_to_data + '/{cohort}/tmp/bam_bin_stats/bin_stats_{sample}_{species}_{chrom}.tsv'),
        filtered = path_to_data + '/{cohort}/results/bam_bin_stats/removed_bins_{sample}_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome'][wildcards.species],
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        clean('''
        Rscript src/R/bin_stats.R
            -b {input}
            -g {params.bsgenome}
            -c {wildcards.chrom}
            -o {output.binstat}
            --filtered {output.filtered}
        ''')

# Preload chromosome values for below rule
chromosomes = {}
for c in config['data']['cohorts']:
    if config['data']['cohorts'][c]['active']:
        chr_data = get_cohort_config(c)['chromosomes']
        chromosome_tuples = [(species, chrom) for species in chr_data for chrom in chr_data[species].split(',')]
        chromosomes[c] = chromosome_tuples

# Merge the bin stats across all chromosomes.
rule merge_bin_stats:
    input:
        lambda wildcards: [path_to_data + '/{{cohort}}/tmp/bam_bin_stats/bin_stats_{{sample}}_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in chromosomes[wildcards.cohort]],
    output:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    params:
        paths = lambda wildcards, input: ','.join(input)
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/row_bind_tables.R -p "{params.paths}" -o {output} --in-tsv --out-feather --omit-paths'

# ----------------------------------------------------------------------------- #
#  Fit cfMeDIP-seq coverage stats to infer absolute methylation using MedReMix  #
# ----------------------------------------------------------------------------- #

# This is the currently used negative binomial GLM approach to fitting
rule cfmedip_nbglm:
    input:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        fit=path_to_data + '/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        model=path_to_data + '/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm_model.Rds',
    resources: cpus=1, time_min='5-00:00:00', mem_mb=lambda wildcards, attempt: 16000 if attempt == 1 else 30000
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'

# Below is the full bayesian approach for fitting, which remains under development
rule fit_bin_stats:
    input:
        path_to_data + '/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        path_to_data + '/{cohort}/results/fit_bin_stats/{sample}_fit_{method}.tsv'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/fit_cpg_bins.R -i {input} -o {output} --method {wildcards.method}'

# --------------------------------------------- #
#  MeDEStrand Method of methylation correction  #
# --------------------------------------------- #

rule medestrand:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        path_to_data + '/{cohort}/results/medestrand/{sample}_medestrand_output.feather',
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    params: bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome']['human']
    shell:
        'Rscript src/R/run_medestrand.R -b {{input}} -o {{output}} -m {} --bsgenome {{params.bsgenome}}'.format(
            config['paths']['dependencies']['medestrand_path']
        )

# ----------------------------------------------------------- #
#  Assemble results overlapping specific regions of interest  #
# ----------------------------------------------------------- #
# Signatures define regions of interest from which we
# would like to capture quantitative methylation information.
#
# The results are similar to those from the matrices/binsets,
# except that signatures don't have to match bins in the data,
# and the scripts here will automatically find all bins
# overlapping the signature regions and output data from them.

# Each type of matrix is composed from specific files
# and specific column within each file category.
# For example, a matrix of posteriors pulls data from
# column 8 of the output files of cfmedip_nbglm.

data_sources = {
    'coverage': {
        'path': path_to_data + '/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        'filetype': 'feather',
        'colname': 'coverage',
        'clip': '_fit_nbglm.feather'
    },
    'posterior': {
        'path': path_to_data + '/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        'filetype': 'feather',
        'colname': 'methylated_posterior',
        'clip': '_fit_nbglm.feather'
    },
    'medestrand': {
        'path': path_to_data + '/{cohort}/results/medestrand/{sample}_medestrand_output.feather',
        'filetype': 'feather',
        'colname': 'binMethyl',
        'clip': '_medestrand_output.feather'
    }
}

def get_matrix_inputs(wildcards):
    """Get paths of all files used to assemble a specific matrix.
    """
    samples = get_cohort_data(wildcards.cohort).sample_name.unique().tolist()
    paths = expand(data_sources[wildcards.valuetype]['path'], cohort = wildcards.cohort, sample = samples)
    return(paths)

# Assembles paths of files for assembling a matrix and
# write them to a file, one path per line
rule signature_matrix_paths:
    input:
        get_matrix_inputs
    output:
        temp(path_to_data + '/{cohort}/tmp/signature_matrix_paths/matrix.{valuetype}.paths.txt')
    run:
        with open(output[0], 'w') as out:
            out.write('\n'.join(input))

# Using the files assembled in signature_matrix_paths,
# and signature regions specified in a file pointed to in
# config > data > cohorts > signatures or the defaults
# specified in config > data > defaults > signatures.
rule extract_signature_matrix:
    input:
        signature = lambda wildcards: get_cohort_config(wildcards.cohort)['signatures'][wildcards.signature],
        paths = path_to_data + '/{cohort}/tmp/signature_matrix_paths/matrix.{valuetype}.paths.txt'
    output:
        path_to_data + '/{cohort}/results/extract_signature_matrix/{cohort}-{signature}-{valuetype}.parquet'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    params:
        colname = lambda wildcards: data_sources[wildcards.valuetype]['colname'],
        filetype = lambda wildcards: data_sources[wildcards.valuetype]['filetype'],
        clip = lambda wildcards: data_sources[wildcards.valuetype]['clip']
    shell:
        clean("""
        Rscript src/R/extract_signature_matrix.R 
            -s {input.signature} 
            -p {input.paths} 
            -c {params.colname} 
            -f {params.filetype}
            -o {output} 
            --clip '{params.clip}'
        """)

# -------------- #
#  Legacy rules  #
# -------------- #

#include: 'legacy.smk'
# This line includes some rules from legacy code, which are
# gradually being retired but still accessible if needed.


